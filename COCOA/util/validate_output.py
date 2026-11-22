#!/usr/bin/env python3
"""Check a COCOA output file against output schema version 2.

Checks the structure (table lengths, links, run-level trees) and the physics
consistency of the truth record described in docs/output_schema.md, and prints
a short summary of its content. Exits with status 1 if any check fails.

Usage:
    python3 util/validate_output.py OUTPUT.root [--tree Out_Tree]

Requires uproot, awkward and numpy (all available in the LCG views).
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter

import awkward as ak
import numpy as np
import uproot

TABLES = ("particle", "conv_vtx", "conv_el", "nucl_vtx", "nucl_dau", "node")


class Checker:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.passed = 0

    def check(self, condition: bool, message: str) -> None:
        if condition:
            self.passed += 1
        else:
            self.failures.append(message)
            print(f"  FAIL: {message}")


def table_branches(keys: list[str], table: str) -> list[str]:
    branches = [k for k in keys if k.startswith(table + "_")]
    if table == "particle":
        # particle_* also matches no other table prefix.
        return branches
    return branches


def valid_index(index: ak.Array, target_length: ak.Array, allow_minus_one: bool = True) -> bool:
    lower = -1 if allow_minus_one else 0
    flat_ok = ak.all((index >= lower) & (index < target_length), axis=None)
    return bool(flat_ok) if len(index) else True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file")
    parser.add_argument("--tree", default="Out_Tree")
    args = parser.parse_args()

    c = Checker()
    f = uproot.open(args.file)
    tree = f[args.tree]
    keys = list(tree.keys())
    a = tree.arrays(library="ak")
    n_events = tree.num_entries
    print(f"{args.file}: {n_events} events, {len(keys)} branches")

    # ---- Run-level trees ----------------------------------------------------
    print("Run-level trees")
    for name in ("RunInfo", "Materials", "Volumes"):
        c.check(name in [k.split(";")[0] for k in f.keys()], f"tree {name} present")
    run_info = f["RunInfo"].arrays(library="np")
    c.check(int(run_info["schema_version"][0]) == 2, "RunInfo.schema_version == 2")
    c.check(int(run_info["n_events_processed"][0]) == n_events, "RunInfo.n_events_processed equals the number of entries")
    c.check(len(str(run_info["geometry_fingerprint"][0])) == 32, "RunInfo.geometry_fingerprint is an MD5 hex digest")
    c.check(str(run_info["cocoa_git_commit"][0]) != "unknown", "RunInfo.cocoa_git_commit known")
    materials = f["Materials"].arrays(library="np")
    volumes = f["Volumes"].arrays(library="np")
    c.check(np.array_equal(materials["material_idx"], np.arange(len(materials["material_idx"]))), "Materials.material_idx is 0..N-1")
    c.check(np.array_equal(volumes["volume_idx"], np.arange(len(volumes["volume_idx"]))), "Volumes.volume_idx is 0..N-1")
    n_materials = len(materials["material_idx"])
    n_volumes = len(volumes["volume_idx"])
    r_max = float(run_info["tracker_region_r_max_mm"][0])
    z_max = float(run_info["tracker_region_abs_z_max_mm"][0])
    print(f"  commit {run_info['cocoa_git_describe'][0]} (dirty={run_info['cocoa_git_dirty'][0]}), "
          f"seed {run_info['run_seed'][0]}, tracker region r < {r_max:.0f} mm, |z| < {z_max:.0f} mm")

    # ---- Event-level and table structure --------------------------------------
    print("Structure")
    c.check("event_id" in keys and "run_seed" in keys, "event_id and run_seed present")
    if "event_id" in keys:
        c.check(np.array_equal(ak.to_numpy(a["event_id"]), np.arange(n_events)), "event_id is 0..N-1")
    removed = [k for k in keys if k.startswith(("primary_photon_", "nucl_int_")) or k in ("conv_el_primary_photon_idx", "conv_el_prod_x")]
    c.check(not removed, f"schema-1 branches removed (found: {removed[:3]})")
    length = {}
    for table in TABLES:
        branches = table_branches(keys, table)
        if table == "node":
            branches = [b for b in branches if not b.startswith("node_edge_")]
        if not branches:
            c.check(table == "node", f"table {table}_* present")
            continue
        counts = {b: ak.to_numpy(ak.num(a[b], axis=1)) for b in branches}
        reference = next(iter(counts.values()))
        mismatched = [b for b, n in counts.items() if not np.array_equal(n, reference)]
        c.check(not mismatched, f"{table}_* branches have equal lengths (mismatch: {mismatched[:3]})")
        length[table] = ak.num(a[branches[0]], axis=1)

    # ---- Links ------------------------------------------------------------------
    print("Links")
    n_particle = length["particle"]
    c.check(valid_index(a["particle_track_idx"], ak.num(a["track_pdgid"], axis=1)), "particle_track_idx -> track_*")
    c.check(valid_index(a["particle_end_material_idx"], n_materials), "particle_end_material_idx -> Materials")
    c.check(valid_index(a["particle_end_volume_idx"], n_volumes), "particle_end_volume_idx -> Volumes")
    if "node" in length:
        c.check(valid_index(a["particle_truth_node_idx"], length["node"]), "particle_truth_node_idx -> node_*")
        c.check(valid_index(a["node_final_state_idx"], n_particle), "node_final_state_idx -> particle_*")
        fs = a["node_final_state_idx"]
        node_rows = ak.local_index(fs)
        back = a["particle_truth_node_idx"][ak.where(fs >= 0, fs, 0)]
        c.check(bool(ak.all((fs < 0) | (back == node_rows))), "node_final_state_idx and particle_truth_node_idx are inverse")
        pdg_match = a["node_pdg_id"][fs >= 0] == a["particle_pdgid"][fs[fs >= 0]]
        c.check(bool(ak.all(pdg_match)), "linked node and particle have the same PDG ID")
        edges = list(zip(ak.to_list(a["node_edge_start"]), ak.to_list(a["node_edge_end"])))
        duplicates = sum(len(s) - len(set(zip(s, e))) for s, e in edges)
        c.check(duplicates == 0, f"truth-graph edges unique ({duplicates} duplicates)")

    for vtx, out, legs in (("conv_vtx", "conv_el", True), ("nucl_vtx", "nucl_dau", False)):
        c.check(valid_index(a[f"{vtx}_in_particle_idx"], n_particle), f"{vtx}_in_particle_idx -> particle_*")
        c.check(valid_index(a[f"{vtx}_in_ancestor_idx"], n_particle), f"{vtx}_in_ancestor_idx -> particle_*")
        c.check(valid_index(a[f"{vtx}_material_idx"], n_materials, False), f"{vtx}_material_idx -> Materials")
        c.check(valid_index(a[f"{vtx}_volume_idx"], n_volumes, False), f"{vtx}_volume_idx -> Volumes")
        c.check(valid_index(a[f"{out}_vtx_idx"], length[vtx], False), f"{out}_vtx_idx -> {vtx}_*")
        primary = a[f"{vtx}_in_particle_idx"]
        is_primary = primary >= 0
        c.check(bool(ak.all(a["particle_pdgid"][primary[is_primary]] == a[f"{vtx}_in_pdgid"][is_primary])),
                f"{vtx}_in_pdgid equals the PDG ID of the linked primary")
        c.check(bool(ak.all(a[f"{vtx}_in_ancestor_idx"][is_primary] == primary[is_primary])),
                f"{vtx}_in_ancestor_idx of a primary is the primary itself")
        if legs:
            for leg, sign in (("ep", 1), ("em", -1)):
                idx = a[f"conv_vtx_{leg}_idx"]
                c.check(valid_index(idx, length["conv_el"], False), f"conv_vtx_{leg}_idx -> conv_el_* (always set)")
                c.check(bool(ak.all(a["conv_el_vtx_idx"][idx] == ak.local_index(idx))), f"conv_el_vtx_idx[conv_vtx_{leg}_idx] points back")
                c.check(bool(ak.all(np.sign(a["conv_el_q"][idx]) == sign)), f"conv_vtx_{leg}_idx has charge {sign:+d}")

    # ---- Physics consistency ---------------------------------------------------
    print("Physics consistency")
    for vtx in ("conv_vtx", "nucl_vtx"):
        r = np.hypot(a[f"{vtx}_x"], a[f"{vtx}_y"])
        c.check(bool(ak.all((r <= r_max + 1e-3) & (abs(a[f"{vtx}_z"]) <= z_max + 1e-3))), f"{vtx} inside the tracker region")
        c.check(bool(ak.all(a[f"{vtx}_in_x0"] >= 0)) and bool(ak.all(a[f"{vtx}_in_lambda"] >= 0)), f"{vtx}_in_x0/lambda >= 0")
    c.check(bool(ak.all(a["conv_vtx_in_pdgid"] == 22)), "conv_vtx_in_pdgid == 22")
    secondary = a["conv_vtx_in_particle_idx"] < 0
    e_min = float(run_info["truth_record_conversion_secondary_photon_e_min_mev"][0])
    c.check(bool(ak.all(a["conv_vtx_in_e"][secondary] >= e_min * (1 - 1e-6))), "secondary conversions above the recording threshold")
    pair_e = a["conv_el_e"][a["conv_vtx_ep_idx"]] + a["conv_el_e"][a["conv_vtx_em_idx"]]
    rel = ak.flatten(abs(pair_e - a["conv_vtx_in_e"]) / a["conv_vtx_in_e"])
    c.check(bool(ak.all(rel < 1e-3)), f"e+ + e- energy equals the photon energy at the conversion (max rel. diff {ak.max(rel) if len(rel) else 0:.2e})")

    # Primary photons that converted inside the tracker: particle_end_* is the vertex.
    p = a["conv_vtx_in_particle_idx"]
    prim = p >= 0
    for axis in "xyz":
        diff = ak.flatten(abs(a[f"particle_end_{axis}"][p[prim]] - a[f"conv_vtx_{axis}"][prim]))
        c.check(bool(ak.all(diff < 1e-3)), f"particle_end_{axis} equals conv_vtx_{axis} for primary photons")
    diff = ak.flatten(abs(a["particle_end_e"][p[prim]] - a["conv_vtx_in_e"][prim]))
    c.check(bool(ak.all(diff < 1e-3)), "particle_end_e equals conv_vtx_in_e for primary photons")
    converted = (a["particle_pdgid"] == 22) & (a["particle_end_process_type"] == 2) & (a["particle_end_process_subtype"] == 14) & (a["particle_end_in_tracker"] == 1)
    c.check(int(ak.sum(converted)) == int(ak.sum(prim)), "every primary photon converting in the tracker has one conv_vtx row")
    photon = a["particle_pdgid"] == 22
    c.check(bool(ak.all(a["particle_end_e"][photon] <= a["particle_e"][photon] * (1 + 1e-6))), "photon energy at its end <= production energy")
    straight = photon & (a["particle_n_scatters_tracker"] == 0)
    c.check(bool(ak.all(abs(a["particle_end_e"][straight] - a["particle_e"][straight]) <= 1e-3 * a["particle_e"][straight])),
            "photons without scatters keep their energy")
    c.check(bool(ak.all(a["particle_x0_tracker"] >= 0)), "particle_x0_tracker >= 0")

    for out in ("conv_el", "nucl_dau"):
        reached = a[f"{out}_calo_reached"] == 1
        r = np.hypot(a[f"{out}_calo_x"][reached], a[f"{out}_calo_y"][reached])
        on_boundary = (abs(r - r_max) < 1e-2) | (abs(abs(a[f"{out}_calo_z"][reached]) - z_max) < 1e-2)
        c.check(bool(ak.all(on_boundary)), f"{out}_calo_* lies on the tracker-region boundary")
        c.check(bool(ak.all(np.isnan(a[f"{out}_calo_e"][~reached]))), f"{out}_calo_e is NaN when not reached")
        c.check(bool(ak.all(a[f"{out}_calo_e"][reached] <= a[f"{out}_e"][reached] * (1 + 1e-6))), f"{out}_calo_e <= production energy")
    n_stored = ak.Array([np.bincount(ak.to_numpy(v), minlength=n) for v, n in zip(a["nucl_dau_vtx_idx"], length["nucl_vtx"])])
    c.check(bool(ak.all(n_stored <= a["nucl_vtx_n_dau_total"])), "stored daughters <= nucl_vtx_n_dau_total")
    hadron_ekin_min = float(run_info["truth_record_nuclear_secondary_hadron_ekin_min_mev"][0])
    del hadron_ekin_min  # kinetic energy is not stored; threshold checked in C++ only

    # ---- Summary ---------------------------------------------------------------
    def total(branch):
        return int(ak.sum(ak.num(a[branch], axis=1)))

    print("Content")
    print(f"  particles {total('particle_pdgid')}, of which photons {int(ak.sum(photon))}")
    print(f"  conversions {total('conv_vtx_x')} (primary photons {int(ak.sum(prim))}), legs {total('conv_el_q')}")
    creators = Counter(zip(ak.to_list(ak.flatten(a["conv_vtx_in_creator_type"][secondary])),
                           ak.to_list(ak.flatten(a["conv_vtx_in_creator_subtype"][secondary]))))
    print(f"  secondary-photon conversions by creator (type, subtype): {dict(creators.most_common(6))}")
    nprim = int(ak.sum(a["nucl_vtx_in_particle_idx"] >= 0))
    print(f"  nuclear interactions {total('nucl_vtx_x')} (primary hadrons {nprim}), daughters {total('nucl_dau_q')}")
    print(f"  legs reaching the calorimeter {int(ak.sum(a['conv_el_calo_reached']))} / {total('conv_el_q')}")

    print(f"\n{c.passed} checks passed, {len(c.failures)} failed")
    return 1 if c.failures else 0


if __name__ == "__main__":
    sys.exit(main())
