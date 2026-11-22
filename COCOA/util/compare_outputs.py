#!/usr/bin/env python3
"""Compare two COCOA output files produced with the same configuration and seed.

Regression check for code changes: branches present in both files must agree
event by event (floats exactly, NaN == NaN). Truth-graph edges are compared as
sets. If the reference file uses output schema version 1 (primary_photon_*
branches) and the candidate uses version 2, the former photon branches are
also compared with the corresponding particle_* columns, and the former
conversion legs with the primary-photon legs of conv_el_*.

Usage:
    python3 util/compare_outputs.py REFERENCE.root CANDIDATE.root [--tree Out_Tree] [--ignore BRANCH ...]

Exits with status 1 if a compared quantity differs.
"""

from __future__ import annotations

import argparse
import sys

import awkward as ak
import numpy as np
import uproot

# Schema 1 photon branch -> schema 2 particle_* column (indexed by primary_photon_particle_idx).
LEGACY_PHOTON_MAP = {
    "energy": "particle_e",
    "eta": "particle_eta",
    "phi": "particle_phi",
    "prod_x": "particle_prod_x",
    "prod_y": "particle_prod_y",
    "prod_z": "particle_prod_z",
    "end_x": "particle_end_x",
    "end_y": "particle_end_y",
    "end_z": "particle_end_z",
    "path_length": "particle_track_length",
    "termination_process_type": "particle_end_process_type",
    "termination_process_subtype": "particle_end_process_subtype",
    "track_id": "particle_g4_track_id",
}


def equal(x: ak.Array, y: ak.Array) -> bool:
    if len(x) != len(y):
        return False
    if x.ndim > 1 and not ak.all(ak.num(x, axis=1) == ak.num(y, axis=1)):
        return False
    same = x == y
    if np.issubdtype(np.asarray(ak.flatten(x, axis=None)).dtype, np.floating):
        same = same | (np.isnan(x) & np.isnan(y))
    return bool(ak.all(same))


def edge_sets(arrays: ak.Array) -> list[set]:
    return [set(zip(s, e)) for s, e in zip(ak.to_list(arrays["node_edge_start"]), ak.to_list(arrays["node_edge_end"]))]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference")
    parser.add_argument("candidate")
    parser.add_argument("--tree", default="Out_Tree")
    parser.add_argument("--ignore", nargs="*", default=[], help="branches not to compare")
    args = parser.parse_args()

    ref_tree = uproot.open(args.reference)[args.tree]
    new_tree = uproot.open(args.candidate)[args.tree]
    ref = ref_tree.arrays(library="ak")
    new = new_tree.arrays(library="ak")
    failures = []

    if ref_tree.num_entries != new_tree.num_entries:
        print(f"Different number of events: {ref_tree.num_entries} vs {new_tree.num_entries}")
        return 1

    common = sorted(set(ref_tree.keys()) & set(new_tree.keys()) - set(args.ignore))
    edges = [k for k in common if k.startswith("node_edge_")]
    compared = [k for k in common if k not in edges]
    differing = [k for k in compared if not equal(ref[k], new[k])]
    print(f"{len(compared)} common branches compared, {len(compared) - len(differing)} identical")
    for k in differing:
        print(f"  differs: {k}")
    failures += differing

    if len(edges) == 2:
        same_edges = edge_sets(ref) == edge_sets(new)
        n_ref = int(ak.sum(ak.num(ref["node_edge_start"])))
        n_new = int(ak.sum(ak.num(new["node_edge_start"])))
        print(f"truth-graph edge sets identical: {same_edges} ({n_ref} vs {n_new} stored edges)")
        if not same_edges:
            failures.append("node_edge_*")

    only_ref = sorted(set(ref_tree.keys()) - set(new_tree.keys()))
    only_new = sorted(set(new_tree.keys()) - set(ref_tree.keys()))
    print(f"{len(only_ref)} branches only in the reference, {len(only_new)} only in the candidate")

    if "primary_photon_particle_idx" in ref_tree.keys() and "particle_end_x" in new_tree.keys():
        print("Schema 1 -> 2 mapping")
        idx = ref["primary_photon_particle_idx"]
        for old, column in LEGACY_PHOTON_MAP.items():
            ok = equal(ref["primary_photon_" + old], new[column][idx])
            print(f"  primary_photon_{old} == {column}: {ok}")
            if not ok:
                failures.append(f"primary_photon_{old}")
        converted = (new["particle_end_process_type"][idx] == 2) & (new["particle_end_process_subtype"][idx] == 14)
        ok = equal(ref["primary_photon_converted"], ak.values_astype(converted, np.int32))
        print(f"  primary_photon_converted == ended by conversion: {ok}")
        if not ok:
            failures.append("primary_photon_converted")

        primary_leg = new["conv_vtx_in_particle_idx"][new["conv_el_vtx_idx"]] >= 0
        for c in ("q", "px", "py", "pz"):
            ok = equal(ref["conv_el_" + c], new["conv_el_" + c][primary_leg])
            print(f"  conv_el_{c} == primary-photon legs: {ok}")
            if not ok:
                failures.append(f"conv_el_{c} (legacy)")
        ok = equal(ref["conv_el_primary_photon_idx"],
                   new["conv_vtx_in_particle_idx"][new["conv_el_vtx_idx"]][primary_leg])
        print(f"  conv_el_primary_photon_idx == conv_vtx_in_particle_idx: {ok}")
        for c in "xyz":
            ok_c = equal(ref[f"conv_el_prod_{c}"], new[f"conv_vtx_{c}"][new["conv_el_vtx_idx"]][primary_leg])
            print(f"  conv_el_prod_{c} == conv_vtx_{c}: {ok_c}")
            ok = ok and ok_c
        if not ok:
            failures.append("conversion vertex (legacy)")

    print("\nRESULT:", "identical" if not failures else f"{len(failures)} differences")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
