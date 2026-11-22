# COCOA output schema (version 2)

This document describes the truth-level content of the COCOA ROOT output and the
run-level metadata written alongside it. It is the reference for analysis code
reading COCOA files. Reconstruction-level tables (`cell_*`, `track_*`, `topo_*`,
jets, graph edges, particle flow, superclusters) are unchanged by schema version 2
and are only listed briefly.

## Conventions

- **File layout.** One ROOT file per run. It contains the event tree `Out_Tree`
  (`Low_Tree` and `High_Tree` when `Use_high_granularity` is on) and the run-level
  trees `RunInfo`, `Materials` and `Volumes`.
- **One entry per event.** `Out_Tree` has exactly one entry per simulated event,
  including events in which nothing was recorded. The number of entries is
  therefore the number of simulated events.
- **Tables.** A branch prefix defines a table, for example `particle_*`. All
  branches of a table are `std::vector`s with the same length in a given event;
  element *i* of every branch describes row *i* of that table.
- **Links.** A column named `<table>_<target>_idx` holds a row index into the
  table `<target>` *of the same event*. `-1` means "no such row". Geant4 track IDs
  are stored for information only (`*_g4_track_id`) and are never used as links.
- **Units.** Energies and momenta in MeV, lengths and positions in mm, angles in
  rad, material amounts in units of the radiation length X0 or the nuclear
  interaction length λ_I (dimensionless). Energies `e` are total energies.
- **Floating point.** Floats are stored as `float32`. `NaN` appears only where
  documented as "undefined".
- **Process codes.** Processes are stored as raw Geant4 process type and subtype
  integers. The most relevant ones are listed in [Process codes](#process-codes).
- **Tracker region.** "Inside the tracker" means inside the Geant4 region
  `TrackerRegion`, i.e. the inner-detector envelope `Inner_LV`: a cylinder with
  radius `Inner_calorimeter_layer` (1500 mm by default) and half-length
  `Inner_calorimeter_layer / tan(θ_min)`, where θ_min corresponds to
  `Max_eta_of_barrel_region` (3194 mm by default). Its outer surface is the front
  face of the calorimeter. The extent used in a run is stored in `RunInfo`.

## Event-level branches

| Branch | Type | Description |
|---|---|---|
| `event_id` | int | Geant4 event ID within the run (0, 1, ...). |
| `run_seed` | int64 | Random seed of the run (`--seed`, or the time-based default). Together with `event_id` this identifies an event uniquely across merged files. |

## `particle_*`: Geant4 primary particles

One row per Geant4 primary track, i.e. per generator particle passed to Geant4
after the fiducial cuts (`Fiducial_cuts`). Only primaries appear in this table,
independent of any other configuration option.

| Branch | Type | Description |
|---|---|---|
| `particle_pdgid` | int | PDG ID. |
| `particle_q` | float | Charge in units of e. |
| `particle_g4_track_id` | int | Geant4 track ID (information only). |
| `particle_pt`, `particle_eta`, `particle_phi`, `particle_e` | float | Generator momentum at production: pT, η, φ, total energy. |
| `particle_prod_x`, `_y`, `_z` | float | Production vertex. |
| `particle_isIso` | int | 1 if the scalar pT sum of other charged primaries within ΔR < min(10 GeV / pT, 0.2) is below 6 % of the particle's pT; 0 for neutral particles. |
| `particle_truth_node_idx` | int | Row in `node_*` of this particle; -1 without truth graph. |
| `particle_track_idx` | int | Row in `track_*` of the emulated track; -1 if none. |
| `particle_eta_extrap_calo`, `particle_phi_extrap_calo` | float | η, φ of the highest-energy same-PDG descendant step in the calorimeter. |
| `particle_eta_extrap_its`, `particle_phi_extrap_its` | float | Same in the outermost inner-detector layer. |
| `particle_dep_energy` | float | Energy deposited in topo-clustered cells. |
| `particle_end_x`, `_y`, `_z` | float | Position where the track ended (end of its last step). |
| `particle_end_e` | float | Total energy at the start of the last step. For a photon this is the energy at its final interaction (e.g. the conversion). |
| `particle_end_px`, `_py`, `_pz` | float | Momentum at the start of the last step. |
| `particle_end_process_type`, `particle_end_process_subtype` | int | Process that ended the track. Transportation (1, 91) means the particle left the world. |
| `particle_end_in_tracker` | int | 1 if the track ended inside the tracker region. |
| `particle_end_material_idx` | int | Row in `Materials` of the material where the track ended. |
| `particle_end_volume_idx` | int | Row in `Volumes` of the logical volume where the track ended. |
| `particle_track_length` | float | Total track length. |
| `particle_x0_tracker`, `particle_lambda_tracker` | float | Material traversed inside the tracker region along the true path, Σ step / X0 and Σ step / λ_I, up to the end of the track or its exit from the region. |
| `particle_n_scatters_tracker` | int | Number of Compton, Rayleigh and hadronic elastic interactions inside the tracker region. 0 for a photon means it travelled on a straight line until it ended or left the tracker. |

## `conv_vtx_*`: photon conversions

One row per γ → e⁺e⁻ conversion (EM subtype 14) whose interaction point lies
inside the tracker region and whose photon is either

- a Geant4 primary (always recorded), or
- a secondary photon with energy ≥ `Truth_record.Conversion_secondary_photon_e_min_mev`
  (recorded only if `Truth_record.Record_secondary_conversions` is true).

| Branch | Type | Description |
|---|---|---|
| `conv_vtx_x`, `_y`, `_z` | float | Conversion point. |
| `conv_vtx_in_*` | | Incoming photon, see [Incoming-particle columns](#incoming-particle-columns). |
| `conv_vtx_material_idx`, `conv_vtx_volume_idx` | int | Material and logical volume at the conversion point. |
| `conv_vtx_ep_idx`, `conv_vtx_em_idx` | int | Rows in `conv_el_*` of the e⁺ and the e⁻. |

## `conv_el_*`: conversion legs

One row per e⁺ or e⁻ created by a conversion recorded in `conv_vtx_*`.

The calorimeter cell branch `cell_conv_el_idx` points into this table. Energy
deposited by a particle is attributed to the *outermost* recorded conversion
leg among its ancestors (the leg itself included), and `cell_conv_el_idx` is
the leg with the largest attributed energy in the cell, or −1 if no energy in
the cell comes from a recorded leg.

- The shower of a primary-photon leg is attributed to that leg as a whole,
  including the showers of recorded secondary conversions inside it. The energy
  attributed to it is the same as in schema 1 and does not depend on
  `Conversion_secondary_photon_e_min_mev`.
- Legs of a conversion inside the shower of another recorded leg receive no
  energy and are never the target of a cell.
- A secondary conversion outside the shower of any other recorded leg (e.g. a
  photon from a π⁰ produced in a hadronic interaction, or a bremsstrahlung
  photon of a primary electron) is attributed to its own legs. Where such a
  shower deposits more energy in a cell than a primary-photon leg, the cell
  points to the secondary leg.

To attribute a cell to the primary particle, use
`conv_vtx_in_ancestor_idx[conv_el_vtx_idx[cell_conv_el_idx]]`.

Columns: see [Outgoing-particle columns](#outgoing-particle-columns), with prefix `conv_el_`.

## `nucl_vtx_*`: nuclear interactions

One row per hadronic inelastic interaction, capture or fission (hadronic
subtypes 121, 131, 141) of a hadron (Geant4 particle type baryon, meson or
nucleus) inside the tracker region, where the hadron is either

- a Geant4 primary (always recorded), or
- a secondary hadron with kinetic energy ≥ `Truth_record.Nuclear_secondary_hadron_ekin_min_mev`
  (recorded only if `Truth_record.Record_secondary_nuclear_interactions` is true).

| Branch | Type | Description |
|---|---|---|
| `nucl_vtx_x`, `_y`, `_z` | float | Interaction point. |
| `nucl_vtx_in_*` | | Incoming hadron, see [Incoming-particle columns](#incoming-particle-columns). |
| `nucl_vtx_process_subtype` | int | Hadronic process subtype (121 inelastic, 131 capture, 141 fission). |
| `nucl_vtx_material_idx`, `nucl_vtx_volume_idx` | int | Material and logical volume at the interaction point. |
| `nucl_vtx_n_dau_total` | int | Number of secondaries produced by the interaction, before any storage threshold. |
| `nucl_vtx_n_dau_charged` | int | Number of charged secondaries among them. |

A hadron can survive an inelastic interaction in Geant4 and continue as the same
track. Every interaction inside the tracker region is a separate row, whereas
`particle_end_*` describes only where the track finally ended, which may be
outside the tracker region.

## `nucl_dau_*`: nuclear-interaction daughters

One row per secondary of a recorded nuclear interaction with kinetic energy ≥
`Truth_record.Nuclear_daughter_ekin_min_mev` (0 by default, i.e. all). Nuclear
fragments use the PDG nuclear code 100ZZZAAAI.

Columns: see [Outgoing-particle columns](#outgoing-particle-columns), with prefix `nucl_dau_`.

## Shared column layouts

### Incoming-particle columns

Used by `conv_vtx_*` and `nucl_vtx_*` with the prefix `<vtx>_in_`. They describe
the incoming particle immediately before the interaction.

| Suffix | Type | Description |
|---|---|---|
| `in_pdgid` | int | PDG ID (22 for conversions). |
| `in_e` | float | Total energy. |
| `in_px`, `in_py`, `in_pz` | float | Momentum. |
| `in_particle_idx` | int | Row in `particle_*` if the incoming particle is a primary, else -1. |
| `in_ancestor_idx` | int | Row in `particle_*` of the primary it descends from (itself if primary). |
| `in_creator_type`, `in_creator_subtype` | int | Process that created the incoming particle; -1 for primaries. Distinguishes e.g. bremsstrahlung, decay and hadronic photons. |
| `in_prod_x`, `in_prod_y`, `in_prod_z` | float | Production vertex of the incoming particle. |
| `in_x0`, `in_lambda` | float | Material traversed by the incoming particle inside the tracker region from its production up to the interaction, in X0 and λ_I. |
| `in_g4_track_id` | int | Geant4 track ID (information only). |

### Outgoing-particle columns

Used by `conv_el_*` and `nucl_dau_*`.

| Suffix | Type | Description |
|---|---|---|
| `vtx_idx` | int | Row in the parent vertex table (`conv_vtx_*` or `nucl_vtx_*`). |
| `pdgid` | int | PDG ID. |
| `q` | float | Charge in units of e. |
| `e` | float | Total energy at production. |
| `px`, `py`, `pz` | float | Momentum at production. |
| `g4_track_id` | int | Geant4 track ID (information only). |
| `calo_reached` | int | 1 if the particle left the tracker region, i.e. reached the calorimeter front face. |
| `calo_e` | float | Total energy when leaving the tracker region; NaN if `calo_reached` is 0. |
| `calo_x`, `calo_y`, `calo_z` | float | Position when leaving the tracker region; NaN if `calo_reached` is 0. |

## `node_*`: generator truth graph

Written when `Save_truth_particle_graph` is true. One row per generator particle
of interest (final-state particles and their relevant ancestors).

| Branch | Type | Description |
|---|---|---|
| `node_idx` | int | Position of the node in the internal list of interesting particles (not a link). |
| `node_pdg_id` | int | PDG ID. |
| `node_pt`, `node_eta`, `node_phi`, `node_m` | float | Generator kinematics (MeV). η is NaN for particles along the beam axis. |
| `node_isfinal` | int | Generator status code. |
| `node_final_state_idx` | int | Row in `particle_*` if the node was passed to Geant4, else -1. |
| `node_prodx`, `_prody`, `_prodz` | float | Production vertex; NaN if none (incoming beams). |
| `node_decx`, `_decy`, `_decz` | float | Decay vertex; NaN if none (final-state particles). |
| `node_edge_start`, `node_edge_end` | int | Parent → child edges between `node_*` rows. Each edge appears once. |

## Reconstruction-level tables (unchanged)

`cell_*` (calorimeter cells; `cell_parent_idx` → `particle_*`, `cell_conv_el_idx` →
`conv_el_*`), `track_*` and `track_{x,y,z}_layer_<i>` (emulated tracks;
`track_parent_idx` → `particle_*`), `topo_*`, `true_jet_*`, `topo_jet_*`,
graph-edge branches, and, if enabled, `pflow_*`, `pflow_jet_*` and `supercluster_*`.

With `Do_superclustering` enabled, conversion legs are also passed to the track
emulation. Their tracks have `track_parent_idx ≥ len(particle_*)`; the leg's row
in `conv_el_*` is `track_parent_idx - len(particle_*)`.

## Run-level trees

### `RunInfo` (one entry)

| Branch | Type | Description |
|---|---|---|
| `schema_version` | int | 2 for this document. |
| `cocoa_git_commit`, `cocoa_git_describe` | string | Commit of the COCOA build. |
| `cocoa_git_dirty` | int | 1 if tracked files had uncommitted changes at build time. |
| `geant4_version`, `root_version`, `pythia8_version` | string | Library versions. |
| `physics_list` | string | Geant4 physics list. |
| `command_line` | string | Full command line. |
| `config_path`, `config_json` | string | Path and verbatim content of the JSON config. |
| `macro_path`, `macro_text` | string | Path and verbatim content of the macro. |
| `output_path` | string | Output ROOT file as given. |
| `run_seed` | int64 | Random seed. |
| `n_events_requested` | int | `--nevents`; -1 if taken from the macro. |
| `n_events_processed` | int | Number of events in this file. |
| `type_of_running` | string | `Type_of_running` from the config. |
| `geometry_fingerprint` | string | MD5 over all placed volumes (path, placement, solid parameters, material). Equal fingerprints mean identical geometries. |
| `tracker_region_r_max_mm`, `tracker_region_abs_z_max_mm` | float | Extent of the tracker region. |
| `truth_record_secondary_conversions`, `truth_record_conversion_secondary_photon_e_min_mev` | int, float | Effective conversion recording settings. |
| `truth_record_secondary_nuclear`, `truth_record_nuclear_secondary_hadron_ekin_min_mev`, `truth_record_nuclear_daughter_ekin_min_mev` | int, float, float | Effective nuclear-interaction recording settings. |
| `fiducial_pt_min_gev`, `fiducial_e_min_gev`, `fiducial_eta_max`, `fiducial_dr_cut` | float | Effective generator-level fiducial cuts. |
| `magnetic_field_tesla` | float | Solenoid field. |
| `use_inner_detector`, `use_detailed_tracker_material` | int | Geometry switches. |
| `start_time_utc`, `end_time_utc`, `host` | string | Run timestamps (ISO 8601) and host name. |

### `Materials` (one entry per Geant4 material)

| Branch | Type | Description |
|---|---|---|
| `material_idx` | int | Index used by `*_material_idx` columns. |
| `name` | string | Material name. |
| `density_g_cm3` | float | Density. |
| `x0_mm`, `lambda_mm` | float | Radiation length and nuclear interaction length. |
| `z_eff`, `a_eff_g_mol` | float | Mass-fraction-weighted mean Z and A (the same definition as ROOT's `TGeoMixture`). |

### `Volumes` (one entry per logical volume)

| Branch | Type | Description |
|---|---|---|
| `volume_idx` | int | Index used by `*_volume_idx` columns. |
| `name` | string | Logical volume name. |
| `material_idx` | int | Row in `Materials`. |
| `in_tracker` | int | 1 if the volume belongs to the tracker region. |

## Configuration keys

```jsonc
"Fiducial_cuts": {
    "pt_min_gev": 0.5,   // generator particles below this pT are not simulated
    "e_min_gev": 0.0,    // generator particles below this energy are not simulated
    "eta_max": 3.0,
    "dR_cut": -1
},
"Truth_record": {
    "Record_secondary_conversions": true,
    "Conversion_secondary_photon_e_min_mev": 100,
    "Record_secondary_nuclear_interactions": true,
    "Nuclear_secondary_hadron_ekin_min_mev": 300,
    "Nuclear_daughter_ekin_min_mev": 0
}
```

The values shown for `Truth_record` are the defaults used when the block is absent.

## Geometry export

With `Save_geometry_json` / `Save_geometry_gdml`, the geometry is written to
`Geometry_json_path` / `Geometry_gdml_path` if set, otherwise next to the output
ROOT file as `<output without .root>.geometry.json` / `.geometry.gdml`. The
output path includes a `--output` override. If neither path is available COCOA
stops with an error. An existing GDML file at the target path is replaced.

The geometry JSON (`format: cocoa_geometry_material_dump_v3`) lists every placed
volume with its global placement (object rotation, i.e. the rotation applied to
the volume) and its material, including `z_eff`, `a_eff_g_mol` and the element
composition.

## Process codes

| Type | Subtype | Process |
|---|---|---|
| 1 Transportation | 91 | left the world volume |
| 2 Electromagnetic | 11 | Rayleigh scattering |
| 2 Electromagnetic | 12 | photoelectric effect |
| 2 Electromagnetic | 13 | Compton scattering |
| 2 Electromagnetic | 14 | γ → e⁺e⁻ conversion |
| 2 Electromagnetic | 2, 3 | ionisation, bremsstrahlung |
| 4 Hadronic | 111 | hadronic elastic |
| 4 Hadronic | 121 | hadronic inelastic |
| 4 Hadronic | 131 | capture |
| 4 Hadronic | 141 | fission |
| 4 Hadronic | 151 | capture at rest |
| 6 Decay | 201 | decay |

## Common selections

- **Primary photons:** `particle_pdgid == 22`.
- **Primary photons that converted (anywhere):** additionally
  `particle_end_process_type == 2 && particle_end_process_subtype == 14`; the
  conversion point is `particle_end_{x,y,z}`, the photon energy there `particle_end_e`.
- **Conversions inside the tracker with their legs:** rows of `conv_vtx_*`; the
  photon is a primary if `conv_vtx_in_particle_idx >= 0`; the legs are
  `conv_el_*[conv_vtx_ep_idx]` and `conv_el_*[conv_vtx_em_idx]`.
- **Nuclear interactions of primary hadrons:** rows of `nucl_vtx_*` with
  `nucl_vtx_in_particle_idx >= 0`; their daughters are the `nucl_dau_*` rows with
  `nucl_dau_vtx_idx` equal to the vertex row.

## Changes from schema version 1

| Version 1 | Version 2 |
|---|---|
| `primary_photon_*` (16 branches) | `particle_*` rows with `particle_pdgid == 22`; see the mapping below. |
| `nucl_int_*` (one row per daughter of a primary hadron) | `nucl_vtx_*` (one row per interaction) and `nucl_dau_*` (daughters). |
| `conv_el_primary_photon_idx` | `conv_vtx_in_particle_idx[conv_el_vtx_idx]` |
| `conv_el_prod_{x,y,z}` | `conv_vtx_{x,y,z}[conv_el_vtx_idx]` |
| `conv_el_*` only for primary photons | also secondary photons above threshold |
| `cell_conv_el_idx` only primary-photon legs | outermost recorded leg, which can also be a secondary-photon leg (see `conv_el_*`) |
| truth-graph edges with duplicates | each edge once |
| no run metadata | `RunInfo`, `Materials`, `Volumes`, `event_id`, `run_seed` |

Mapping of the former photon branches:

| `primary_photon_*` | Version 2 |
|---|---|
| `particle_idx` | the `particle_*` row itself |
| `track_id` | `particle_g4_track_id` |
| `energy`, `eta`, `phi`, `prod_{x,y,z}` | `particle_e`, `particle_eta`, `particle_phi`, `particle_prod_{x,y,z}` |
| `end_{x,y,z}`, `path_length` | `particle_end_{x,y,z}`, `particle_track_length` |
| `converted` | `particle_end_process_type == 2 && particle_end_process_subtype == 14` |
| `termination_process_type`, `_subtype` | `particle_end_process_type`, `_subtype` |
| `track_status` | dropped (always 2, "stop and kill", at the end of a track) |

## Validation tools

- `util/validate_output.py FILE.root` checks a file against this document:
  table lengths, every link, the run-level trees, and physics consistency of
  the truth record (e.g. e⁺ + e⁻ energy equals the photon energy at the
  conversion, the conversion point of a primary photon equals its
  `particle_end_*`, calorimeter-face points lie on the tracker-region boundary).
- `util/compare_outputs.py REFERENCE.root CANDIDATE.root` compares two files made
  with the same configuration and seed, branch by branch. Against a schema-1
  reference it also checks the former `primary_photon_*` and `conv_el_*` branches
  against their schema-2 replacements, and that a cell with a schema-1
  `cell_conv_el_idx` keeps it or now points to a secondary-photon leg. Expected
  differences between schema 1 and 2 for the same seed: `conv_el_*`
  (additional legs of secondary photons) and `cell_conv_el_idx` (cells that
  were −1, or where a secondary-photon shower now outweighs the primary-photon
  leg); pass them with `--ignore`.
