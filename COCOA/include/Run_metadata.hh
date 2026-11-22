#ifndef __RUN_METADATA_H__
#define __RUN_METADATA_H__

//
// Run-level trees written next to the event tree (docs/output_schema.md):
//   RunInfo    one entry: provenance, versions, settings, geometry fingerprint
//   Materials  one entry per Geant4 material  (target of *_material_idx)
//   Volumes    one entry per logical volume   (target of *_volume_idx)
//

class Run_metadata
{
public:
    // Records the start time of the run.
    static void BeginRun();
    // Writes the three trees into the current ROOT directory.
    static void Write(int n_events_processed);

    // Schema version of the output, see docs/output_schema.md.
    static constexpr int kSchemaVersion = 2;
};

#endif // __RUN_METADATA_H__
