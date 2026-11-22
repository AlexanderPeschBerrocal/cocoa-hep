#ifndef __TRUTH_RECORD_DATA_H__
#define __TRUTH_RECORD_DATA_H__

//
// Truth record of interactions inside the tracker region (schema version 2,
// see docs/output_schema.md):
//
//   conv_vtx_*  photon conversions           conv_el_*   their e+/e- legs
//   nucl_vtx_*  hadronic inelastic/capture/  nucl_dau_*  their secondaries
//               fission interactions
//
// plus the per-track bookkeeping needed to fill these tables and the end-state
// columns of particle_*: material traversed inside the tracker, number of
// scatters, primary ancestor of every track, and the state of selected
// secondaries when they leave the tracker region (calorimeter front face).
//
// Geant4 tracks one track at a time and stacks the secondaries of a track only
// after that track has ended, so
//   - per-track accumulators can be kept for the current track only, and
//   - a vertex created while stepping its incoming particle always exists
//     before its outgoing particles start to be tracked.
//

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "G4ThreeVector.hh"
#include "TTree.h"

class G4LogicalVolume;
class G4Region;
class G4Step;
class G4StepPoint;
class G4Track;
class G4VProcess;
struct FullTrajectoryInfo;

// Incoming particle immediately before an interaction.
struct Truth_incoming
{
    int   pdgid           = 0;
    float e               = 0.f;
    float px              = 0.f;
    float py              = 0.f;
    float pz              = 0.f;
    int   particle_idx    = -1;
    int   ancestor_idx    = -1;
    int   creator_type    = -1;
    int   creator_subtype = -1;
    float prod_x          = 0.f;
    float prod_y          = 0.f;
    float prod_z          = 0.f;
    float x0              = 0.f;
    float lambda          = 0.f;
    int   g4_track_id     = -1;
};

struct Truth_vertex
{
    float          x = 0.f;
    float          y = 0.f;
    float          z = 0.f;
    Truth_incoming in;
    int            material_idx = -1;
    int            volume_idx   = -1;
    // Conversions only.
    int ep_idx = -1;
    int em_idx = -1;
    // Nuclear interactions only.
    int process_subtype = -1;
    int n_dau_total     = 0;
    int n_dau_charged   = 0;
};

// Particle created at a recorded vertex.
struct Truth_outgoing
{
    int   vtx_idx      = -1;
    int   pdgid        = 0;
    float q            = 0.f;
    float e            = 0.f;
    float px           = 0.f;
    float py           = 0.f;
    float pz           = 0.f;
    int   g4_track_id  = -1;
    int   calo_reached = 0;
    float calo_e;
    float calo_x;
    float calo_y;
    float calo_z;
    Truth_outgoing();
};

class Truth_record_data
{
public:
    static Truth_record_data &GetInstance()
    {
        static Truth_record_data instance;
        return instance;
    }

    // Per event.
    void clear();
    void set_tree_branches(TTree *outTree);
    void set_event_tree_branches(TTree *outTree);
    void fill_var();
    void SetEventInfo(int event_id, long long run_seed);

    // Called from the user actions.
    void RegisterPrimary(int g4_track_id, int particle_row);
    // Returns the conv_el_* row created for this track, or -1.
    int  BeginTrack(const G4Track *track);
    void Step(const G4Step *step);
    void EndPrimaryTrack(const G4Track *track, FullTrajectoryInfo &primary);

    // Geometry helpers, also used for the run-level metadata.
    const G4Region *TrackerRegion();
    bool InTracker(const G4LogicalVolume *lv);
    int  VolumeIndex(const G4LogicalVolume *lv);

    std::vector<Truth_vertex>   conv_vtx;
    std::vector<Truth_outgoing> conv_el;
    std::vector<Truth_vertex>   nucl_vtx;
    std::vector<Truth_outgoing> nucl_dau;

private:
    Truth_record_data() = default;

    Truth_incoming MakeIncoming(const G4Track *track, const G4StepPoint *pre) const;
    void MakeConversionVertex(const G4Step *step);
    void MakeNuclearVertex(const G4Step *step, const G4VProcess *process);
    static void FillOutgoing(const G4Track *track, int vtx_idx, Truth_outgoing &out);

    enum class Watched { kConvEl, kNuclDau };

    struct Pending_nuclear
    {
        int               row;
        const G4VProcess *process;
        G4ThreeVector     position;
    };

    // Current track.
    int    m_current_track_id = -1;
    double m_x0               = 0.0;
    double m_lambda           = 0.0;
    int    m_n_scatters       = 0;

    // Per event, keyed by Geant4 track ID.
    std::unordered_map<int, int> m_primary_row;   // primary -> particle_* row
    std::unordered_map<int, int> m_ancestor_row;  // any track -> particle_* row of its primary
    std::unordered_map<int, int> m_pending_conv;  // converted photon -> conv_vtx_* row
    std::unordered_map<int, std::vector<Pending_nuclear>> m_pending_nucl; // hadron -> its recorded nuclear vertices
    std::unordered_map<int, std::pair<Watched, int>> m_exit_watch; // outgoing particle -> table row

    // Per run.
    const G4Region *m_tracker_region = nullptr;
    std::unordered_map<const G4LogicalVolume *, int> m_volume_index;

    // Event-level branches.
    int       m_event_id = -1;
    long long m_run_seed = 0;

    // Branch buffers.
    struct Vertex_branches
    {
        std::vector<float> x, y, z;
        std::vector<int>   in_pdgid;
        std::vector<float> in_e, in_px, in_py, in_pz;
        std::vector<int>   in_particle_idx, in_ancestor_idx, in_creator_type, in_creator_subtype;
        std::vector<float> in_prod_x, in_prod_y, in_prod_z, in_x0, in_lambda;
        std::vector<int>   in_g4_track_id;
        std::vector<int>   material_idx, volume_idx;
        void clear();
        void fill(const std::vector<Truth_vertex> &rows);
        void set_branches(TTree *tree, const std::string &prefix);
    };
    struct Outgoing_branches
    {
        std::vector<int>   vtx_idx, pdgid;
        std::vector<float> q, e, px, py, pz;
        std::vector<int>   g4_track_id, calo_reached;
        std::vector<float> calo_e, calo_x, calo_y, calo_z;
        void clear();
        void fill(const std::vector<Truth_outgoing> &rows);
        void set_branches(TTree *tree, const std::string &prefix);
    };

    Vertex_branches   m_conv_vtx_branches;
    std::vector<int>  m_conv_vtx_ep_idx, m_conv_vtx_em_idx;
    Outgoing_branches m_conv_el_branches;
    Vertex_branches   m_nucl_vtx_branches;
    std::vector<int>  m_nucl_vtx_process_subtype, m_nucl_vtx_n_dau_total, m_nucl_vtx_n_dau_charged;
    Outgoing_branches m_nucl_dau_branches;
};

#endif // __TRUTH_RECORD_DATA_H__
