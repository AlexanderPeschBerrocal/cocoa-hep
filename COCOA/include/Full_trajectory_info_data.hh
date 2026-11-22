#ifndef __FULL_TRAJECTORY_INFO_VAR_H__
#define __FULL_TRAJECTORY_INFO_VAR_H__

#include "vector"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include "TTree.h"
#include "Tracks_data.hh"
#include "Track_var.hh"
#include "fastjet/ClusterSequence.hh"


struct FullTrajectoryInfo {
	G4int fTrackID = -1;
	G4int fParentID = -1;
	G4int fPDGCode = 0;
	G4int fGeneratorIdx = -1; // index in TruthRecordGraph::m_final_state_particles, -1 if none
	G4double fPDGCharge = 0.0;
	G4ThreeVector fMomentum;
	G4ThreeVector fMomentumDir;
	G4double    fEnergy = 0.0;
	G4double    fMass = 0.0;
	G4ThreeVector fVertexPosition;
	G4double fGlobalTime = 0.0;

	// State at the end of the track (filled for Geant4 primaries).
	G4ThreeVector fEndPosition;
	G4double fEndEnergy = 0.0;
	G4ThreeVector fEndMomentum;
	G4int fEndProcessType = -1;
	G4int fEndProcessSubType = -1;
	G4int fEndInTracker = 0;
	G4int fEndMaterialIdx = -1;
	G4int fEndVolumeIdx = -1;
	G4double fTrackLength = 0.0;
	G4double fX0Tracker = 0.0;
	G4double fLambdaTracker = 0.0;
	G4int fNScattersTracker = 0;

	float caloExtrapolMaxEkin = 0.f;
	float caloExtrapolEta = 0.f;
	float caloExtrapolPhi = 0.f;

	float idExtrapolMaxEkin = 0.f;
	float idExtrapolEta = 0.f;
	float idExtrapolPhi = 0.f;

	std::vector<G4int> vTrackID;
	std::vector<G4int> vParentID;
	std::vector<G4ThreeVector> vTrackMomentumDir;
	std::vector<G4ThreeVector> vTrackPos;
	std::vector<G4double> vTrackTime;
	std::vector<G4double> vTrackPDGID;

	bool is_conversion_track = false;

	// fastjet::PseudoJet particle_object;
};

class Full_trajectory_info_data
{
public:
    Full_trajectory_info_data();
	void make_pseudo_jet_particles();
    void clear();
void fill_var();
	void set_tree_branches(TTree *outTree);
    static Full_trajectory_info_data &GetInstance()
    {
        static Full_trajectory_info_data trajectories;
        return trajectories;
    };
    // Geant4 primaries, one entry per particle_* row.
    std::vector < FullTrajectoryInfo> fAllTrajectoryInfo;
    // Conversion legs, one entry per conv_el_* row (see Truth_record_data).
    // Used for the cell attribution (cell_conv_el_idx) and superclustering.
    std::vector < FullTrajectoryInfo> fAllConvElectrons;
	std::vector <int>	particle_to_track;
	std::vector <fastjet::PseudoJet> jets_objects;
    void SetParticleDepEnergy( const std::vector<float> &_particle_dep_energies);
    int  DeltaR_iso(float px, float py, float pz ,size_t idx_m,int particle_loop_pdgid);
private:
	std::vector <int>	particle_pdgid;
    std::vector <int>   particleisIso;
	std::vector <float> particle_pt;
	std::vector <float> particle_eta;
	std::vector <float> particle_phi;
	std::vector <float> particle_e;
	std::vector <float> particle_prod_x;
	std::vector <float> particle_prod_y;
	std::vector <float> particle_prod_z;
    std::vector<float>  particle_dep_energy;

    // Schema version 2 additions, see docs/output_schema.md.
    std::vector<float> particle_q;
    std::vector<int>   particle_g4_track_id;
    std::vector<int>   particle_truth_node_idx;
    std::vector<float> particle_end_x;
    std::vector<float> particle_end_y;
    std::vector<float> particle_end_z;
    std::vector<float> particle_end_e;
    std::vector<float> particle_end_px;
    std::vector<float> particle_end_py;
    std::vector<float> particle_end_pz;
    std::vector<int>   particle_end_process_type;
    std::vector<int>   particle_end_process_subtype;
    std::vector<int>   particle_end_in_tracker;
    std::vector<int>   particle_end_material_idx;
    std::vector<int>   particle_end_volume_idx;
    std::vector<float> particle_track_length;
    std::vector<float> particle_x0_tracker;
    std::vector<float> particle_lambda_tracker;
    std::vector<int>   particle_n_scatters_tracker;

    std::vector<float> caloExtrapolEta;     // pseudorapidity of the same-pdgid daughter track of maximum kinetic energy in the calorimeter
    std::vector<float> caloExtrapolPhi;     // phi of the same-pdgid daughter track of maximum kinetic energy in the calorimeter
    std::vector<float> idExtrapolEta;     // pseudorapidity of the same-pdgid daughter track of maximum kinetic energy in the outermost id layer
    std::vector<float> idExtrapolPhi;     // phi of the same-pdgid daughter track of maximum kinetic energy in the outermost id layer

};

#endif // __FULL_TRAJECTORY_INFO_VAR_H__
