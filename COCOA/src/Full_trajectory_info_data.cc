#include "Full_trajectory_info_data.hh"
#include "DetectorConstruction.hh"
#include "TruthRecordGraph.hh"
#include <algorithm>
#include <cassert>

Full_trajectory_info_data::Full_trajectory_info_data(/* args */)
{;}


void Full_trajectory_info_data::make_pseudo_jet_particles()
{
    int size_trajectories_list = fAllTrajectoryInfo.size();
    for (int iparticle = 0; iparticle < size_trajectories_list; iparticle++)
    {
        fastjet::PseudoJet particle(fAllTrajectoryInfo.at(iparticle).fMomentum.x(),
                                                  fAllTrajectoryInfo.at(iparticle).fMomentum.y(),
                                                  fAllTrajectoryInfo.at(iparticle).fMomentum.z(),
                                                  fAllTrajectoryInfo.at(iparticle).fEnergy);
        particle.set_user_index(iparticle);
        jets_objects.push_back(particle);
    }
}

int Full_trajectory_info_data::DeltaR_iso(float px, float py, float pz ,size_t idx_m,int pdgid) {

    size_t n_particles = fAllTrajectoryInfo.size();
    
    float ptot_origin = sqrt(px*px+py*py+pz*pz);
    float theta_origin = acos(pz/ptot_origin);
    float eta_origin = -1*log(tan(theta_origin/2.));
    float phi_origin = atan(py/px);

    float pT =  sqrt(px*px+py*py) / 1e3;
    float iso_radius = std::min(10./pT, 0.2);

    int isIso = 0;
    //consider all neutral to be non-isolated
    if(fabs(pdgid) == 22 || fabs(pdgid) == 310 || fabs(pdgid) == 3122 || fabs(pdgid) == 3322 || fabs(pdgid) == 2112) return 0;

    float sum_pt = 0;
    for(size_t part_i=0; part_i < n_particles; part_i++){
	if (idx_m == part_i) continue;
	int pdg_id_loop = fAllTrajectoryInfo.at(part_i).fPDGCode;
         //exclude from the computation neutral particles
         if(fabs(pdg_id_loop) == 22 || fabs(pdg_id_loop) == 310 || fabs(pdg_id_loop) == 3122 || fabs(pdg_id_loop) == 3322 || fabs(pdg_id_loop) == 2112) continue;
	float px_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.x();
	float py_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.y();
	float pz_loop   = fAllTrajectoryInfo.at(part_i).fMomentum.z();
	float phi       = atan(py_loop/px_loop);
	float ptot_loop = sqrt(px_loop*px_loop+py_loop*py_loop+pz_loop*pz_loop);
	float theta     = acos(pz_loop/ptot_loop);
	float eta       = -1*log(tan(theta/2.));
	float dphi      = acos(cos(phi_origin-phi));
	float deta      = eta_origin-eta;
	float dr        = sqrt(dphi*dphi+deta*deta);
	if(dr<iso_radius){
	    sum_pt = sqrt(px_loop*px_loop+py_loop*py_loop) /1e3 + sum_pt;
	}
    }
    if(sum_pt/pT < 0.06) isIso = 1;
    return isIso;
    
}

void Full_trajectory_info_data::fill_var(){
    // Link particle_* rows and the truth graph; this also turns
    // node_final_state_idx into a particle_* row index.
    std::vector<int> generator_idx;
    generator_idx.reserve(fAllTrajectoryInfo.size());
    for (const FullTrajectoryInfo &trajectory : fAllTrajectoryInfo)
        generator_idx.push_back(trajectory.fGeneratorIdx);
    TruthRecordGraph::GetInstance().link_to_particles(generator_idx, particle_truth_node_idx);

    int size_trajectories_list = fAllTrajectoryInfo.size();
    for (int iparticle = 0; iparticle < size_trajectories_list; iparticle++)
    {
        const FullTrajectoryInfo &trajectory = fAllTrajectoryInfo.at(iparticle);
        float px    = trajectory.fMomentum.x();
        float py    = trajectory.fMomentum.y();
        float pz    = trajectory.fMomentum.z();
        float p     = sqrt( px * px + py * py + pz * pz ); // Have some 4-momentum object at some point (but would like to circumvent ROOT). PR.
        float theta = acos( pz / p );
        float eta   = -1.0 * log( tan( 0.5 * theta ) );
        float phi   = GetPhi( px, py );

        particle_pdgid.push_back(trajectory.fPDGCode);
        particleisIso.push_back(DeltaR_iso(px, py, pz,iparticle,int(trajectory.fPDGCode)));
        particle_pt.push_back(sqrtf(px*px + py*py));
        particle_eta.push_back(eta);
        particle_phi.push_back(phi);
        particle_e.push_back(trajectory.fEnergy);
        particle_prod_x.push_back(trajectory.fVertexPosition.x());
        particle_prod_y.push_back(trajectory.fVertexPosition.y());
        particle_prod_z.push_back(trajectory.fVertexPosition.z());

        particle_q.push_back(trajectory.fPDGCharge);
        particle_g4_track_id.push_back(trajectory.fTrackID);
        particle_end_x.push_back(trajectory.fEndPosition.x());
        particle_end_y.push_back(trajectory.fEndPosition.y());
        particle_end_z.push_back(trajectory.fEndPosition.z());
        particle_end_e.push_back(trajectory.fEndEnergy);
        particle_end_px.push_back(trajectory.fEndMomentum.x());
        particle_end_py.push_back(trajectory.fEndMomentum.y());
        particle_end_pz.push_back(trajectory.fEndMomentum.z());
        particle_end_process_type.push_back(trajectory.fEndProcessType);
        particle_end_process_subtype.push_back(trajectory.fEndProcessSubType);
        particle_end_in_tracker.push_back(trajectory.fEndInTracker);
        particle_end_material_idx.push_back(trajectory.fEndMaterialIdx);
        particle_end_volume_idx.push_back(trajectory.fEndVolumeIdx);
        particle_track_length.push_back(trajectory.fTrackLength);
        particle_x0_tracker.push_back(trajectory.fX0Tracker);
        particle_lambda_tracker.push_back(trajectory.fLambdaTracker);
        particle_n_scatters_tracker.push_back(trajectory.fNScattersTracker);

        caloExtrapolEta.push_back(trajectory.caloExtrapolEta);
        caloExtrapolPhi.push_back(trajectory.caloExtrapolPhi);

        idExtrapolEta.push_back(trajectory.idExtrapolEta);
        idExtrapolPhi.push_back(trajectory.idExtrapolPhi);

        particle_to_track.push_back(-1);
    }

    Tracks_data &tracks_obj = Tracks_data::GetLow();
    int size_tracks = tracks_obj.Tracks_list.size();
    for (int itrack = 0; itrack < size_tracks; itrack++)
    {
	    int trackParticleId = tracks_obj.Tracks_list[itrack].nFinal_State_Particles;
	    // Tracks of conversion legs (superclustering) do not belong to a primary.
	    if (trackParticleId < (int)particle_to_track.size())
	        particle_to_track.at(trackParticleId) = itrack;
    }
}



void Full_trajectory_info_data::clear()
{
    fAllTrajectoryInfo.clear();
    fAllConvElectrons.clear();
    for (auto *v : {&particle_pdgid, &particleisIso, &particle_g4_track_id, &particle_truth_node_idx,
                    &particle_end_process_type, &particle_end_process_subtype, &particle_end_in_tracker,
                    &particle_end_material_idx, &particle_end_volume_idx, &particle_n_scatters_tracker,
                    &particle_to_track})
        v->clear();
    for (auto *v : {&particle_pt, &particle_eta, &particle_phi, &particle_e,
                    &particle_prod_x, &particle_prod_y, &particle_prod_z, &particle_q,
                    &particle_end_x, &particle_end_y, &particle_end_z, &particle_end_e,
                    &particle_end_px, &particle_end_py, &particle_end_pz,
                    &particle_track_length, &particle_x0_tracker, &particle_lambda_tracker,
                    &caloExtrapolEta, &caloExtrapolPhi, &idExtrapolEta, &idExtrapolPhi,
                    &particle_dep_energy})
        v->clear();
    jets_objects.clear();
}

void Full_trajectory_info_data::set_tree_branches(TTree *outTree)
{
    outTree->Branch("particle_pdgid"   ,           "vector<int>", &particle_pdgid);
    outTree->Branch("particle_isIso"   ,           "vector<int>", &particleisIso);
    outTree->Branch("particle_pt"  ,               "vector<float>", &particle_pt);
    outTree->Branch("particle_eta"  ,              "vector<float>", &particle_eta);
    outTree->Branch("particle_phi"  ,              "vector<float>", &particle_phi);
    outTree->Branch("particle_e"   ,               "vector<float>", &particle_e);
    outTree->Branch("particle_prod_x"  ,           "vector<float>", &particle_prod_x);
    outTree->Branch("particle_prod_y"  ,           "vector<float>", &particle_prod_y);
    outTree->Branch("particle_prod_z"  ,           "vector<float>", &particle_prod_z);
    outTree->Branch("particle_track_idx",          "vector<int>", &particle_to_track);
    outTree->Branch("particle_eta_extrap_calo",    "vector<float>", &caloExtrapolEta);
    outTree->Branch("particle_phi_extrap_calo",    "vector<float>", &caloExtrapolPhi);
    outTree->Branch("particle_eta_extrap_its",     "vector<float>", &idExtrapolEta);
    outTree->Branch("particle_phi_extrap_its",     "vector<float>", &idExtrapolPhi);
    outTree->Branch("particle_dep_energy"  ,       "vector<float>", &particle_dep_energy);

    outTree->Branch("particle_q",                   "vector<float>", &particle_q);
    outTree->Branch("particle_g4_track_id",         "vector<int>",   &particle_g4_track_id);
    outTree->Branch("particle_truth_node_idx",      "vector<int>",   &particle_truth_node_idx);
    outTree->Branch("particle_end_x",               "vector<float>", &particle_end_x);
    outTree->Branch("particle_end_y",               "vector<float>", &particle_end_y);
    outTree->Branch("particle_end_z",               "vector<float>", &particle_end_z);
    outTree->Branch("particle_end_e",               "vector<float>", &particle_end_e);
    outTree->Branch("particle_end_px",              "vector<float>", &particle_end_px);
    outTree->Branch("particle_end_py",              "vector<float>", &particle_end_py);
    outTree->Branch("particle_end_pz",              "vector<float>", &particle_end_pz);
    outTree->Branch("particle_end_process_type",    "vector<int>",   &particle_end_process_type);
    outTree->Branch("particle_end_process_subtype", "vector<int>",   &particle_end_process_subtype);
    outTree->Branch("particle_end_in_tracker",      "vector<int>",   &particle_end_in_tracker);
    outTree->Branch("particle_end_material_idx",    "vector<int>",   &particle_end_material_idx);
    outTree->Branch("particle_end_volume_idx",      "vector<int>",   &particle_end_volume_idx);
    outTree->Branch("particle_track_length",        "vector<float>", &particle_track_length);
    outTree->Branch("particle_x0_tracker",          "vector<float>", &particle_x0_tracker);
    outTree->Branch("particle_lambda_tracker",      "vector<float>", &particle_lambda_tracker);
    outTree->Branch("particle_n_scatters_tracker",  "vector<int>",   &particle_n_scatters_tracker);
}

void Full_trajectory_info_data::SetParticleDepEnergy( const std::vector<float> &_particle_dep_energies) {

    particle_dep_energy = _particle_dep_energies;
    
}
