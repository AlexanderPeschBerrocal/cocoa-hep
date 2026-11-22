#ifndef __H02TRACKINGACTION_H__
#define __H02TRACKINGACTION_H__

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file runAndEvent/H02/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
#include "OutputRunAction.hh"
// #include "TrackInformation.hh"
#include "TrackingAction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Trajectory.hh"
#include "G4RunManager.hh"
#include "G4String.hh"
// #include "DataStorage.hh"

#include "SteppingAction.hh"
//#include "FullTrajectoryInfo.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
using namespace std;


#include <cassert>

#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "PrimaryParticleInformation.hh"
#include "Truth_record_data.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction()
	: G4UserTrackingAction()
{
	;
}

void TrackingAction::PreUserTrackingAction(const G4Track*aTrack) 
{
	Full_trajectory_info_data &trajectories = Full_trajectory_info_data::GetInstance();
	Truth_record_data &truth = Truth_record_data::GetInstance();
	G4int ParentID = aTrack->GetParentID();
	if (aTrack->GetParentID() == 0)
	{
		const G4PrimaryParticle *primary = aTrack->GetDynamicParticle()->GetPrimaryParticle();
	        FullTrajectoryInfo trjInfo;
		trjInfo.is_conversion_track = false;
		trjInfo.fParentID = primary->GetTrackID();
		trjInfo.fTrackID = aTrack->GetTrackID();
		trjInfo.fPDGCharge = primary->GetCharge();
		trjInfo.fPDGCode = aTrack->GetDefinition()->GetPDGEncoding();
		trjInfo.fMomentum = primary->GetMomentum();
		trjInfo.fMomentumDir = primary->GetMomentumDirection();
		trjInfo.fEnergy = primary->GetTotalEnergy();
		trjInfo.fMass = primary->GetMass();
		if (const auto *info = dynamic_cast<const PrimaryParticleInformation *>(primary->GetUserInformation()))
			trjInfo.fGeneratorIdx = info->GetGeneratorIdx();

		trjInfo.fVertexPosition = aTrack->GetVertexPosition();
		trjInfo.fEndPosition = aTrack->GetPosition();
		trjInfo.fGlobalTime = aTrack->GetGlobalTime();
		
		trjInfo.caloExtrapolMaxEkin = 0.0;
		trjInfo.caloExtrapolEta     = trjInfo.fMomentum.getEta();
		trjInfo.caloExtrapolPhi     = GetPhi( trjInfo.fMomentum.x(),
						      trjInfo.fMomentum.y() );
		
		trjInfo.idExtrapolMaxEkin = trjInfo.caloExtrapolMaxEkin;
		trjInfo.idExtrapolEta     = trjInfo.caloExtrapolEta;
		trjInfo.idExtrapolPhi     = trjInfo.caloExtrapolPhi;
		
		trjInfo.vTrackMomentumDir.push_back(aTrack->GetMomentum());
		trjInfo.vTrackID.push_back(aTrack->GetTrackID());
		trjInfo.vParentID.push_back(aTrack->GetParentID());
		trjInfo.vTrackPos.push_back(aTrack->GetPosition());
		trjInfo.vTrackTime.push_back(aTrack->GetGlobalTime());
		trjInfo.vTrackPDGID.push_back(aTrack->GetDefinition()->GetPDGEncoding());

		truth.RegisterPrimary(aTrack->GetTrackID(), (int)trajectories.fAllTrajectoryInfo.size());
		trajectories.fAllTrajectoryInfo.push_back(trjInfo);
	} // if(trj->GetParentID() == 0)
	else
	{
	    for ( std::vector < FullTrajectoryInfo>* _trajectories : { &trajectories.fAllTrajectoryInfo, &trajectories.fAllConvElectrons } )
		{
		    bool foundTraj(false);
		    int mTraj(-1); //, mParent(-1);
		    for (int iTraj = (int)_trajectories->size() - 1; iTraj >= 0; iTraj--)
			{
			    for (int iParent = (int)_trajectories->at(iTraj).vTrackID.size() - 1; iParent >= 0; iParent--)
				{
				    if (ParentID == _trajectories->at(iTraj).vTrackID.at(iParent))
					{
					    foundTraj = true;
					    mTraj = iTraj;
					    break;
					    
					} // if( _trajectories->at(iTraj).vParentID.at(iParent) == ParentID )
				    
				} // for(int iParent = 0; iParent < (int)_trajectories->at(iTraj).vParentID.size(); iParent++  )
			    
			    if (foundTraj)
				break;
			    
			} // for(int iTraj = 0; iTraj < (int)_trajectories->size(); iTraj++)
		    if (foundTraj)
			{
			    _trajectories->at(mTraj).vTrackMomentumDir.push_back(aTrack->GetMomentum());
			    _trajectories->at(mTraj).vParentID.push_back(aTrack->GetParentID());
			    _trajectories->at(mTraj).vTrackID.push_back(aTrack->GetTrackID());
			    _trajectories->at(mTraj).vTrackPos.push_back(aTrack->GetPosition());
			    _trajectories->at(mTraj).vTrackTime.push_back(aTrack->GetGlobalTime());
			    _trajectories->at(mTraj).vTrackPDGID.push_back(aTrack->GetDefinition()->GetPDGEncoding());
			} //  if(foundTraj)
		}
	}	  //if(aTrack->GetParentID() != 0)

	// Truth record (docs/output_schema.md). Returns the conv_el_* row if this
	// track is a leg of a recorded conversion.
	const int conv_el_row = truth.BeginTrack(aTrack);
	if (conv_el_row >= 0)
	{
		// Keep fAllConvElectrons aligned with conv_el_*: the cell attribution
		// (cell_conv_el_idx) and superclustering index this list.
		assert(conv_el_row == (int)trajectories.fAllConvElectrons.size());

		const Truth_outgoing &leg = truth.conv_el[conv_el_row];
		const Truth_vertex &vertex = truth.conv_vtx[leg.vtx_idx];

		FullTrajectoryInfo conv_el_tr;
		conv_el_tr.is_conversion_track = true;
		conv_el_tr.fPDGCharge = aTrack->GetDynamicParticle()->GetCharge();
		conv_el_tr.fMomentumDir = aTrack->GetDynamicParticle()->GetMomentumDirection();
		conv_el_tr.fEnergy = aTrack->GetDynamicParticle()->GetTotalEnergy();
		conv_el_tr.fMass = aTrack->GetDynamicParticle()->GetMass();
		conv_el_tr.fTrackID        = aTrack->GetTrackID();
		conv_el_tr.fParentID       = vertex.in.ancestor_idx;
		conv_el_tr.fPDGCode        = 22; //Choose to label conv. electron track as a photon //aTrack->GetDefinition()->GetPDGEncoding();
		conv_el_tr.fMomentum       = aTrack->GetMomentum();
		conv_el_tr.caloExtrapolMaxEkin = 0.0;
		conv_el_tr.caloExtrapolEta     = conv_el_tr.fMomentum.getEta();
		conv_el_tr.caloExtrapolPhi     = GetPhi( conv_el_tr.fMomentum.x(),conv_el_tr.fMomentum.y() );
		conv_el_tr.idExtrapolMaxEkin = conv_el_tr.caloExtrapolMaxEkin;
		conv_el_tr.idExtrapolEta     = conv_el_tr.caloExtrapolEta;
		conv_el_tr.idExtrapolPhi     = conv_el_tr.caloExtrapolPhi;
		conv_el_tr.fVertexPosition = aTrack->GetVertexPosition();
		conv_el_tr.fGlobalTime     = aTrack->GetGlobalTime();
		conv_el_tr.vTrackMomentumDir.push_back(aTrack->GetMomentum());
		conv_el_tr.vParentID.push_back(aTrack->GetParentID());
		conv_el_tr.vTrackID.push_back(aTrack->GetTrackID());
		conv_el_tr.vTrackPos.push_back(aTrack->GetPosition());
		conv_el_tr.vTrackTime.push_back(aTrack->GetGlobalTime());
		conv_el_tr.vTrackPDGID.push_back(aTrack->GetDefinition()->GetPDGEncoding());

		trajectories.fAllConvElectrons.push_back( conv_el_tr );
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
	if (!aTrack || aTrack->GetParentID() != 0)
		return;

	// End state of a Geant4 primary (particle_end_*, particle_x0_tracker, ...).
	Full_trajectory_info_data &trajectories = Full_trajectory_info_data::GetInstance();
	for (FullTrajectoryInfo &primary : trajectories.fAllTrajectoryInfo)
	{
		if (primary.fTrackID != aTrack->GetTrackID())
			continue;
		Truth_record_data::GetInstance().EndPrimaryTrack(aTrack, primary);
		break;
	}
}

#endif // __H02TRACKINGACTION_H__
