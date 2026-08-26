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
// #include "DataStorage.hh"

#include "SteppingAction.hh"
//#include "FullTrajectoryInfo.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
using namespace std;

#include "G4VProcess.hh"
#include "G4ProcessType.hh"

#include <algorithm>
#include <string>

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4EmProcessSubType.hh"

#include <cstdlib>

#include "G4LogicalVolume.hh"
#include "G4Region.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
TrackingAction::TrackingAction()
	: G4UserTrackingAction()
{
	;
}

void TrackingAction::PreUserTrackingAction(const G4Track*aTrack) 
{
	auto runGeneratorAction = static_cast<const PrimaryGeneratorAction *>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
	G4String GeneratorName = runGeneratorAction->GetGeneratorName();
	Full_trajectory_info_data &trajectories = Full_trajectory_info_data::GetInstance();
	G4int ParentID = aTrack->GetParentID();
	if (aTrack->GetParentID() == 0)
	{
	        FullTrajectoryInfo trjInfo;
		trjInfo.is_conversion_track = false;
		trjInfo.fParentID = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTrackID();
		trjInfo.fTrackID = aTrack->GetTrackID();
		trjInfo.fPDGCharge = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetCharge();
		trjInfo.fPDGCode = aTrack->GetDefinition()->GetPDGEncoding();
		trjInfo.fMomentum = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentum();
		trjInfo.fMomentumDir = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMomentumDirection();
		trjInfo.fEnergy = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
		trjInfo.fMass = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetMass();

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

		trajectories.fAllTrajectoryInfo.push_back(trjInfo);
		//!Pythia8
		if (GeneratorName == "pythia8")
		{
			if (aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetCharge() == 0)
			{
				float Enu = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
				(void) Enu;
			}
			else
			{
				float Ech = aTrack->GetDynamicParticle()->GetPrimaryParticle()->GetTotalEnergy();
				(void) Ech;
			}
			// trajectories.AddTrueEnergy(Ech, Enu);
		}
		//!Pythia8 END
		//break;
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
	if ( IsConversionElectron( aTrack ) &&
		 HasPrimaryPhotonParent( aTrack ) &&
	     IsVertexInTrackerRegion( aTrack ) ) {

	    FullTrajectoryInfo conv_el_tr;
		conv_el_tr.is_conversion_track = true;
		conv_el_tr.fPDGCharge = aTrack->GetDynamicParticle()->GetCharge();
		conv_el_tr.fMomentumDir = aTrack->GetDynamicParticle()->GetMomentumDirection();
		conv_el_tr.fEnergy = aTrack->GetDynamicParticle()->GetTotalEnergy();
		conv_el_tr.fMass = aTrack->GetDynamicParticle()->GetMass();
		conv_el_tr.fTrackID        = aTrack->GetTrackID();
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

	    for( size_t i_primary_tr = 0; i_primary_tr < trajectories.fAllTrajectoryInfo.size(); ++i_primary_tr ) {
		if ( trajectories.fAllTrajectoryInfo[i_primary_tr].fTrackID == aTrack->GetParentID() ) {
		    conv_el_tr.fParentID = i_primary_tr;
		    break;
		}		    
	    }

	    trajectories.fAllConvElectrons.push_back( conv_el_tr );
	}

	if (IsNuclearInteractionDaughter(aTrack) && HasPrimaryParticleParent(aTrack) && IsVertexInTrackerRegion(aTrack)) {
		FullTrajectoryInfo nucl_int_tr;

		nucl_int_tr.is_conversion_track = false;
		nucl_int_tr.fPDGCharge = aTrack->GetDynamicParticle()->GetCharge();
		nucl_int_tr.fMomentumDir = aTrack->GetDynamicParticle()->GetMomentumDirection();
		nucl_int_tr.fEnergy = aTrack->GetDynamicParticle()->GetTotalEnergy();
		nucl_int_tr.fMass = aTrack->GetDynamicParticle()->GetMass();

		nucl_int_tr.fTrackID = aTrack->GetTrackID();
		nucl_int_tr.fPDGCode = aTrack->GetDefinition()->GetPDGEncoding();
		nucl_int_tr.fMomentum = aTrack->GetMomentum();

		nucl_int_tr.caloExtrapolMaxEkin = 0.0;
		nucl_int_tr.caloExtrapolEta = nucl_int_tr.fMomentum.getEta();
		nucl_int_tr.caloExtrapolPhi = GetPhi(nucl_int_tr.fMomentum.x(), nucl_int_tr.fMomentum.y());

		nucl_int_tr.idExtrapolMaxEkin = nucl_int_tr.caloExtrapolMaxEkin;
		nucl_int_tr.idExtrapolEta = nucl_int_tr.caloExtrapolEta;
		nucl_int_tr.idExtrapolPhi = nucl_int_tr.caloExtrapolPhi;

		nucl_int_tr.fVertexPosition = aTrack->GetVertexPosition();
		nucl_int_tr.fGlobalTime = aTrack->GetGlobalTime();

		nucl_int_tr.vTrackMomentumDir.push_back(aTrack->GetMomentum());
		nucl_int_tr.vParentID.push_back(aTrack->GetParentID());
		nucl_int_tr.vTrackID.push_back(aTrack->GetTrackID());
		nucl_int_tr.vTrackPos.push_back(aTrack->GetPosition());
		nucl_int_tr.vTrackTime.push_back(aTrack->GetGlobalTime());
		nucl_int_tr.vTrackPDGID.push_back(aTrack->GetDefinition()->GetPDGEncoding());

		nucl_int_tr.fParentID = FindPrimaryAncestorIndex(aTrack);

		trajectories.fAllNuclearInteractionDaughters.push_back(nucl_int_tr);
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
	if (!aTrack || aTrack->GetParentID() != 0 ||
		aTrack->GetDefinition()->GetPDGEncoding() != 22)
		return;

	auto& primaryParticles =
		Full_trajectory_info_data::GetInstance().fAllTrajectoryInfo;
	for (FullTrajectoryInfo& primaryParticle : primaryParticles) {
		if (primaryParticle.fTrackID != aTrack->GetTrackID() ||
			primaryParticle.fPDGCode != 22 ||
			primaryParticle.is_conversion_track)
			continue;

		primaryParticle.fTrackLength = aTrack->GetTrackLength();
		primaryParticle.fEndPosition = aTrack->GetPosition();
		primaryParticle.fTrackStatus = static_cast<G4int>(aTrack->GetTrackStatus());

		const G4Step* finalStep = aTrack->GetStep();
		const G4VProcess* terminatingProcess = nullptr;
		if (finalStep && finalStep->GetPostStepPoint())
			terminatingProcess =
				finalStep->GetPostStepPoint()->GetProcessDefinedStep();

		if (terminatingProcess) {
			primaryParticle.fTerminationProcessType =
				static_cast<G4int>(terminatingProcess->GetProcessType());
			primaryParticle.fTerminationProcessSubType =
				terminatingProcess->GetProcessSubType();
			primaryParticle.fConverted =
				terminatingProcess->GetProcessType() == fElectromagnetic &&
				terminatingProcess->GetProcessSubType() == fGammaConversion;
		}
		break;
	}
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool TrackingAction::IsPrimaryPhotonDaughter(const G4Track* aTrack) const {

    if ( fabs( aTrack->GetDefinition()->GetPDGEncoding() ) != 11 )
	return false;

    G4int parentID                = aTrack->GetParentID();
    for ( const FullTrajectoryInfo& primaryParticle : Full_trajectory_info_data::GetInstance().fAllTrajectoryInfo ) {
	if ( primaryParticle.fTrackID == parentID &&
	     primaryParticle.fPDGCode == 22 )
	    return true;
    }

    return false;
    
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool TrackingAction::IsInnerDetectorTrack(const G4Track* aTrack) const {

    G4String logicalVolumeName = aTrack->GetVolume()->GetName();
    logicalVolumeName.toLower();
    
    return logicalVolumeName.substr( 0, 5 ) == "inner";
    
}

bool TrackingAction::IsVertexInTrackerRegion(
    const G4Track* aTrack) const
{
    if (!aTrack)
        return false;

    const G4LogicalVolume* vertexVolume =
        aTrack->GetLogicalVolumeAtVertex();

    if (!vertexVolume)
        return false;

    const G4Region* region = vertexVolume->GetRegion();

    return region && region->GetName() == "TrackerRegion";
}

bool TrackingAction::IsConversionElectron(const G4Track* aTrack) const
{
    if (!aTrack || !aTrack->GetDefinition())
        return false;

    // Conventionally includes both the electron and positron.
    if (std::abs(aTrack->GetDefinition()->GetPDGEncoding()) != 11)
        return false;

    const G4VProcess* creator = aTrack->GetCreatorProcess();

    if (!creator)
        return false;

    return creator->GetProcessType() == fElectromagnetic &&
           creator->GetProcessSubType() == fGammaConversion;
}

bool TrackingAction::HasPrimaryPhotonParent(const G4Track* aTrack) const
{
    if (!aTrack || aTrack->GetParentID() == 0)
        return false;

    const G4int parentID = aTrack->GetParentID();
    const auto& primaryParticles =
        Full_trajectory_info_data::GetInstance().fAllTrajectoryInfo;

    for (const FullTrajectoryInfo& primaryParticle : primaryParticles) {
        if (primaryParticle.fTrackID == parentID &&
            primaryParticle.fPDGCode == 22) {
            return true;
        }
    }

    return false;
}

bool TrackingAction::IsNuclearInteractionDaughter(const G4Track* aTrack) const {
    if (!aTrack) return false;

    const G4VProcess* creator = aTrack->GetCreatorProcess();
    if (!creator) return false;

    const G4String& processName = creator->GetProcessName();

    // Require a Geant4 hadronic/nuclear creation process.
    // This rejects decays and ordinary EM processes.
    if (creator->GetProcessType() != fHadronic) return false;

    // Strict inelastic/capture definition.
    // This accepts typical names such as pi+Inelastic, pi-Inelastic,
    // kaon+Inelastic, protonInelastic, neutronInelastic, nCapture, nFission.
    const bool isNuclearProcess =
        processName.find("Inelastic") != std::string::npos ||
        processName == "nCapture" ||
        processName == "nFission";

    return isNuclearProcess;
}

bool TrackingAction::HasPrimaryParticleParent(const G4Track* aTrack) const
{
    if (!aTrack || aTrack->GetParentID() == 0)
        return false;

    const G4int parentID = aTrack->GetParentID();
    const auto& primaries = Full_trajectory_info_data::GetInstance().fAllTrajectoryInfo;

    for (const FullTrajectoryInfo& primary : primaries)
    {
        // A direct match means the incident parent itself was primary.
        if (primary.fTrackID != parentID)
            continue;

        const G4ParticleDefinition* definition =
            G4ParticleTable::GetParticleTable()->FindParticle(primary.fPDGCode);

        if (!definition)
            return false;

        const G4String& type = definition->GetParticleType();

        return type == "baryon" ||
               type == "meson"  ||
               type == "nucleus";
    }

    return false;
}

int TrackingAction::FindPrimaryAncestorIndex(const G4Track* aTrack) const {
    if (!aTrack) return -1;

    const G4int parentID = aTrack->GetParentID();
    const auto& primaries = Full_trajectory_info_data::GetInstance().fAllTrajectoryInfo;

    for (size_t ip = 0; ip < primaries.size(); ++ip) {
        const FullTrajectoryInfo& primary = primaries[ip];

        if (primary.fTrackID == parentID) {
            return static_cast<int>(ip);
        }

        if (std::find(primary.vTrackID.begin(),
                      primary.vTrackID.end(),
                      parentID) != primary.vTrackID.end()) {
            return static_cast<int>(ip);
        }
    }

    return -1;
}

#endif // __H02TRACKINGACTION_H__
