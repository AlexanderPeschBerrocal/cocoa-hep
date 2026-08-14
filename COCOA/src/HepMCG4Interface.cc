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
/// \file eventgenerator/HepMC/HepMCEx01/src/HepMCG4Interface.cc
/// \brief Implementation of the HepMCG4Interface class
//
//
// #include "OutputRunAction.hh"

#include "HepMCG4Interface.hh"
#include "OutputRunAction.hh"

#include "G4RunManager.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <vector>


HepMCG4Interface::HepMCG4Interface()
	: hepmcEvent(0)
{
	;
}

HepMCG4Interface::~HepMCG4Interface()
{
	delete hepmcEvent;
}


G4bool HepMCG4Interface::CheckVertexInsideWorld(const G4ThreeVector &pos) const
{
	G4Navigator *navigator = G4TransportationManager::GetTransportationManager()
								 ->GetNavigatorForTracking();

	G4VPhysicalVolume *world = navigator->GetWorldVolume();
	G4VSolid *solid = world->GetLogicalVolume()->GetSolid();
	EInside qinside = solid->Inside(pos);

	if (qinside != kInside)
		return false;
	else
		return true;
}

void HepMCG4Interface::HepMC2G4(const HepMC::GenEvent *hepmcevt,
									G4Event *g4event)
{
	float eta_primary=-100;
	float phi_primary=-100;
	const double momentum_to_mev = HepMC::Units::conversion_factor(
		hepmcevt->momentum_unit(), HepMC::Units::MEV);
	const double length_to_mm = HepMC::Units::conversion_factor(
		hepmcevt->length_unit(), HepMC::Units::MM);
	std::vector<HepMC::GenParticle *> selected_particles;

	//* loop for vertex
	for (HepMC::GenEvent::vertex_const_iterator vitr = hepmcevt->vertices_begin();
		vitr != hepmcevt->vertices_end(); ++vitr)
	{

		//* is vertex real ?
		G4bool qvtx = false;
		for (HepMC::GenVertex::particle_iterator
				pitr = (*vitr)->particles_begin(HepMC::children);
				pitr != (*vitr)->particles_end(HepMC::children); ++pitr)
		{
			if (!(*pitr)->end_vertex() && (*pitr)->status() == 1)
			{
				qvtx = true;
				break;
			}
		}
		if (!qvtx)
		{
			for (HepMC::GenVertex::particle_iterator
				vpitr= (*vitr)->particles_begin(HepMC::children);
				vpitr != (*vitr)->particles_end(HepMC::children); ++vpitr)
			{ 
				if((*vpitr)->status()  > 20 && (*vpitr)->status() < 30)
				{
					eta_primary =  ( (*vpitr)->momentum().eta() );
					phi_primary =  ( (*vpitr)->momentum().phi() );
				}
			}
            continue;
        }

		for (HepMC::GenVertex::particle_iterator
				vpitr = (*vitr)->particles_begin(HepMC::children);
				vpitr != (*vitr)->particles_end(HepMC::children); ++vpitr)
		{

			if ((*vpitr)->status() != 1)
				continue;
			if ((*vpitr)->momentum().perp() * momentum_to_mev <
				config_json_var.fiducial_cuts.pt_min * GeV)
			{
				continue;
			}
			if (fabs((*vpitr)->momentum().eta()) > config_json_var.fiducial_cuts.eta_max)
			{
				continue;
			}

			// Cut away everything separated from the primary eta,phi by more than dR (used for single-jet data)
			if(config_json_var.fiducial_cuts.dR_cut > 0)
			{
				float eta_cut = (*vpitr)->momentum().eta();
				float phi_cut = (*vpitr)->momentum().phi();
				float dphi = acos(cos(phi_cut-phi_primary));
				float deta = eta_cut-eta_primary;
				float dR = sqrt(dphi*dphi+deta*deta);
				if(dR > config_json_var.fiducial_cuts.dR_cut)
					continue;
			}

			// if the particle decayed "too far" into the detector, replace it with its parent. otherwise this function just returns the original particle.
			HepMC::GenParticle *pptr = m_truthrecordgraph.check_prod_location(*vpitr);

			if (config_json_var.Save_truth_particle_graph)
			{
				m_truthrecordgraph.add_to_vector(pptr, m_truthrecordgraph.m_interesting_particles);
				m_truthrecordgraph.add_all_moving_parents(
					pptr, m_truthrecordgraph.m_interesting_particles);
			}

			G4int pdgcode = (pptr)->pdg_id();
			// skip neutrinos
			if (abs(pdgcode) == 12 || abs(pdgcode) == 14 || abs(pdgcode) == 16)
			{
				continue;
			}

			//add to selected particles (selected to be passed to GEANT), avoid duplicates.
			m_truthrecordgraph.add_to_vector(pptr, selected_particles);
		}
	}

	for (auto *particle : selected_particles)
	{
		const HepMC::GenVertex *production_vertex = particle->production_vertex();
		if (!production_vertex)
			continue;

		const HepMC::FourVector position = production_vertex->position();
		const G4ThreeVector g4_position(position.x() * length_to_mm * mm,
			position.y() * length_to_mm * mm,
			position.z() * length_to_mm * mm);
		if (!CheckVertexInsideWorld(g4_position))
			continue;

		auto *g4_vertex = new G4PrimaryVertex(
			g4_position.x(), g4_position.y(), g4_position.z(),
			position.t() * length_to_mm * mm / c_light);
		const HepMC::FourVector momentum = particle->momentum();
		auto *g4_particle = new G4PrimaryParticle(
			particle->pdg_id(),
			momentum.px() * momentum_to_mev * MeV,
			momentum.py() * momentum_to_mev * MeV,
			momentum.pz() * momentum_to_mev * MeV);
		g4_vertex->SetPrimary(g4_particle);
		g4event->AddPrimaryVertex(g4_vertex);
		if (config_json_var.Save_truth_particle_graph)
			m_truthrecordgraph.m_final_state_particles.push_back(particle);
	}

	if (config_json_var.Save_truth_particle_graph)
		m_truthrecordgraph.fill_truth_graph();
}

HepMC::GenEvent *HepMCG4Interface::GenerateHepMCEvent()
{
	HepMC::GenEvent *aevent = new HepMC::GenEvent();
	return aevent;
}

void HepMCG4Interface::GeneratePrimaryVertex(G4Event *anEvent)
{
	// The truth graph stores non-owning HepMC pointers, so clear them before
	// deleting the event that owns the particles.
	m_truthrecordgraph.clear();
	// delete previous event object
	delete hepmcEvent;
	hepmcEvent = nullptr;

	// generate next event
	hepmcEvent = GenerateHepMCEvent();
	if (!hepmcEvent)
	{

		G4RunManager::GetRunManager()->AbortRun();
		return;
	}
	HepMC2G4(hepmcEvent, anEvent);
}
