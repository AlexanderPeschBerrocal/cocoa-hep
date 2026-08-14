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
/// \file eventgenerator/HepMC/HepMCEx03/src/HepMCG4Pythia8Interface.cc
/// \brief Implementation of the HepMCG4PythiaInterface class for Pythia8
//

#ifdef G4LIB_USE_PYTHIA8
// ! Defenetly active when we use Pythia
#include "HepMCG4Pythia8Interface.hh"
#include "HepMCG4Pythia8Messenger.hh"
// #include "OutputRunAction.hh"
#include "HepMC/GenEvent.h"
#include "Randomize.hh"
#include <memory>
// #include "Pythia8/Pythia.h"

using namespace Pythia8;
//using namespace std;

// ================================= //

bool MyUserHooks::canVetoProcessLevel()
{
	return true;
}

bool MyUserHooks::canModifySigma()
{
	return true;
}

double MyUserHooks::multiplySigmaBy(const Pythia8::SigmaProcess *sigmaProcessPtr, const Pythia8::PhaseSpace *phaseSpacePtr, bool /*inEvent*/)
{
	// All events should be 2 -> 2, but kill them if not.
	// if (sigmaProcessPtr->nFinal() != 2) return 0.;
	double weight = 0;
	double Hmass = sigmaProcessPtr->m(3);
	if (Hmass > 130 or Hmass < 120)
		weight = 1;
	else if (sqrt(phaseSpacePtr->sHat()) > 800)
		weight = 1;
	return weight;
}

bool MyUserHooks::doVetoProcessLevel(Pythia8::Event &process)
{
	bool pass = false;
	int numL = 0;
	for (int i = 0; i < process.size(); ++i)
	{
		if (process.at(i).status() > 22 && process.at(i).pT() > 1)
		{
			if (fabs(process.at(i).eta()) > 3.0)
				pass = true; //3.0
			if (fabs(process.at(process.at(i).mother1()).id()) == 24 && (fabs(process.at(i).id()) == 11 || fabs(process.at(i).id()) == 13))
				numL++;
			if (fabs(process.at(process.at(i).mother1()).id()) == 24 && fabs(process.at(i).id()) == 15)
				pass = true;
		}
	}
	if (numL != 1)
		pass = true;
	return pass;
}
MyUserHooks::MyUserHooks()
{
	slowJet = new Pythia8::SlowJet(-1, 0.7, 10., 5.);
	eff = 0;
	p = 0;
	f = 0;
}

MyUserHooks::~MyUserHooks()
{
	delete slowJet;
}

HepMCG4Pythia8Interface::HepMCG4Pythia8Interface() : verbose(0)
{
	messenger = new HepMCG4Pythia8Messenger(this);
}

HepMCG4Pythia8Interface::~HepMCG4Pythia8Interface()
{
	delete messenger;
}

bool HepMCG4Pythia8Interface::CallPythiaReadString(const G4String &par)
{
	return pythia.readString(par);
}

bool HepMCG4Pythia8Interface::CallPythiaInit()
{

// 	auto myUserHooks = make_shared<MyUserHooks>();
// 	pythia.setUserHooksPtr(myUserHooks);
	//MyUserHooks myUserHooks = MyUserHooks();
	//pythia.setUserHooksPtr(myUserHooks);
	return pythia.init();
}

bool HepMCG4Pythia8Interface::ConfigureBeams(G4int beamPdg,
		G4int targetPdg, G4double eCMGeV)
{
	return pythia.readString("Beams:idA = " + std::to_string(beamPdg)) &&
		pythia.readString("Beams:idB = " + std::to_string(targetPdg)) &&
		pythia.readString("Beams:eCM = " + std::to_string(eCMGeV)) &&
		pythia.init();
}

void HepMCG4Pythia8Interface::CallPythiaStat()
{
	pythia.stat();
}

Pythia8::Event HepMCG4Pythia8Interface::GetPythiaObject()
{
	return sum_events;
}

bool HepMCG4Pythia8Interface::SetRandomSeed(G4int iseed)
{
	if (iseed < 1 || iseed > 900000000)
		return false;
	if (!pythia.readString("Random:setSeed = on"))
		return false;
	ostringstream Seed;
	Seed << "Random:seed = " << iseed;
	return pythia.readString(Seed.str());
}

void HepMCG4Pythia8Interface::PrintRandomStatus(std::ostream &ostr) const //
{
	ostr << "RandomStatus: " << G4endl;
}
void HepMCG4Pythia8Interface::PrintRandomStatus() const //
{
	;
}

void HepMCG4Pythia8Interface::SetUserParameters()
{
	G4cout << "set user parameters of PYTHIA common." << G4endl
		   << "nothing to be done in default."
		   << G4endl;
}

HepMC::GenEvent *HepMCG4Pythia8Interface::GenerateHepMCEvent()
{
	//* Pile-up
	// double nPileupAvg = 2.5;
	// int nPileup = 0;//poisson(nPileupAvg, pythia.rndm);

	// HepMC::GenEvent* hepmcevt;
	// hepmcevt = new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM);
	// //signall event
	// pythia.next();
	// sum_events = pythia.event;
	// ToHepMC.fill_next_event( pythia, hepmcevt,-1, false );
	// if(verbose>0)
	//    {./Sc
	//       hepmcevt-> print();
	//    }
	// //pile-up
	// for (int iPileup = 0; iPileup < nPileup; ++iPileup)
	// {
	//    pythia.next();
	//    sum_events += pythia.event;
	//    ToHepMC.fill_next_event( pythia, hepmcevt,-1, false );
	//    if(verbose>0)
	//    {
	//       hepmcevt-> print();
	//    }

	// }

	// return hepmcevt;
	//*  Pile-up end
	int id = messenger->GQparticle;
	std::unique_ptr<HepMC::GenEvent> hepmcevt(
		new HepMC::GenEvent(HepMC::Units::MEV, HepMC::Units::MM));
	if (id == 22122212) //* pp colision
	{
		constexpr int maxAttempts = 10;
		bool generated = false;
		for (int attempt = 0; attempt < maxAttempts && !generated; ++attempt)
			generated = pythia.next();
		if (!generated)
		{
			G4cerr << "Pythia8 event generation failed after " << maxAttempts
				   << " attempts" << G4endl;
			return nullptr;
		}

		ToHepMC.fill_next_event(pythia, hepmcevt.get());

		if (verbose > 1)
			hepmcevt->print();

		sum_events = pythia.event;
	}
	else //*  Parton or particle with colour singlet
	{
		Pythia8::Event &event = pythia.event;
		if (messenger->MinEnergy <= 0. ||
			messenger->MaxEnergy < messenger->MinEnergy ||
			messenger->MaxEta < messenger->MinEta)
		{
			G4cerr << "Invalid particle-gun energy or eta range" << G4endl;
			return nullptr;
		}

		constexpr int maxAttempts = 5;
		bool generated = false;
		for (int attempt = 0; attempt < maxAttempts && !generated; ++attempt)
		{
			if (!pythia.readString("ProcessLevel:all = off") ||
				!pythia.readString("Next:numberShowInfo = 0") ||
				!pythia.readString("Next:numberShowProcess = 0") ||
				!pythia.readString("Next:numberShowEvent = 0") ||
				!SetRandomSeed(1 + CLHEP::RandFlat::shootInt(899999999)) ||
				!pythia.init())
			{
				G4cerr << "Pythia8 initialization failed" << G4endl;
				return nullptr;
			}

			double minenergy = messenger->MinEnergy;
			double maxenergy = messenger->MaxEnergy;
			double mineta = messenger->MinEta;
			double maxeta = messenger->MaxEta;

			double ee = CLHEP::RandFlat::shoot(minenergy, maxenergy); //100;

			double eta = CLHEP::RandFlat::shoot(mineta, maxeta);

			double phi = CLHEP::RandFlat::shoot(0., 2 * M_PI);
			Pythia8::ParticleData &pdt = pythia.particleData;
			// Reset event record to allow for new event.
			event.reset();
			if (id == 21) //* Gluon
			{
				double pt = ee / cosh(eta);
				double pz = pt * sinh(eta);
				double px = pt * cos(phi);
				double py = pt * sin(phi);
				event.append(id, 23, 101, 102, px, py, pz, ee);
				event.append(id, 23, 102, 101, -px, -py, -pz, ee);
			}
			else if (id == 2) //*  Quark
			{
				double mass_0 = pdt.m0(id);
				double pp = Pythia8::sqrtpos(ee * ee - mass_0 * mass_0);
				double pt = pp / cosh(eta);
				double pz = pt * sinh(eta);
				double px = pt * cos(phi);
				double py = pt * sin(phi);
				event.append(id, 23, 101, 0, px, py, pz, ee, mass_0);
				event.append(-id, 23, 0, 101, -px, -py, -pz, ee, mass_0);
			}
			else if (id == 22 || id == 11 || id == 130 || id == 13 || id == 211) //* Photon, Electron, K0long, Muon, Charged Pion
			{
				fillParticle(id, ee, eta, phi, 23, event, pdt, false);
			}
			else //* Any particle with colour singlet
			{
				int sign = 0;
				if (CLHEP::RandFlat::shoot(0.0, 1.0) > 0.5) sign = 1;
				else sign = -1;
				eta = sign * eta;
				if (config_var.Type_of_running == "TopoClustering_debug" || config_var.Type_of_running == "Detector_Response")
				{
					fillParticle(211, ee, eta, phi, 24, event, pdt, false);
					fillParticle(111, ee, eta, phi + M_PI, 23, event, pdt, false);
				}
				else
				{
					fillParticle(id, ee, eta, phi, 2, event, pdt, false);
					// fillParticle(-id, ee, -eta, phi, 2, event, pdt, false);
				}
			}
			generated = pythia.next();
		}
		if (!generated)
		{
			G4cerr << "Pythia8 event generation failed after " << maxAttempts
				   << " attempts" << G4endl;
			return nullptr;
		}
		ToHepMC.fill_next_event(pythia, hepmcevt.get(), -1, false);
		sum_events = event;
		if (verbose > 1)
		{
			hepmcevt->print();
		}
	}
	return hepmcevt.release();
}

void HepMCG4Pythia8Interface::fillParticle(int pdgid, double ee, double eta_init, double phi_init, int status, Pythia8::Event &event, Pythia8::ParticleData &pdt, bool atRest)
{
	// Select particle mass; where relevant according to Breit-Wigner.
	double mm = pdt.m0(pdgid);
	double pp = Pythia8::sqrtpos(ee * ee - mm * mm);

	// Special case when particle is supposed to be at rest.
	if (atRest)
	{
		ee = mm;
		pp = 0.;
	}
	double pt = pp / cosh(eta_init);
	double pz = pt * sinh(eta_init);
	double px = pt * cos(phi_init);
	double py = pt * sin(phi_init);

	// Store the particle in the event record.
	event.append(pdgid, status, 0, 0, px, py, pz, ee, mm);
}

void HepMCG4Pythia8Interface::Print() const
{
	;
}

//* for pile-up
int HepMCG4Pythia8Interface::poisson(double nAvg, Pythia8::Rndm &rndm)
{

	// Set maximum to avoid overflow.
	const int NMAX = 100;

	// Random number.
	double rPoisson = rndm.flat() * exp(nAvg);

	// Initialize.
	double rSum = 0.;
	double rTerm = 1.;

	// Add to sum and check whether done.
	for (int i = 0; i < NMAX;)
	{
		rSum += rTerm;
		if (rSum > rPoisson)
			return i;

		// Evaluate next term.
		++i;
		rTerm *= nAvg / i;
	}
	return NMAX;
}

#endif
