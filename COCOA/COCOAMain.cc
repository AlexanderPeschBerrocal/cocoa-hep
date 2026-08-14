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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string> 
#include <stdexcept>
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4Types.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Config_reader_func.hh"
#include "Config_reader_var.hh"

#include "RunTest.hh"

static void show_usage(std::string name)
{
	std::cerr << "Usage: \n" << name << " <option(s)> "
			  << "Options:\n"
			  << "\t--config (-c) <str>\t path to json configuration file\n"
			  << "\t--macro (-m) <str>\t path to Geant4, Pythia8, or HepMC macro file for event generation (can be set in json configuration file)\n"
			  << "\t--output (-o) <str>\t path (incl. name) of output ROOT file to be written (can be set in json configuration file)\n"
			  << "\t--input (-i) <str>\t path to HepMC (.hmc) input file (overrides the default path set in the HepMC macro file)\n"
			  << "\t--seed (-s) <int>\t set random seed\n"
			  << "\t--nevents (-n) <int>\t number of events to generate (default is taken from macro).\n"
			  << "\t--help (-h)\t show this message\n"
			  << "no <option(s)> will call UI interactive command submission\n" 
			  << std::endl;
}

static bool apply_command(G4UImanager *ui_manager, const G4String &command)
{
	const G4int status = ui_manager->ApplyCommand(command);
	if (status != 0)
	{
		G4cerr << "Geant4 command failed (status " << status << "): "
			   << command << G4endl;
		return false;
	}
	return true;
}


int main(int argc, char **argv)
{
	std::string path_to_config = "./config/config_default.json";
	time_t systime = time(NULL);
	G4long seed = (long)systime;
	G4UIExecutive *ui = nullptr;
	std::string root_file_path = "";
	std::string macro_file_path = "";
	std::string input_file_path = "";
	int nEvents = -1;

	if (argc == 1)
	{
		ui = new G4UIExecutive(argc, argv);
	}
	else if (argc == 2)
	{
		std::string arg = argv[1];
		if (arg == "--help" || arg == "-h")
		{
			show_usage(argv[0]);
			return 0;
		}
		else
		{
			std::cerr << "Option "<< arg <<" does not exist!" << std::endl;
			show_usage(argv[0]);
		}
		return 1;
	}
	else
	{
		for (int i = 1; i < argc; ++i)
		{
			std::string arg = argv[i];
			if (arg == "--output" || arg == "-o")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					root_file_path = argv[i];
				}
				else
				{
					std::cerr << "--output option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--macro" || arg == "-m")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					macro_file_path = argv[i];
				}
				else
				{
					std::cerr << "--macro option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--seed" || arg == "-s")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					try
					{
						seed = std::stol(argv[i]);
					}
					catch (const std::exception &)
					{
						std::cerr << "--seed requires a valid integer." << std::endl;
						return 1;
					}
				}
				else
				{
					std::cerr << "--seed option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--nevents" || arg == "-n")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					try
					{
						nEvents = std::stoi(argv[i]);
					}
					catch (const std::exception &)
					{
						std::cerr << "--nevents requires a valid integer." << std::endl;
						return 1;
					}
					if (nEvents <= 0)
					{
						std::cerr << "--nevents must be positive." << std::endl;
						return 1;
					}
				}
				else
				{
					std::cerr << "--nevents option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--config" || arg == "-c")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					path_to_config = argv[i];
				}
				else
				{
					std::cerr << "--config option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--input" || arg == "-i")
			{
				if (i + 1 < argc) // Make sure we aren't at the end of argv
				{
					i++;
					input_file_path = argv[i];
				}
				else
				{
					std::cerr << "--input option requires one argument." << std::endl;
					return 1;
				}
			}
			else if (arg == "--help" || arg == "-h")
			{
				show_usage(argv[0]);
				return 0;
			}
			else
			{
				std::cerr << "Option "<< arg <<" does not exist!" << std::endl;
				show_usage(argv[0]);
				return 1;
			}
		}
	}
	Config_reader_var &config_var = Config_reader_var::GetInstance();
	try
	{
		Config_reader_func config_json_func(path_to_config, config_var);
	}
	catch (const std::exception &error)
	{
		G4cerr << "Configuration error: " << error.what() << G4endl;
		return 1;
	}
	if (seed <= 0)
	{
		G4cerr << "Seed must be a positive integer" << G4endl;
		return 1;
	}
	
	//* choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
	CLHEP::HepRandom::setTheSeed(seed);
	// }
	Geometry_definition geometry = config_var.low_resolution;
	if (config_var.Use_high_granularity)
	{
		geometry = config_var.high_resolution;
	}
	if (root_file_path.empty())
	{
		root_file_path = config_var.Output_file_path;
		if (root_file_path.empty())
		{
			G4cerr << "root_file_path is not given!" << G4endl;
			return 1;
		}
	}

	G4RunManager *runManager = new G4RunManager;

	// User Initialization classes (mandatory)
	//
	runManager->SetUserInitialization(new DetectorConstruction(geometry));

	//
	runManager->SetUserInitialization(new FTFP_BERT);
	runManager->SetUserInitialization(new ActionInitialization(
		geometry, root_file_path, config_var.Save_truth_particle_graph));

	runManager->Initialize();

	G4VisManager *visManager = new G4VisExecutive;
	visManager->Initialize();

	//get the pointer to the User Interface manager
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	if (!ui)
	{ // batch mode

		if ((macro_file_path == ""))
		{
			if (config_var.Macro_file_path!="")
				macro_file_path = config_var.Macro_file_path;
			else 
			{
				G4cout<<"macro_file_path is not given!"<<G4endl;
				return 1;
			}
		}
		
		visManager->SetVerboseLevel("quiet");
		const G4long pythia_seed = 1 + (seed % 899999999L);
		if (!apply_command(UImanager, "/generator/pythia8/setSeed " +
			std::to_string(pythia_seed)))
			return 1;
		const int number_of_events = nEvents > 0 ? nEvents : config_var.Number_of_events;
		if (!apply_command(UImanager, "/control/alias numberOfEvents " +
			std::to_string(number_of_events)))
			return 1;

		if (input_file_path.empty())
		{
			if (!apply_command(UImanager, G4String("/control/execute ") + macro_file_path))
				return 1;
		}
		else
		{
			std::ifstream filestream(macro_file_path);
			if (!filestream.is_open())
			{
				G4cerr << "Cannot open macro file: " << macro_file_path << G4endl;
				return 1;
			}
			std::string line;
			while (std::getline(filestream, line))
			{
				const std::size_t first = line.find_first_not_of(" \t\r");
				if (first == std::string::npos || line[first] == '#')
					continue;
				G4String command = line.substr(first);
				if (command.find("/generator/hepmcAscii/open") == 0)
					command = "/generator/hepmcAscii/open " + input_file_path;
				if (!apply_command(UImanager, command))
					return 1;
			}
		}
	}
	else
	{ // interactive mode : define UI session
		if (!apply_command(UImanager, "/control/execute init_vis.mac"))
			return 1;
		// if (ui->IsGUI())
		// {
		// 	UImanager->ApplyCommand("/control/execute gui.mac");
		// }

		ui->SessionStart();
		delete ui;
	}

	// Free the store: user actions, physics_list and detector_description are
	//                 owned and deleted by the run manager, so they should not
	//                 be deleted in the main() program !

	////////////////////////////
	// Test output if desired //
	////////////////////////////

	int exit_code = 0;

	test_t test_type = test_t::kNONE;
	if ( config_var.run_hadron_test )
	    test_type = kHADRON;
	if ( config_var.run_piZero_test )
	    test_type = kPI_ZERO;

	if ( test_type )
	    exit_code = RunTest( test_type );
	
	
	delete visManager;
	delete runManager;
	//  delete UImanager;

	return exit_code;
}
