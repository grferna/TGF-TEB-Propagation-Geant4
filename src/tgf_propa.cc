////////////////////////////////////////////////////////////////////////////////

// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
//
// //
// // ********************************************************************
// // * License and Disclaimer                                           *
// // *                                                                  *
// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
// // * conditions of the Geant4 Software License,  included in the file *
// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
// // * include a list of copyright holders.                             *
// // *                                                                  *
// // * Neither the authors of this software system, nor their employing *
// // * institutes,nor the agencies providing financial support for this *
// // * work  make  any representation or  warranty, express or implied, *
// // * regarding  this  software system or assume any liability for its *
// // * use.  Please see the license in the file  LICENSE  and URL above *
// // * for the full disclaimer and the limitation of liability.         *
// // *                                                                  *
// // * This  code  implementation is the result of  the  scientific and *
// // * technical work of the GEANT4 collaboration.                      *
// // * By using,  copying,  modifying or  distributing the software (or *
// // * any work based  on the software)  you  agree  to acknowledge its *
// // * use  in  resulting  scientific  publications,  and indicate your *
// // * acceptance of all terms of the Geant4 Software license.          *
// // ********************************************************************
////////////////////////////////////////////////////////////////////////////////

#include <DetectorConstruction.hh>
#include <PhysicsList.hh>
#include "G4UImanager.hh"

#ifdef G4VIS_USE

# include "G4VisExecutive.hh"

#endif // ifdef G4VIS_USE

#ifdef G4UI_USE

# include "G4UIExecutive.hh"

#endif // ifdef G4UI_USE

#include <RunAction.hh>

#include <EventAction.hh>
#include <SteppingAction.hh>

#include "Bline_tracer/G4BlineTracer.hh"

#ifndef _WIN32
#include <chrono>
#endif

using namespace std;

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	Settings *settings = Settings::getInstance();

	// random seed, different each time code is run
	// std::chrono::high_resolution_clock m_clock;

#ifdef _WIN32
	long start = time(0);
#else
	std::chrono::high_resolution_clock m_clock;
	long start = (long)std::abs(std::chrono::duration_cast<std::chrono::nanoseconds>(m_clock.now().time_since_epoch()).count());
#endif

	G4cout << start << " ns" << G4endl;

	settings->RANDOM_SEED = start;

	G4String nb_to_get_per_run = "1000";

	settings->record_altitude.push_back(400.); // default value

	G4String Mode = "run";

	G4String number_st = "1000";

	// input parameters
	//     std::cout << "Argument list" << std::endl;
	//     for (int i = 0; i < argc; ++i) {
	//         std::cout << argv[i] << std::endl;
	//     }

	if (argc == 8)
	{
		Mode = "run";
		nb_to_get_per_run = argv[1];
		settings->SOURCE_ALT = std::stod(argv[2]);
		settings->OPENING_ANGLE = std::stod(argv[3]);
		settings->TILT_ANGLE = std::stod(argv[4]);
		settings->BEAMING_TYPE = argv[5];
		settings->SOURCE_SIGMA_TIME = std::stod(argv[6]);

		settings->record_altitude.clear();
		settings->record_altitude.push_back(std::stod(argv[7]));
	}
	else if (argc == 1)
	{
		// default values can be seen in src/src/Settings.cc
		Mode = "visu";
		nb_to_get_per_run = "1000";
	}
	else {
		G4cout << "Wrong number of Arguments passed (should be none, or 7). Aborting" << G4endl;
		std::abort;
	}

	// choose the Random engine and give seed
	//    G4Random::setTheEngine(new CLHEP::MTwistEngine);
	G4Random::setTheSeed(settings->RANDOM_SEED);

	// Construct the default run manager
	auto *runManager = new G4RunManager;

	// Construct the helper class to manage the electric field &
	// the parameters for the propagation of particles in it.

	// set mandatory initialization classes
	runManager->SetUserInitialization(new TGFDetectorConstruction);

	runManager->SetUserInitialization(new TGF_PhysicsList); // for using defined physics list

	// for using Reference physics list

	//   G4PhysListFactory*physListFactory= new G4PhysListFactory();
	// G4VUserPhysicsList *physicsList=physListFactory->GetReferencePhysList("LBE_LIV"); // low and high energy physics with livermore
	//                                                               // _PEN for PENELOPE
	//  runManager->SetUserInitialization(physicsList);
	// std::vector<G4String> v =physListFactory->AvailablePhysLists();

	// set mandatory user action class
	runManager->SetUserAction(new PrimaryGeneratorAction);
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new EventAction);
	runManager->SetUserAction(new SteppingAction);

	//    G4BlineTracer* theBlineTool = new G4BlineTracer();

	Analysis *analysis = Analysis::getInstance();

	// Initialize G4 kernel
	runManager->Initialize();
	G4cout << G4endl << "Initialization OK" << G4endl;

	// get the pointer to the User Interface manager
	G4UImanager *UImanager = G4UImanager::GetUIpointer();

	if (Mode == "visu")
	{
#ifdef G4VIS_USE
		G4VisManager *visManager = new G4VisExecutive;
		visManager->Initialize();
#endif // ifdef G4VIS_USE

#ifdef G4UI_USE
		G4UIExecutive *ui = new G4UIExecutive(argc, argv);
# ifdef G4VIS_USE
		UImanager->ApplyCommand("/control/execute vis.mac");
# endif // ifdef G4VIS_USE
		ui->SessionStart();
		delete ui;
#endif // ifdef G4UI_USE
#ifdef G4VIS_USE
		delete visManager;
#endif // ifdef G4VIS_USE
	}
	else if (Mode == "run")
	{
		while (std::stoi(nb_to_get_per_run) > analysis->get_NB_RECORDED())
			UImanager->ApplyCommand("/run/beamOn " + number_st);
	}

#ifdef _WIN32
	system("pause");
#endif

	delete runManager;

	return 0;
} // main
