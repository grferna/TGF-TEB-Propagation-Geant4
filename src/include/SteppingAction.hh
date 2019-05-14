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

#pragma once

#include <Settings.hh>
#include <vector>

#include <Analysis.hh>
#include "G4UserSteppingAction.hh"

class G4Step;

struct record_coords
{
	G4ThreeVector position;
	G4double time;
};

class SteppingAction : public G4UserSteppingAction
{
public:

	//        G4double Get_Altitude(const G4double &x,
	//                              const G4double &y,
	//                              const G4double &z);

	//        G4double Get_dist_rad(const G4double &lat,
	//                              const G4double &lon,
	//                              const G4double &alt_tmp);

	SteppingAction();


	~SteppingAction() override;

	void UserSteppingAction(const G4Step *aStep) override;

private:

	Settings *settings = Settings::getInstance();

	G4int current_NB_EVENT = -10; // just initilization

	std::vector<int> ID_list;

	bool IDpart_not_recorded_yet_elec(G4int ID);

	void analyze_number_electrons_per_primary(const G4Step *step);

	bool is_new_event();

	G4int nb_produced_electrons_from_primary_photons = 0;

	void output_produced_electrons_spec_time(const G4Step *step);

	void analyze_produced_photons();

	G4String asciiFileName = "initilization_value";
	G4String asciiFileName_phot = "initilization_value";

	void output_produced_bremstrahlung_photons_spec_time(const G4Step *step);

	const G4int PDG_positron = -11;
	const G4int PDG_electron = 11;
	const G4int PDG_photon = 22;
};
