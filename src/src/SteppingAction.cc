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

#include <RunAction.hh>
#include <SteppingAction.hh>
#include "G4Track.hh"
#include "G4RunManager.hh"

SteppingAction::SteppingAction()
{

    //    asciiFileName = "./extra_outputs/created_electrons";
    //    std::ofstream asciiFile00(asciiFileName, std::ios::trunc); // to clean the output file
    //    asciiFile00.close();

    //    asciiFileName_phot = "./extra_outputs/created_photons";
    //    std::ofstream asciiFile11(asciiFileName_phot, std::ios::trunc); // to clean the output file
    //    asciiFile11.close();

}

SteppingAction::~SteppingAction() = default;

void SteppingAction::UserSteppingAction(const G4Step *aStep)
{
    // Optional fucntions can be used to get information on particles, e.g.
    // - number of electrons created by primary particle
    // - produced electron energy and time distribution,
    // - bremsstrahlung photon energy and time distribution

    //    if (is_new_event()) {
    ////      G4cout << double(nb_produced_electrons_from_primary_photons) / double(current_NB_EVENT) << G4endl;
    //      ID_list.clear();
    ////      nb_produced_electrons_from_primary_photons = 0;
    //    }

    //    output_produced_electrons_spec_time(aStep);

    //    output_produced_bremstrahlung_photons_spec_time(aStep);

    //    analyze_produced_electrons(aStep);
    G4Track *track = aStep->GetTrack();
    G4double global_time = track->GetGlobalTime();

    if (global_time > settings->TIME_LIMIT)
    {
        track->SetTrackStatus(fStopAndKill);
    }

    if (track->GetKineticEnergy() < settings->MIN_ENERGY_OUTPUT)
    {
        track->SetTrackStatus(fStopAndKill);
    }

}

bool SteppingAction::is_new_event()
{

    if (settings->NB_EVENT != current_NB_EVENT)
    {

        current_NB_EVENT = settings->NB_EVENT;
        return true;
    }

    return false;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::analyze_number_electrons_per_primary(const G4Step *step)
{
    G4Track *track = step->GetTrack();

    G4String creator_process = "initial";

    if (track->GetCreatorProcess())
    {
        creator_process = track->GetCreatorProcess()->GetProcessName();
    }

    if (track->GetParentID())
    {
        G4int parent_ID = track->GetParentID();
        G4int ID = track->GetTrackID();

        if (parent_ID == 1) // if primary particle (necessarily photon...)
        {
            {
                //            if (track->GetParticleDefinition()->GetParticleName()=="e-")
                //            {

                if (IDpart_not_recorded_yet_elec(ID))
                {
                    ID_list.push_back(ID);
                    nb_produced_electrons_from_primary_photons++;
                }

                //            }
            }
        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::output_produced_electrons_spec_time(const G4Step *step)
{
    G4Track *track = step->GetTrack();

    if (track->GetTrackID())
    {
        G4int ID = track->GetTrackID();

        if (IDpart_not_recorded_yet_elec(ID) && ID != 1)
        {

            ID_list.push_back(ID);

            if (track->GetParticleDefinition()->GetParticleName() == "e-")
            {
                if (track->GetKineticEnergy() > 10.0 * keV)
                {

                    std::ofstream asciiFile00;
                    asciiFile00.open(asciiFileName, std::ios::app);

                    if (asciiFile00.is_open())
                    {
                        asciiFile00 << track->GetKineticEnergy() / keV << "  " << track->GetGlobalTime() / microsecond << G4endl;
                        asciiFile00.close();
                    }

                }
            }

        }
    }
}

void SteppingAction::output_produced_bremstrahlung_photons_spec_time(const G4Step *step)
{
    G4Track *track = step->GetTrack();

    if (track->GetTrackID())
    {
        G4int ID = track->GetTrackID();

        if (IDpart_not_recorded_yet_elec(ID) && ID != 1)
        {

            ID_list.push_back(ID);

            if (track->GetParticleDefinition()->GetParticleName() == "gamma")
            {
                G4String creator_process = track->GetCreatorProcess()->GetProcessName();

                //                    G4cout << creator_process << G4endl;
                if (creator_process == "eBrem") // check if it is  produced by bremsstrahlung (positron annihiliation is also possible but not wanted to record)
                {

                    if (track->GetKineticEnergy() > 10.0 * keV)
                    {

                        std::ofstream asciiFile11;
                        asciiFile11.open(asciiFileName_phot, std::ios::app);

                        if (asciiFile11.is_open())
                        {
                            asciiFile11 << track->GetKineticEnergy() / keV << "  " << track->GetGlobalTime() / microsecond << G4endl;
                            asciiFile11.close();
                        }
                    }
                }
            }

        }
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::analyze_produced_photons()
{

}

bool SteppingAction::IDpart_not_recorded_yet_elec(G4int ID)
{
    if (ID_list.empty())
    {
        return true;
    }

    return !(std::find(ID_list.begin(), ID_list.end(), ID) != ID_list.end());
}
