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

#include <Settings.hh>
#include <Analysis.hh>
#include <fstream>
#include <iomanip>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"

// singleton pattern initialization
TGFAnalysis *TGFAnalysis::instance = nullptr;

// constructor
TGFAnalysis::TGFAnalysis() {

    const G4double ALT_MAX_RECORDED = *std::max_element(Settings::record_altitudes.begin(),
                                                        Settings::record_altitudes.end());

    const G4String output_filename_second_part = std::to_string(Settings::RANDOM_SEED) + "_"
                                                 + std::to_string(int(ALT_MAX_RECORDED)) + "_" +
                                                 std::to_string(int(Settings::SOURCE_ALT)) + "_"
                                                 + std::to_string(int(Settings::OPENING_ANGLE)) + "_" +
                                                 Settings::BEAMING_TYPE + "_"
                                                 + std::to_string(int(Settings::SOURCE_SIGMA_TIME)) + ".out";

    asciiFileName2 = "./output/detParticles_" + output_filename_second_part;

    std::ofstream asciiFile00(asciiFileName2, std::ios::trunc); // to clean the output file
    asciiFile00.close();

    output_lines.clear();
}

G4int TGFAnalysis::get_NB_RECORDED() const {
    return NB_RECORDED;
}

TGFAnalysis::~TGFAnalysis() = default;

TGFAnalysis *
TGFAnalysis::getInstance() {
    if (instance == 0)
        instance = new TGFAnalysis;

    return instance;
}

void TGFAnalysis::save_in_output_buffer(const G4int PDG_NB, const G4double &time, const G4double &energy,
                                        const G4double &dist_rad, const G4int ID, const G4double &ecef_x,
                                        const G4double &ecef_y,
                                        const G4double &ecef_z) {
    // Write to buffer; only gamma can get here

    geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x * 1000.0, ecef_y * 1000.0, ecef_z * 1000.0, lat, lon,
                                                alt);
    alt = alt / 1000.0;

    std::stringstream buffer;
    buffer << std::scientific << std::setprecision(7); // scientific notation with
    // 5 significant digits
    //   asciiFile << name;
    //   asciiFile << ' ';
    buffer << Settings::RANDOM_SEED;
    buffer << ' ';
    buffer << Settings::SOURCE_ALT;
    buffer << ' ';
    buffer << Settings::OPENING_ANGLE;
    buffer << ' ';
    buffer << Settings::TILT_ANGLE;
    buffer << ' ';
    buffer << Settings::NB_EVENT;
    buffer << ' ';
    buffer << ID;
    buffer << ' ';
    buffer << PDG_NB;
    buffer << ' ';
    buffer << time;
    buffer << ' ';
    buffer << energy;
    buffer << ' ';
    buffer << alt;
    buffer << ' ';
    buffer << lat;
    buffer << ' ';
    buffer << lon;
    buffer << ' ';

    NB_RECORDED++;

    if (Settings::OUTPUT_RadDist) {
        buffer << dist_rad;
    }

    if (Settings::OUTPUT_ECEF_COORDS) {
        buffer << ecef_x;
        buffer << ' ';
        buffer << ecef_y;
        buffer << ' ';
        buffer << ecef_z;
    }

    //    buffer << ' ';
    //    buffer << creator_process;
    buffer << G4endl;

    output_lines.push_back(buffer.str());

    write_in_output_file();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void TGFAnalysis::write_in_output_file() {
    if (output_lines.size() <= output_buffer_size)
        return;

    std::ofstream asciiFile2;
    asciiFile2.open(asciiFileName2, std::ios::app);

    if (asciiFile2.is_open()) {
        for (G4String &line : output_lines) {
            asciiFile2 << line;
        }

        asciiFile2.close();
        output_lines.clear();
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void TGFAnalysis::write_in_output_file_endOfRun() {
    if (output_lines.size() <= 0)
        return;

    std::ofstream asciiFile1;
    asciiFile1.open(asciiFileName2, std::ios::app);

    if (asciiFile1.is_open()) {
        for (G4String &line : output_lines) {
            asciiFile1 << line;
        }

        asciiFile1.close();
        output_lines.clear();
    }
}
