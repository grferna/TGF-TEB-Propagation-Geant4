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

#include "globals.hh"
#include "G4ios.hh"

#include "geodetic_converter.hh"
#include "Settings.hh"
// following a singleton pattern

class G4Track;

class TGFAnalysis {
public:

    ~TGFAnalysis();

    void
    save_in_output_buffer(const G4int PDG_NB, const G4double &time, const G4double &energy,
                          const G4double &dist_rad, const G4int ID, const G4double &ecef_x, const G4double &ecef_y,
                          const G4double &ecef_z);

    static TGFAnalysis *
    getInstance();

    void
    write_in_output_file();

    void
    write_in_output_file_endOfRun();

    G4int
    get_NB_RECORDED() const;

private:

    TGFAnalysis();

    static TGFAnalysis *instance;

    std::vector<G4String> output_lines;

    G4String asciiFileName2;

    const unsigned int output_buffer_size = 1;

    G4int NB_RECORDED = 0;

    G4double lat = 0, lon = 0, alt = 0;
};
