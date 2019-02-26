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
// implementation of the IGRF 12 magnetic field that is faster than EarthMagField.hh/cc
// because it uses ECEF cartesian coordinates all the time, so does not require multiple coordinate conversions
////////////////////////////////////////////////////////////////////////////////
#pragma once

#include <vector>

#include "globals.hh"
#include "G4MagneticField.hh"
#include <vector>
#include <geodetic_converter.hh>
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

class EarthMagField_alt : public G4MagneticField {
public:

    EarthMagField_alt();

    ~EarthMagField_alt() override;

    void
    GetFieldValue(const double Point[3], double *Bfield) const override;

private:

    mutable int IENTY = 2;
    mutable float date = 2018.0;
    mutable double alt = 0;
    mutable double f = 0, lat = 0, lon = 0;

    mutable float Bfield_ecef_x = 0;
    mutable float Bfield_ecef_y = 0;
    mutable float Bfield_ecef_z = 0;

    const double earthradius = 6371.2*kilometer;
    mutable float Bfield_mag;
    mutable float XCORD;
    mutable float YCORD;
    mutable float ZCORD;

};
