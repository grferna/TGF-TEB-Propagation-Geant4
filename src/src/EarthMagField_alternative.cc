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

#include "PrimaryGeneratorAction.hh"
#include "EarthMagField_alternative.hh"

using namespace std;

// interface avec code fortran igrf12 provenance du site : http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
extern "C" {
void cofrm_(float *date);

void feldg_(int *IENTY, float *XCORD, float *YCORD, float *ZCORD, float *Bfield_ecef_x, float *Bfield_ecef_y, float *Bfield_ecef_z, float *Bfield_mag);
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_alt::EarthMagField_alt()
{

    cofrm_(&date);

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField_alt::~EarthMagField_alt() = default;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField_alt::GetFieldValue(const double Point[3], double *Bfield) const
{
    //  geodetic_converter::GeodeticConverter g_geodetic_converter;

    XCORD = static_cast<float>(Point[0] / earthradius);
    YCORD = static_cast<float>(Point[1] / earthradius);
    ZCORD = static_cast<float>(Point[2] / earthradius);

    feldg_(&IENTY, &XCORD, &YCORD, &ZCORD, &Bfield_ecef_x, &Bfield_ecef_y, &Bfield_ecef_z, &Bfield_mag);


    Bfield[0] = static_cast<double>(Bfield_ecef_x) * gauss;
    Bfield[1] = static_cast<double>(Bfield_ecef_y) * gauss;
    Bfield[2] = static_cast<double>(Bfield_ecef_z) * gauss;

    /////// Simple dipole (very fast)
    //   const G4double B0 = 3.07e-5*tesla;
    //
    //   G4double x =Point[0] /m;
    //   G4double y =Point[1] /m;
    //   G4double z =Point[2] /m;
    //
    //   const G4double Re = earthRadius/m ; // in meters
    //
    //   G4double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
    //
    //   G4double fact = -B0 * pow(Re,3) / pow(r,5);
    //
    //   // Dipole magnetic field
    //
    //   // magnetic field
    //   Bfield[0]=fact*3.*x*z;
    //   Bfield[1]=fact*3.*y*z;
    //   Bfield[2]=fact*(2.*pow(z,2)-pow(x,2)-pow(y,2));

    // G4cout << Bvec[0] << Bvec[1] << Bvec[2]<< G4endl;
} // EarthMagField::GetFieldValue
