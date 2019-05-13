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

#include "PrimaryGeneratorAction.hh"
#include "EarthMagField.hh"

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField::EarthMagField() {

    /* Memory allocation */

    strncpy(VersionDate, VERSIONDATE_LARGE + 39, 11);
    VersionDate[11] = '\0';

    char filename[] = "WMM.COF";

    if (!MAG_robustReadMagModels(filename, &MagneticModels, epochs)) {
        printf("\n WMM.COF not found.  Aborting... \n ");
        std::abort();
    }

    if (nMax < MagneticModels[0]->nMax) nMax = MagneticModels[0]->nMax;
    NumTerms = ((nMax + 1) * (nMax + 2) / 2);
    TimedMagneticModel = MAG_AllocateModelMemory(NumTerms); /* For storing the time modified WMM Model parameters */

    MAG_SetDefaults(&Ellip, &Geoid); /* Set default values and constants */

    Geoid.Geoid_Initialized = 1;
    /* Set EGM96 Geoid parameters END */

    UserDate.Year = settings->year;
    UserDate.Day = settings->day;
    UserDate.Month = settings->month;
    UserDate.DecimalYear = double(UserDate.Year) + UserDate.Month / 12.0 + UserDate.Day / 365.24;

    MAG_TimelyModifyMagneticModel(UserDate, MagneticModels[0],
                                  TimedMagneticModel); /* Time adjust the coefficients, Equation 19, WMM Technical report */

//    typedef struct {
//        double lambda; /* longitude */
//        double phi; /* geodetic latitude */
//        double HeightAboveEllipsoid; /* height above the ellipsoid (HaE) */
//        double HeightAboveGeoid; /* (height above the EGM96 geoid model ) */
//        int UseGeoid;
//    } MAGtype_CoordGeodetic;


}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EarthMagField::~EarthMagField() {}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EarthMagField::GetFieldValue(const double Point[3], double *Bfield) const {
    //  geodetic_converter::GeodeticConverter g_geodetic_converter;

    xx = Point[0] / m;
    yy = Point[1] / m;
    zz = Point[2] / m;

    // conversion ECEF to geodetic

    // input x y z in meters, output altitude in km, output lat lon in degrees
    geod_conv::GeodeticConverter::ecef2Geodetic(xx, yy, zz, lat, lon, alt);

    if (alt < 45000.0) {
        Bfield[0] = 0;
        Bfield[1] = 0;
        Bfield[2] = 0;
        return;
    }

    elong = lon; // degrees
    colat = (180.0 / 2.0 - lat); // degrees
    alt_km = alt / 1000.0;

    //
    // // Get altitude from latitude
    // G4double altt = p*coslat + (z + e2*Nphi*sinlat)*sinlat - Nphi;

    //     itype = 1 if geodetic (spheroid)
    //     itype = 2 if geocentric (sphere)
    //     alt   = height in km above sea level if itype = 1
    //           = distance from centre of Earth in km if itype = 2 (>3485 km)                                              // degrees

    // IGRF magnetic field
//    igrf12syn_(&isv, &date, &itype, &alt_km, &colat, &elong, &Bx, &By, &Bz, &f, &ier); // gives NED (North East Down) components, in nT

    CoordGeodetic.lambda = lon;
    CoordGeodetic.phi = lat;
    CoordGeodetic.HeightAboveEllipsoid = alt_km;

//    if (MAG_GetUserInput(MagneticModels[0], &Geoid, &CoordGeodetic, &UserDate) == 1) /*Get User Input */
//    {

    MAG_GeodeticToSpherical(Ellip, CoordGeodetic,
                            &CoordSpherical); /*Convert from geodetic to Spherical Equations: 17-18, WMM Technical report*/

    MAG_Geomag(Ellip, CoordSpherical, CoordGeodetic, TimedMagneticModel,
               &GeoMagneticElements); /* Computes the geoMagnetic field elements and their time change*/
//    }

    Bx = GeoMagneticElements.X; // X Y Z output in nanoTesla
    By = GeoMagneticElements.Y;
    Bz = GeoMagneticElements.Z;

    //     G4cout << Bx << " " << By << " " << Bz << G4endl;

    // NED (North East Down) to ECEF convertion

    // NOT USING EIGEN TO GET MORE PERFORMANCE

    // Local north direction      //Local east direction     // Local vertical (down)
    //        MAT_enu_ecef << -coslon*sinlat,            -sinlon,                  -coslon*coslat,
    //                        -sinlon*sinlat,             coslon,                  -sinlon*coslat,
    //                  coslat,                      0.   ,                     -sinlat;

    //        Bvec=MAT_enu_ecef*Bvec;

    lat = lat * degree;   // degree to radians
    lon = lon * degree;   // degree to radians

#ifdef _WIN32
	sinlat = sin(lat);
	coslat = cos(lat);
	sinlon = sin(lon);
	coslon = cos(lon);
#else
	sincos(lon, &sinlon, &coslon); // important for performance
	sincos(lat, &sinlat, &coslat);
#endif

    Bfield_ecef_x = -coslon * sinlat * Bx - sinlon * By - coslon * coslat * Bz;
    Bfield_ecef_y = -sinlon * sinlat * Bx + coslon * By - sinlon * coslat * Bz;
    Bfield_ecef_z = coslat * Bx - sinlat * Bz;

    Bfield[0] = Bfield_ecef_x * nano_tesla_to_G4;
    Bfield[1] = Bfield_ecef_y * nano_tesla_to_G4;
    Bfield[2] = Bfield_ecef_z * nano_tesla_to_G4;

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
