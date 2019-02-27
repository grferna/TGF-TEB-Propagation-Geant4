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
#include <fortran.hh>

using namespace std;

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//C     INPUT VARIABLES:
//C        IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
//C              (Year ignored in current model)
//C        SEC - UT(SEC)
//C        ALT - ALTITUDE(KM)
//C        GLAT - GEODETIC LATITUDE(DEG)
//C        GLONG - GEODETIC LONGITUDE(DEG)
//C        STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
//C        F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
//C        F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
//C        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
//C           - ARRAY CONTAINING:
//C             (1) DAILY AP
//C             (2) 3 HR AP INDEX FOR CURRENT TIME
//C             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
//C             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
//C             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
//C             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
//C                    TO CURRENT TIME
//C             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
//C                    TO CURRENT TIME
//C        MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
//C                 CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
//C                 MASS 17 IS Anomalous O ONLY.)
//C
//C     NOTES ON INPUT VARIABLES:
//C        UT, Local Time, and Longitude are used independently in the
//C        model and are not of equal importance for every situation.
//C        For the most physically realistic calculation these three
//C        variables should be consistent (STL=SEC/3600+GLONG/15).
//C        The Equation of Time departures from the above formula
//C        for apparent local time can be included if available but
//C        are of minor importance.
//c
//C        F107 and F107A values used to generate the model correspond
//C        to the 10.7 cm radio flux at the actual distance of the Earth
//C        from the Sun rather than the radio flux at 1 AU. The following
//C        site provides both classes of values:
//C        ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
//C
//C        F107, F107A, and AP effects are neither large nor well
//C        established below 80 km and these parameters should be set to
//C        150., 150., and 4. respectively.
//C
//C     OUTPUT VARIABLES:
//C        D(1) - HE NUMBER DENSITY(CM-3)
//C        D(2) - O NUMBER DENSITY(CM-3)
//C        D(3) - N2 NUMBER DENSITY(CM-3)
//C        D(4) - O2 NUMBER DENSITY(CM-3)
//C        D(5) - AR NUMBER DENSITY(CM-3)
//C        D(6) - TOTAL MASS DENSITY(GM/CM3)
//C        D(7) - H NUMBER DENSITY(CM-3)
//C        D(8) - N NUMBER DENSITY(CM-3)
//C        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
//C        T(1) - EXOSPHERIC TEMPERATURE
//C        T(2) - TEMPERATURE AT ALT

// IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,MASS,D,T

// extrernal fortran subroutine to get MSIS atmospheric densities
extern "C" {
void gtd7_(INTEGER &IYD, // YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
           REAL &SEC, // UT(SEC)
           REAL &ALT, // ALTITUDE(KM)
           REAL &GLAT, // GEODETIC LATITUDE(DEG)
           REAL &GLONG, // GEODETIC LONGITUDE(DEG)
           REAL &STL, // LOCAL APPARENT SOLAR TIME
           REAL &F107A, // 81 day AVERAGE OF F10.7 FLUX (centered on day DDD
           REAL &F107, // DAILY F10.7 FLUX FOR PREVIOUS DAY
           REAL &AP,  // MAGNETIC INDEX(DAILY)
           INTEGER &MASS, // MASS NUMBER
           REAL *D, REAL *T); // OUTPUT VARIABLES temperatures
}

TGFDetectorConstruction::TGFDetectorConstruction()
{


    logicalWorld = nullptr;
    physicalWorld = nullptr;

    globalfieldMgr = nullptr;

    if (settings->OUTPUT_ALT_LAYERS_TO_FILE)
    {
        asciiFile.open("alt_dens_debug.txt", std::ios::trunc);

        asciiFile << "altitude (km) // density (g/cm2)" << G4endl;
    }

}

TGFDetectorConstruction::~TGFDetectorConstruction() = default;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

G4VPhysicalVolume *TGFDetectorConstruction::Construct()
{
    //    G4FieldManager *null_field = nullptr;

    // cleaning geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::Clean();
    G4LogicalVolumeStore::Clean();
    G4SolidStore::Clean();

    // construct field managers (for magnetic fields)
    if (settings->MAG_FIELD_ON)
    {
        Construct_MagField_Managers();
    }

    // generating the world layers of constant density
    calculate_altitudes_list();

    //    for (int jj = 0; jj < altitudes.size(); ++jj) {
    //
    //        G4cout << altitudes[jj] << G4endl;
    //    }

    // World Material : Vaccum
    G4NistManager *man = G4NistManager::Instance();
    vac = man->FindOrBuildMaterial("G4_Galactic");

    // World solid
    G4Sphere *solidWorld;
    solidWorld = new G4Sphere("world_S", settings->earthRadius, (settings->earthRadius + world_max_altitude), 0 * degree, 360 * degree, 0 * degree, 180 * degree); // earth radius + 100 km
    // World logical

    logicalWorld = new G4LogicalVolume(solidWorld,                      // solid
                                       vac,                                                     // material
                                       "world_L");

    if (settings->MAG_FIELD_ON)
    {
        logicalWorld->SetFieldManager(globalfieldMgr, true);
    }

    // Physical volume
    physicalWorld = new G4PVPlacement(nullptr, G4ThreeVector(), "world_P", // name (2nd constructor)
                                      logicalWorld,             // logical volume
                                      nullptr,                      // mother volume
                                      false,              // no boolean operation
                                      0);                          // copy number

    G4VisAttributes *VisAttWorld = new G4VisAttributes(G4Colour(204 / 255., 255 / 255., 255 / 255.));
    logicalWorld->SetVisAttributes(VisAttWorld);

    // setting default (world) region info
    G4Region *defaultRegion = (*(G4RegionStore::GetInstance()))[0]; // the default (world) region is index 0 in the region store
    auto *defaultRInfo = new RegionInformation();
    defaultRInfo->Set_World();
    defaultRegion->SetUserInformation(defaultRInfo);

    // volume de détection (+comblage du monde)
    //         G4Sphere* detec_sol400=new G4Sphere("detec_sol400",
    //                                             (shared_var::earthRadius+shared_var::ALT_DETECT*km),
    //                                             shared_var::earthRadius+world_max_altitude,
    //                                             0*degree,
    //                                             360*degree,
    //                                             0*degree,
    //                                             180*degree );
    //         G4LogicalVolume* detec_log400=new G4LogicalVolume(detec_sol400, vac, "detec_log400", 0, 0, 0);
    //
    //         detec_log400->SetFieldManager(globalfieldMgr, true);
    //
    //         G4VPhysicalVolume* detec_phys400=new G4PVPlacement( 0,G4ThreeVector(),"detec_phys400",detec_log400, physicalWorld,false,0,0);

    // Make Invisible
    //  logicalWorld -> SetVisAttributes(G4VisAttributes::Invisible);

    std::vector<G4Material *> Airs = TGFDetectorConstruction::Construct_Atmos_layers_Materials(altitudes_geodetic); // used MSIS C++ code; creates the Airs[jj]

    // considered atmosphere region (particle flying ou thsi region are killed)

    auto *Reginfo = new RegionInformation();
    Reginfo->set_considered_atmosphere();
    considered_atmos_Region->SetUserInformation(Reginfo);

    // XrayTelDetectorConstruction::ConstructAtmosMats2(); // Use datafile generated from MSIS website
    // !! : significant differences (>10 %) between the two

    G4double innerRad = 0;
    G4double outerRad = 0;

    G4int id_SD = 0;

    // atmosphere construction

    const G4double ALT_MAX_RECORDED = *std::max_element(settings->record_altitudes.begin(), settings->record_altitudes.end());

    for (unsigned int jj = 0; jj < altitudes_geodetic.size() - 1; jj++) // geocentric altitudes
    {
        innerRad = settings->earthRadius + altitudes_geodetic[jj];
        outerRad = settings->earthRadius + altitudes_geodetic[jj + 1];

        atmosLayers_S.push_back(new G4Sphere("atmosLayer_S_" + std::to_string(jj), innerRad, outerRad, 0 * degree, 360 * degree, 0 * degree, 180 * degree));

        if ((innerRad + outerRad) / 2. < (settings->earthRadius + ALT_MAX_RECORDED * km))
        {
            atmosLayers_LV.push_back(new G4LogicalVolume(atmosLayers_S.back(), Airs[jj], "atmosphere_LV_" + std::to_string(jj), nullptr, nullptr, nullptr)); // put vaccum to test if angle sampling is OK
        }
        else                                     // vaccum if altitude > ALT_MAX
        {
            atmosLayers_LV.push_back(new G4LogicalVolume(atmosLayers_S.back(), vac, "atmosphere_LV_" + std::to_string(jj), nullptr, nullptr, nullptr));
        }

        // assigning sensitive detector
        for (G4double rec_alt : settings->record_altitudes)
        {
            if (altitudes_geodetic[jj] == rec_alt * km)
            {
                sens_det_List.push_back(new SensitiveDet("sens_det_" + std::to_string(id_SD), id_SD, altitudes_geodetic[jj] / km));
                atmosLayers_LV.back()->SetSensitiveDetector(sens_det_List.back());
                G4SDManager::GetSDMpointer()->AddNewDetector(sens_det_List.back());
                id_SD++;

                if (settings->USE_STEP_MAX_for_record)
                {
                    atmosLayers_LV.back()->SetUserLimits(stepLimit);
                }
            }
        }

        // assigning magnetic field
        if (settings->MAG_FIELD_ON)
        {
            atmosLayers_LV.back()->SetFieldManager(globalfieldMgr, true);

            // null magnetic field if altitude < 45 km
            if ((innerRad + outerRad) / 2. < (settings->earthRadius + alt_account_Mag_Field))
            {
                atmosLayers_LV.back()->SetFieldManager(Null_FieldManager, true);
            }
        }

        // setting to the 'considered atmosphere' region
        atmosLayers_LV.back()->SetRegion(considered_atmos_Region);
        considered_atmos_Region->AddRootLogicalVolume(atmosLayers_LV.back());

        G4String name_PV = "atmosphere_PV_" + std::to_string(jj);
        atmosLayers_PV.push_back(new G4PVPlacement(nullptr, G4ThreeVector(), name_PV, atmosLayers_LV.back(), physicalWorld, false, 0, false));
    }

    G4cout << G4endl << "Geometry built successfully." << G4endl << G4endl;

    return physicalWorld;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void TGFDetectorConstruction::calculate_altitudes_list()
// fills the vector altitudes
{
    const G4double ALT_MAX_RECORDED = *std::max_element(settings->record_altitudes.begin(), settings->record_altitudes.end());

    const G4double alt_max_construction = min(alt_max_atmosphere, ALT_MAX_RECORDED * km);
    // it is either 200 km, either the maximum detection altitude if smaller than 200 km

    // defining the altitude vector
    for (G4int jj = 0; jj < nb_altitudes; jj++)   // geocentric altitudes
    {
        altitudes_geodetic.push_back(exp(log(alt_min) + (log(alt_max_construction) - log(alt_min)) * double(jj) / double(nb_altitudes - 1)));
    }

    // adding an extra altitude 1km more than the max of the recorded altitude
    altitudes_geodetic.push_back(ALT_MAX_RECORDED * km + 1. * km);

    // adding the record volume : a thin layer volume starting at the record altitude

    for (G4double rec_alt : settings->record_altitudes)
    {
        if (not_contains(rec_alt * km, altitudes_geodetic))
        {
            altitudes_geodetic.push_back(rec_alt * km);
        }

        G4double rec_alt2 = rec_alt * km + 0.01 * km;

        if (not_contains(rec_alt2, altitudes_geodetic))
        {
            altitudes_geodetic.push_back(rec_alt2);
        }
    }

    // sorting in increasing value
    std::sort(altitudes_geodetic.begin(), altitudes_geodetic.end());

    if (hasDuplicates(altitudes_geodetic))
    {
        G4cout << "ERROR : There are duplicates values in the altitude list. Aborting." << G4endl;
        std::abort();
    }
}

bool TGFDetectorConstruction::hasDuplicates(const std::vector<G4double> &arr)
{
    for (uint i = 0; i < arr.size(); ++i)
    {
        for (uint j = i + 1; j < arr.size(); ++j)
        {
            if (arr[i] == arr[j])
            {
                return true;
            }
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
G4bool TGFDetectorConstruction::not_contains(G4double x, const std::vector<G4double> &v)
{
    return !(std::find(v.begin(), v.end(), x) != v.end());
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// Caltulating materials of the atmopsheric layers, based on the MSIS C++ model integrated to this code
// ref : https://ccmc.gsfc.nasa.gov/modelweb/atmos/nrlmsise00.html

std::vector<G4Material *> TGFDetectorConstruction::Construct_Atmos_layers_Materials(const std::vector<G4double> altitudes_)
{
    std::vector<G4Material *> Airs;

    // Vaccum
    G4NistManager *man = G4NistManager::Instance();
    G4Material *vaccum = man->FindOrBuildMaterial("G4_Galactic");

    //    elHe = new G4Element(name = "Helium", symbol = "He", z = 2., He_molarMass);
    //    elH  = new G4Element(name = "Hydrogen", symbol = "H", z = 1., H_molarMass);

    for (uint idx_alt = 0; idx_alt < altitudes_.size() - 1; idx_alt++)
    {
        const double innerAlt = altitudes_[idx_alt];
        const double outerAlt = altitudes_[idx_alt + 1];
        const double altitude_in_km = (innerAlt + outerAlt) / 2. / km; // geocentric altitude

        //            G4cout << altitude_in_km << G4endl;

        if (altitude_in_km > alt_max_atmosphere / km)
        {
            Airs.push_back(vaccum);
        }
        else
        {
            INTEGER input_iyd = 172; // IYD - YEAR AND DAY AS YYDDD
            REAL input_sec = 29000.0;
            REAL input_alt = (REAL) altitude_in_km;
            REAL input_g_lat = (REAL) settings->SOURCE_LAT;
            REAL input_g_long = (REAL) settings->SOURCE_LONG;
            REAL input_lst = 16.0;
            REAL input_f107A = 150.0;
            REAL input_f107 = 150.0;
            REAL input_ap = 4.0;
            INTEGER input_mass = 48;
            REAL output_D[9];
            REAL output_T[2];

            // G4cout << altitude << G4endl;

            gtd7_(input_iyd, input_sec, input_alt, input_g_lat, input_g_long, input_lst, input_f107A, input_f107, input_ap, input_mass, output_D, output_T); // MSIS, fortran function call

            if (std::isnan(output_D[5]) || std::isinf(isnan(output_D[5])))
            {
                G4cout << "ERROR : density from gtd7_ is NaN. Aborting" << G4endl;
                std::abort();
            }

            G4double density_air = output_D[5] * g / cm3; // getting density and converting it to the GEANT4 system of unit

            if (settings->OUTPUT_ALT_LAYERS_TO_FILE)
            {
                asciiFile << altitude_in_km << " " << output_D[5] << G4endl;
            }

            Airs.push_back(man->BuildMaterialWithNewDensity("Air_" + std::to_string(idx_alt), "G4_AIR", density_air));
        }
    }

    G4cout << "Successfully created air materials list" << G4endl;

    return Airs;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

void TGFDetectorConstruction::Construct_MagField_Managers()
{
    /////////// Magnetic field

    myEarthMagField = new EarthMagField_alt;
    //    myEarthMagField = new EarthMagField;

    G4double distanceConst = std::stod(settings->CACHED_LENGTH) * meter;

    myCachedEarthMagField = new G4CachedMagneticField(myEarthMagField, distanceConst);

    globalfieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

    if (settings->CACHED_LENGTH == "0")
    {
        //            globalfieldMgr->CreateChordFinder(myEarthMagField);
        pMagFldEquation = new G4Mag_UsualEqRhs(myEarthMagField);
        fStepper = new G4DormandPrince745(pMagFldEquation);
        fChordFinder = new G4ChordFinder(myEarthMagField, fMinStep, fStepper);
        globalfieldMgr->SetChordFinder(fChordFinder);
        globalfieldMgr->SetDetectorField(myEarthMagField);
    }
    else
    {
        //            globalfieldMgr->CreateChordFinder(myCachedEarthMagField);
        pMagFldEquation = new G4Mag_UsualEqRhs(myCachedEarthMagField);
        fStepper = new G4DormandPrince745(pMagFldEquation);
        fChordFinder = new G4ChordFinder(myCachedEarthMagField, fMinStep, fStepper);
        globalfieldMgr->SetChordFinder(fChordFinder);
        globalfieldMgr->SetDetectorField(myCachedEarthMagField);
    }

    //    G4cout << "CACHED_LENGTH = " << settings->CachedLength() << G4endl;
    //    !!! RQ : avoid G4NystromRK4 , bug with use of G4CachedMagneticField
    //           avoid G4ConstRK4 and G4ImplicitEuler : bug with G4CachedMagneticField
    //               avoid G4DormandPrinceRK78
    //    G4DormandPrince745 is the best, from first tests
    //    CACHED_LENGTH at 1000 m gives bad results
    //   this was tested with G4 10.02 may be fixed with further updates.
    //
    globalfieldMgr->SetMinimumEpsilonStep(minEps);
    globalfieldMgr->SetMaximumEpsilonStep(maxEps);
    globalfieldMgr->SetDeltaOneStep(1000. * cm / 2);
    globalfieldMgr->SetDeltaIntersection(1000. * cm / 2);
    globalfieldMgr->GetChordFinder()->SetDeltaChord(1000. * cm / 2);

    //// NUll field manager
    ///
    Null_FieldManager = new G4FieldManager();
    magField_null = new G4UniformMagField(G4ThreeVector(0., 0., 0.));
    Null_FieldManager->SetDetectorField(magField_null);
    Null_FieldManager->CreateChordFinder(magField_null);
    Null_FieldManager->SetMinimumEpsilonStep(minEps);
    Null_FieldManager->SetMaximumEpsilonStep(maxEps);
    Null_FieldManager->SetDeltaOneStep(1000. * cm / 2);
    Null_FieldManager->SetDeltaIntersection(1000. * cm / 2);
    Null_FieldManager->GetChordFinder()->SetDeltaChord(1000. * cm / 2);

    // set maximum acceptable step everywhere, not used anymore
//    G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(settings->STEP_MAX_VAL);

}

//// Based on MSIS model from NASA website
// void XrayTelDetectorConstruction::ConstructAtmosMats2() {
//
//    G4String name, symbol;
//    G4double z;
//    G4int ncomponents, ncomponents2, natoms;
//    G4double density, densityO, densityN2, densityO2, densityHe, densityAr, densityH, densityN, density_total;
//    G4double NumberDensityO, NumberDensityN2, NumberDensityO2, NumberDensityHe, NumberDensityAr;
//    G4double NumberDensityH, NumberDensityN;
//    G4double proporN2, proporO2, proporAr, proporO, proporHe, proporH, proporN;
//
//    G4double NmolarMass = 14.0067 * g / mole;
//    G4double OmolarMass = 15.9994 * g / mole;
//    G4double ArmolarMass = 39.948 * g / mole;
//    G4double HemolarMass = 4.002602 * g / mole;
//    G4double HmolarMass = 1.00794 * g / mole;
//
//    G4double Avog = 6.02214179e23;
//
//    G4double NmolarMass2 = 14.0067;
//    G4double OmolarMass2 = 15.9994;
//    G4double ArmolarMass2 = 39.948;
//    G4double HemolarMass2 = 4.002602;
//    G4double HmolarMass2 = 1.00794;
//    G4double N2molarMass2 = 28.01340;
//    G4double O2molarMass2 = 31.99880;
//
//    elN = new G4Element(name = "Nitrogen", symbol = "N", z = 7., NmolarMass);
//    elO = new G4Element(name = "Oxygen", symbol = "O", z = 8., OmolarMass);
//    elA = new G4Element(name = "Argon", symbol = "Ar", z = 18., ArmolarMass);
//    elHe = new G4Element(name = "Helium", symbol = "He", z = 2., HemolarMass);
//    elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1., HmolarMass);
//
//// read data file
//
//    ReadInputAtmosFile();
//
//// construct atmosphere layer materials
//    for (unsigned int i = 0; i < altitudes.size() - 1; i++) {
//
//        const double innerAlt = altitudes[i];
//        const double outerAlt = altitudes[i + 1];
//        double altitude = (innerAlt + outerAlt) / 2. / km;
//
//        density = interp1(FileAltitudes, FileDensities, altitude) * g / cm3;
//
//        NumberDensityO = interp1(FileAltitudes, FileO, altitude);
//        NumberDensityN2 = interp1(FileAltitudes, FileN2, altitude);
//        NumberDensityO2 = interp1(FileAltitudes, FileO2, altitude);
//        NumberDensityHe = interp1(FileAltitudes, FileHe, altitude);
//        NumberDensityAr = interp1(FileAltitudes, FileAr, altitude);
//        NumberDensityH = interp1(FileAltitudes, FileH, altitude);
//        NumberDensityN = interp1(FileAltitudes, FileN, altitude);
//
//// G4cout << altitude << " " << density /(g/cm3) << G4endl;
//// look how many elements in Air with non-null densities and compute total and proportions
//
//        densityN2 = NumberDensityN2 * N2molarMass2 / Avog * g / cm3;
//        densityO2 = NumberDensityO2 * O2molarMass2 / Avog * g / cm3;
//        densityAr = NumberDensityAr * ArmolarMass2 / Avog * g / cm3;
//        densityO = NumberDensityO * OmolarMass2 / Avog * g / cm3;
//        densityHe = NumberDensityHe * HemolarMass2 / Avog * g / cm3;
//        densityH = NumberDensityH * HmolarMass2 / Avog * g / cm3;
//        densityN = NumberDensityN * NmolarMass2 / Avog * g / cm3;
//
//        ncomponents2 = 0;
//        density_total = 0.;
//
//        if (NumberDensityHe > 1.) {
//            density_total = density_total + densityHe, ncomponents2++;
//        }
//        if (NumberDensityO > 1.) {
//            density_total = density_total + densityO, ncomponents2++;
//        }
//        if (NumberDensityN2 > 1.) {
//            density_total = density_total + densityN2, ncomponents2++;
//        }
//        if (NumberDensityO2 > 1.) {
//            density_total = density_total + densityO2, ncomponents2++;
//        }
//        if (NumberDensityAr > 1.) {
//            density_total = density_total + densityAr, ncomponents2++;
//        }
//        if (NumberDensityH > 1.) {
//            density_total = density_total + densityH, ncomponents2++;
//        }
//        if (NumberDensityN > 1.) {
//            density_total = density_total + densityN, ncomponents2++;
//        }
//
//        proporN2 = densityN2 / density_total;
//        proporO2 = densityO2 / density_total;
//        proporAr = densityAr / density_total;
//        proporO = densityO / density_total;
//        proporHe = densityHe / density_total;
//        proporH = densityH / density_total;
//        proporN = densityN / density_total;
//
//        //G4cout << proporN2 << " " << proporO2 << " " << proporAr << G4endl;
//
//// to remvoe warning messages
//        if (densityN2 <= 1.e-24) densityN2 = 1.;
//        if (densityO2 <= 1.e-24) densityO2 = 1.;
//        if (densityAr <= 1.e-24) densityAr = 1.;
//        if (densityO <= 1.e-24) densityO = 1.;
//        if (densityHe <= 1.e-24) densityHe = 1.;
//        if (densityH <= 1.e-24) densityH = 1.;
//        if (densityN <= 1.e-24) densityN = 1.;
//
//        N2 = new G4Material(name = "N2_" + std::to_string(i), densityN2, ncomponents = 1);
//        N2->AddElement(elN, natoms = 2);
//
//        O2 = new G4Material(name = "O2_" + std::to_string(i), densityO2, ncomponents = 1);
//        O2->AddElement(elO, natoms = 2);
//
//        Ar = new G4Material(name = "Ar_" + std::to_string(i), densityAr, ncomponents = 1);
//        Ar->AddElement(elA, natoms = 1);
//
//        O = new G4Material(name = "O_" + std::to_string(i), densityO, ncomponents = 1);
//        O->AddElement(elO, natoms = 1);
//
//        He = new G4Material(name = "He_" + std::to_string(i), densityHe, ncomponents = 1);
//        He->AddElement(elHe, natoms = 1);
//
//        H = new G4Material(name = "H_" + std::to_string(i), densityH, ncomponents = 1);
//        H->AddElement(elH, natoms = 1);
//
//        N = new G4Material(name = "N_" + std::to_string(i), densityN, ncomponents = 1);
//        N->AddElement(elH, natoms = 1);
//
//        Airs.push_back(new G4Material(name = "air_" + std::to_string(i), density, ncomponents2));
//
//        if (NumberDensityHe > 1.) {
//            Airs[i]->AddMaterial(He, proporHe);
//        }
//        if (NumberDensityO > 1.) {
//            Airs[i]->AddMaterial(O, proporO);
//        }
//        if (NumberDensityN2 > 1.) {
//            Airs[i]->AddMaterial(N2, proporN2);
//        }
//        if (NumberDensityO2 > 1.) {
//            Airs[i]->AddMaterial(O2, proporO2);
//        }
//        if (NumberDensityAr > 1.) {
//            Airs[i]->AddMaterial(Ar, proporAr);
//        }
//        if (NumberDensityH > 1.) {
//            Airs[i]->AddMaterial(H, proporH);
//        }
//        if (NumberDensityN > 1.) {
//            Airs[i]->AddMaterial(N, proporN);
//        }
//    }
// }
//
//// based NIST air and variations of density according to MSIS NASA website values
// void XrayTelDetectorConstruction::ConstructAtmosMats3() {
//
//    G4double density;
//
//// read data file
//
//    ReadInputAtmosFile();
//
//    G4NistManager* man = G4NistManager::Instance();
//
//// construct atmosphere layer materials
//    for (unsigned int i = 0; i < altitudes.size() - 1; i++) {
//
//        const double innerAlt = altitudes[i];
//        const double outerAlt = altitudes[i + 1];
//        double altitude = (innerAlt + outerAlt) / 2. / km;
//
//        density = interp1(FileAltitudes, FileDensities, altitude) * g / cm3;
//
//        Airs.push_back(man->BuildMaterialWithNewDensity("air_" + std::to_string(i), "G4_AIR", density));
//
//    }
//
// }
//
//// reading the data file generated from website NASA MSIS E 90 model
// void XrayTelDetectorConstruction::ReadInputAtmosFile() {
//    static double sss, d, e, f, gg, h, i, j, k;
//
//    ifstream file("../atmosphere/misise90.dat");
//
//    if (file) {
//        //L'ouverture s'est bien passée, on peut donc lire
//
//        string ligne; //Une variable pour stocker les lignes lues
//
//        while (!file.eof()) //Tant qu'on n'est pas à la fin, on lit
//        {
//
//            file >> sss >> d >> e >> f >> gg >> h >> i >> j >> k;
//
//            FileAltitudes.push_back(sss);
//            FileDensities.push_back(gg);
//            FileN2.push_back(e);
//            FileO2.push_back(f);
//            FileAr.push_back(i);
//            FileO.push_back(d);
//            FileN.push_back(k);
//            FileH.push_back(j);
//            FileHe.push_back(h);
//
//        }
//    } else {
//        cout << "ERROR: Impossible to open file for reading." << endl;
//        std::abort();
//    }
// }

// linear interporlation

//double TGFDetectorConstruction::interp1(vector<double>x, vector<double>y, double x_new)
//{
//    double y_new;

//    std::vector<double> dx, dy, slope, intercept;
//    dx.reserve(x.size());
//    dy.reserve(x.size());
//    slope.reserve(x.size());
//    intercept.reserve(x.size());

//    for (unsigned int i = 0; i < x.size(); ++i)
//        {
//            if (i < x.size() - 1)
//                {
//                    dx.push_back(x[i + 1] - x[i]);
//                    dy.push_back(y[i + 1] - y[i]);
//                    slope.push_back(dy[i] / dx[i]);
//                    intercept.push_back(y[i] - x[i] * slope[i]);
//                }
//            else
//                {
//                    dx.push_back(dx[i - 1]);
//                    dy.push_back(dy[i - 1]);
//                    slope.push_back(slope[i - 1]);
//                    intercept.push_back(intercept[i - 1]);
//                }
//        }

//    int idx = findNearestNeighbourIndex(x_new, x);
//    y_new = slope[idx] * x_new + intercept[idx];

//    return y_new;
//}

//int TGFDetectorConstruction::findNearestNeighbourIndex(double value, vector<double>x)
//{
//    double dist = FLT_MAX;
//    int    idx  = -1;

//    for (unsigned int i = 0; i < x.size(); ++i)
//        {
//            double newDist = value - x[i];

//            if ((newDist > 0) && (newDist < dist))
//                {
//                    dist = newDist;
//                    idx  = i;
//                }
//        }

//    return idx;
//}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
