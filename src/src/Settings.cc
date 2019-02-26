#include "Settings.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Settings {
// All these variables are put here to be shared amounts source file
// (not very c++ but easier to implement)
//

// Earth radius
    const G4double earthRadius = 6378.137 * km;

// parameters : initialization values, can be modified in the main code
    G4long RANDOM_SEED = 12345; // dummy value that will be replaced

    G4int NB_EVENT = 0;

    G4String CACHED_LENGTH = "10";   // for magnetic field solver, in meters

// Source parameters, geodetic coordinates ( = geographic = GPS)
    const G4double SOURCE_LAT = 11.01;     // degree
    const G4double SOURCE_LONG = -95.40;      // degree
    G4double SOURCE_ALT = 15.;       // km

    G4double OPENING_ANGLE = 40.;    // degree
    G4String BEAMING_TYPE = "Uniform";
    G4double TILT_ANGLE = 0.0;

    G4double TIME_LIMIT = 2.0 * second;
    G4double MIN_ENERGY_OUTPUT = 10.0 * keV;

    G4double SOURCE_SIGMA_TIME = 0.; // microsecond

// output altitudes
    std::vector<G4double> record_altitudes; // ! : geodetic altitudes (remark: when building the geometry, geocentric altitudes are used)

    G4bool MAG_FIELD_ON = true;

    G4bool USE_STEP_MAX_for_record = true; // force max step only for layers where particles are recorded
    G4bool USE_STEP_MAX_global = false; // force max step everywhere

    G4bool OUTPUT_ALT_LAYERS_TO_FILE = false; // output list of altitude and densities of layer to file (for debug)

    G4bool RECORD_ELEC_POSI_ONLY = true; // record only electron and positrons
    G4bool RECORD_PHOT_ONLY = false; // record only photons

    G4bool OUTPUT_ECEF_COORDS = false;
    G4bool OUTPUT_RadDist = false; // Radial distance can also be calculated a posteriori form ECEF data (e.g. Matlab routines)
}

////// Global static pointer used to ensure a single instance of the class.
///// (singleton pattern)
//Settings *Settings::m_pInstance = 0;

//void
//Settings::diplay_settings()
//{
//    G4cout << G4endl;
//    G4cout << "*************************************************************" << G4endl;
//    G4cout << " SIMULATION SETTINGS : " << G4endl;
//    G4cout << "    Random Number Seed : " << Rand_seed << G4endl;
//    G4cout << "    Cached length for Integrator : " << CACHED_LENGTH << " m" << G4endl;
//    G4cout << "    Initial altitude : " << SOURCE_ALT << " km" << G4endl;
//    G4cout << "    Beaming : " << G4endl;
//    G4cout << "        Type : " << BEAMING_TYPE << G4endl;
//    G4cout << "        Angle : " << OPENING_ANGLE << " degrees; corresponds to" << G4endl;
//    G4cout << "          - max angle if Type = uniform;" << G4endl;
//    G4cout << "          - sigma if Type = gaussian (normal)" << G4endl;

//    G4cout << "    Source timimg (gaussian sigma) : " << SOURCE_SIGMA_TIME << " microsecond" << G4endl;

//    G4cout << "    Output Altitudes : " << G4endl;
//    G4cout << "        ";

//    for (unsigned int ii = 0; ii < record_altitudes.size(); ++ii)
//        {
//            G4cout << record_altitudes[ii] << " km, ";
//        }

//    G4cout << G4endl;

//    G4cout << "*************************************************************" << G4endl;
//    G4cout << G4endl;
//}

//G4String Settings::BeamingType() const
//{
//    return BEAMING_TYPE;
//}

//void Settings::set_BeamingType(const G4String &value)
//{
//    BEAMING_TYPE = value;
//}

//void Settings::set_AltMax_recorded()
//{
//    G4double max_alt = 0.;

//    for (int ii = 0; ii < nb_altitude_record; ++ii)
//        {
//            if (record_altitudes[ii] > max_alt)
//                {
//                    max_alt = record_altitudes[ii];
//                }
//        }

//    ALT_MAX_RECORDED = max_alt;
//}

//G4bool Settings::USE_STEP_MAX() const
//{
//    return USE_STEP_MAX_;
//}

//void Settings::set_MAG_FIELD_ON(const G4bool MAG_FIELD_BOOL)
//{
//    MAG_FIELD_ON_ = MAG_FIELD_BOOL;
//}

//G4bool Settings::MAG_FIELD_ON() const
//{
//    return MAG_FIELD_ON_;
//}
