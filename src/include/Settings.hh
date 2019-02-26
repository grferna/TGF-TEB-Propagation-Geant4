//////////////////////////////////////////////////////////////////////////////////

//// /* GEANT4 code for propagation of gamma-rays, electron and positrons in Earth's atmosphere */
////
//// //
//// // ********************************************************************
//// // * License and Disclaimer                                           *
//// // *                                                                  *
//// // * The  Geant4 software  is  copyright of the Copyright Holders  of *
//// // * the Geant4 Collaboration.  It is provided  under  the terms  and *
//// // * conditions of the Geant4 Software License,  included in the file *
//// // * LICENSE and available at  http://cern.ch/geant4/license .  These *
//// // * include a list of copyright holders.                             *
//// // *                                                                  *
//// // * Neither the authors of this software system, nor their employing *
//// // * institutes,nor the agencies providing financial support for this *
//// // * work  make  any representation or  warranty, express or implied, *
//// // * regarding  this  software system or assume any liability for its *
//// // * use.  Please see the license in the file  LICENSE  and URL above *
//// // * for the full disclaimer and the limitation of liability.         *
//// // *                                                                  *
//// // * This  code  implementation is the result of  the  scientific and *
//// // * technical work of the GEANT4 collaboration.                      *
//// // * By using,  copying,  modifying or  distributing the software (or *
//// // * any work based  on the software)  you  agree  to acknowledge its *
//// // * use  in  resulting  scientific  publications,  and indicate your *
//// // * acceptance of all terms of the Geant4 Software license.          *
//// // ********************************************************************
//////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Settings {
// Earth radius
    extern const G4double earthRadius;

// parameters : initialization values, can be modified in the main code
    extern G4long RANDOM_SEED;

    extern G4int NB_EVENT;

    extern G4String CACHED_LENGTH;   // for magnetic field solver

// Source parameters, geodetic coordinates ( = geographic = GPS)
    extern const G4double SOURCE_LAT;     // degree
    extern const G4double SOURCE_LONG;      // degree
    extern G4double SOURCE_ALT;       // km

    extern G4double OPENING_ANGLE;    // degree
    extern G4String BEAMING_TYPE;
    extern G4double TILT_ANGLE;

    extern G4double TIME_LIMIT;
    extern G4double MIN_ENERGY_OUTPUT;

    extern G4double SOURCE_SIGMA_TIME; // microsecond

// output altitudes
    extern std::vector<G4double> record_altitudes; // ! : geodetic altitudes (remark: when building the geometry, geocentric altitudes are used)

    extern G4bool MAG_FIELD_ON;

    extern G4bool USE_STEP_MAX_for_record;
    extern G4bool USE_STEP_MAX_global;

    extern G4bool OUTPUT_ALT_LAYERS_TO_FILE;

    extern G4bool RECORD_ELEC_POSI_ONLY;

    extern G4bool RECORD_PHOT_ONLY;

    extern G4bool OUTPUT_ECEF_COORDS; // add ECEF (x,y,z) coordinates to output
    extern G4bool OUTPUT_RadDist; // add radial distance to output
}

//#pragma once

//#include "G4RunManager.hh"
//#include "G4UnitsTable.hh"
//#include "G4SystemOfUnits.hh"
//#include "G4PhysicalConstants.hh"

//#include <vector>

//// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#include <string>

//class Settings
//{
//    private:

//        Settings() {} // Private so that it can not be called

//        Settings(Settings const &) {}

//        // copy constructor is private
//        // assignment operator is private
//        static Settings *m_pInstance;

//        // Earth radius
//        const G4double earthRadius = 6378.137 * km;

//        // parameters : initialization values, can be modified in the main code
//        G4int Rand_seed = 78;

//        G4int NB_EVENT = 0;

//        G4String CACHED_LENGTH = "10";   // for magnetic field solver

//        // Source parameters, geodetic coordinates ( = geographic = GPS)
//        G4double SOURCE_LAT  = -13.;     // degree
//        G4double SOURCE_LONG = 32.;      // degree
//        G4double SOURCE_ALT  = 15;       // km
//        G4double ALT_MAX_RECORDED     = 500.;     // initial value

//        G4double OPENING_ANGLE = 45.;    // degree
//        G4String BEAMING_TYPE  = "Uniform";

//        G4double SOURCE_SIGMA_TIME = 0.; // microsecond

//        // output altitudes
//        G4int nb_altitude_record = 0;
//        std::vector < G4double > record_altitudes; // ! : geodetic altitudes (remark: when building the geometry, geocentric altitudes are used)

//        void set_AltMax_recorded();

//        G4bool MAG_FIELD_ON_ = false;

//        G4bool USE_STEP_MAX_ = false;

//    public:

//        static Settings *get_Instance()
//        {
//            if (m_pInstance == 0)   // Only allow one instance of class to be generated (lazy initialization)
//                {
//                    m_pInstance = new
//                    Settings();
//                }

//            return m_pInstance;
//        }

//        G4double AltMax_recorded() const
//        {
//            return ALT_MAX_RECORDED;
//        }

//        G4double AltMax_recorded_times_km() const
//        {
//            return ALT_MAX_RECORDED * km;
//        }

//        const G4String &CachedLength() const
//        {
//            return CACHED_LENGTH;
//        }

//        void set_CachedLength(const G4String &cachedLength)
//        {
//            CACHED_LENGTH = cachedLength;
//        }

//        G4int NbAltitudeRecord() const
//        {
//            return nb_altitude_record;
//        }

//        G4int NbEvent() const
//        {
//            return NB_EVENT;
//        }

//        void set_NbEvent(G4int nbEvent)
//        {
//            NB_EVENT = nbEvent;
//        }

//        G4double OpeningAngle() const
//        {
//            return OPENING_ANGLE;
//        }

//        void set_OpeningAngle(G4double openingAngle)
//        {
//            OPENING_ANGLE = openingAngle;
//        }

//        const std::vector < G4double > &Output_Altitudes() const
//        {
//            return record_altitudes;
//        }

//        void set_OutputAltitudes(const std::vector < G4double > &outputAltitudes)
//        {
//            record_altitudes = outputAltitudes;
//        }

//        G4int RandSeed() const
//        {
//            return Rand_seed;
//        }

//        void set_RandSeed(G4int randSeed)
//        {
//            Rand_seed = randSeed;
//        }

//        G4double SourceAlt() const
//        {
//            return SOURCE_ALT;
//        }

//        void set_SourceAlt(G4double sourceAlt)
//        {
//            SOURCE_ALT = sourceAlt;
//        }

//        G4double SourceLat() const
//        {
//            return SOURCE_LAT;
//        }

//        void set_SourceLat(G4double sourceLat)
//        {
//            SOURCE_LAT = sourceLat;
//        }

//        G4double SourceLong() const
//        {
//            return SOURCE_LONG;
//        }

//        void set_SourceLong(G4double sourceLong)
//        {
//            SOURCE_LONG = sourceLong;
//        }

//        G4double SourceSigmaTime() const
//        {
//            return SOURCE_SIGMA_TIME;
//        }

//        void set_SourceSigmaTime(G4double sourceSigmaTime)
//        {
//            SOURCE_SIGMA_TIME = sourceSigmaTime;
//        }

//        void add_OutputAltitude(G4double altitude)
//        {
//            record_altitudes.push_back(altitude);
//            nb_altitude_record++;
//            set_AltMax_recorded();
//        }

//        std::vector < G4double > Output_Altitudes()
//        {
//            return record_altitudes;
//        }

//        G4double EarthRadius() const
//        {
//            return earthRadius;
//        }

//        void diplay_settings();

//        G4String BeamingType() const;
//        void     set_BeamingType(const G4String &value);
//        G4bool MAG_FIELD_ON() const;
//        void set_MAG_FIELD_ON(const G4bool MAG_FIELD_BOOL);
//        G4bool USE_STEP_MAX() const;
//};
