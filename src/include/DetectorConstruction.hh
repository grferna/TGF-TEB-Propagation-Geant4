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

#include "EarthMagField_WMM.hh"
#include "EarthMagField_IGRF.hh"

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Sphere.hh"
#include <vector>

#include "G4VPVParameterisation.hh"
#include "G4Element.hh"

extern "C" {
#include <nrlmsise-00.h>
}

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <iostream>
#include <fstream>
#include <string>
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "SD.hh"

#include "RegionInformation.hh"
#include "G4RegionStore.hh"

#include <PrimaryGeneratorAction.hh>
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4String.hh"
#include "G4KM_DummyField.hh"
#include "G4UniformMagField.hh"
#include "G4CachedMagneticField.hh"
#include "G4PVParameterised.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "EarthMagField_WMM.hh"
#include <string>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4IntegrationDriver.hh"
#include "G4VIntegrationDriver.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformElectricField.hh"
#include "G4DormandPrince745.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include <geodetic_converter.hh>
#include <G4ExtendedMaterial.hh>

class G4UserLimits;

using namespace std;

class TGFDetectorConstruction : public G4VUserDetectorConstruction {
public:

	TGFDetectorConstruction();

	~TGFDetectorConstruction() override;

	G4VPhysicalVolume *Construct() override;


private:

	Settings *settings = Settings::getInstance();

	std::vector<SensitiveDet *> sens_det_List;

	std::vector<G4Material *> Construct_Atmos_layers_Materials(const std::vector<G4double> altitudes_);

	//    void ConstructAtmosMats2();
	//    void ConstructAtmosMats3();
	//    void ReadInputAtmosFile();

	//    double interp1(vector < double >,
	//                   vector < double >,
	//                   double);
	//    int    findNearestNeighbourIndex(double,
	//                                     vector < double >);

	void Construct_MagField_Managers();

	G4LogicalVolume *logicalWorld;
	G4VPhysicalVolume *physicalWorld;
	G4Material *vac = nullptr;

	G4Material *N2 = nullptr;
	G4Material *O2 = nullptr;
	G4Material *O = nullptr;
	G4Material *N = nullptr;
	G4Material *H = nullptr;

	std::vector<G4Sphere *> atmosLayers_S;
	std::vector<G4LogicalVolume *> atmosLayers_LV;
	std::vector<G4VPhysicalVolume *> atmosLayers_PV;

	std::vector<G4double> altitudes_geodetic; // (geodetic)altitudes intervals of the layers
	G4int nb_altitudes = 256;
	G4double alt_min = 1. * km; // (geodetic)
	G4double alt_max_atmosphere = 150. * km; // maximum altitude where the atmosphere is not negigible (with margin)
	G4double alt_account_Mag_Field = 45. * km; // altitude above which earth's magnetic field is takem into account

	// RK: the difference between geodetic and geographic altitude is small and probabli negligible

	G4double world_max_altitude = 15000. * km;

	// magnetic field managers
	G4FieldManager *globalfieldMgr = nullptr;
	G4FieldManager *Null_FieldManager = nullptr;
	G4ChordFinder *fChordFinder = nullptr;
	G4UniformMagField *magField_null = nullptr;

	G4double fMinStep = 0.010 * mm;
	G4double minEps = 1.0e-6; //   Minimum & value for smallest steps
	G4double maxEps = 1.0e-5; //   Maximum & value for largest steps

	// regions
	G4Region *considered_atmos_Region = new G4Region("considered_atmos");

	void calculate_altitudes_list();

	G4double maxStep = settings->STEP_MAX_VAL;
	G4UserLimits *stepLimit = new G4UserLimits(maxStep);

	G4bool not_contains(G4double value, const std::vector<G4double> &vec);

	G4MagIntegratorStepper *fStepper = nullptr;
	G4Mag_UsualEqRhs *pMagFldEquation = nullptr;

	// EarthMagField_alt *myEarthMagField = nullptr;
	G4MagneticField *myCachedEarthMagField = nullptr;

	std::ofstream asciiFile;

	bool hasDuplicates(const std::vector<G4double> &arr);

};
