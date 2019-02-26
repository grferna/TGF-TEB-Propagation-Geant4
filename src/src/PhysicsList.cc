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
#include <PhysicsList.hh>
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4RadioactiveDecay.hh"

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TGF_PhysicsList::TGF_PhysicsList() :
        G4VUserPhysicsList() {
    emPhysicsList = new G4EmStandardPhysics_option1();

    this->DumpCutValuesTable();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TGF_PhysicsList::~TGF_PhysicsList() {
    delete emPhysicsList;
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructParticle() {
    emPhysicsList->ConstructParticle();

    //   G4GenericIon::GenericIon();
    //   G4NeutrinoE::NeutrinoEDefinition();
    //   G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    //   G4Alpha::AlphaDefinition();
    //
    //   G4Geantino::GeantinoDefinition();
    //   G4ChargedGeantino::ChargedGeantinoDefinition();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::ConstructProcess() {
    AddTransportation();
    emPhysicsList->ConstructProcess();

    //   G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

    //   radioactiveDecay->SetICM(true);                //Internal Conversion
    //   radioactiveDecay->SetARM(true);               //Atomic Rearangement

    // radioactiveDecay->SetFBeta(true);
    // radioactiveDecay->SetBRBias(true);

    //   G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    //
    //   ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TGF_PhysicsList::SetCuts() {
    defaultCutValue = 1. * nm;

    //
    cutForGamma = defaultCutValue;
    cutForElectron = defaultCutValue;
    cutForPositron = defaultCutValue;

    //
    //
    //
    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");

    // SetCutsWithDefault();

    //    G4EmParameters *param = G4EmParameters::Instance();
    //   param->SetMinEnergy(20.*keV);
    //    param->SetMaxEnergy(100.*MeV);
    //   param->SetLowestElectronEnergy(20.*keV);
    //    param->SetNumberOfBinsPerDecade(16);
    //   param->SetMscRangeFactor(0.01);
    //   param->SetFluo(true);
    //   param->SetAuger(true);
    //   param->SetPixe(true);
    //   param->SetPIXEElectronCrossSectionModel("Penelope");
    //    param->SetLateralDisplacement(true);
    //    param->ActivateAngularGeneratorForIonisation(true);
    //       G4double lowlimit=20.*CLHEP::keV;
    //       G4ProductionCutsTable * aPCTable = G4ProductionCutsTable::GetProductionCutsTable();
    //       aPCTable->SetEnergyRange(lowlimit,100*CLHEP::GeV);

    if (Settings::USE_STEP_MAX_for_record)
        Add_StepMax_for_record_regions();

    if (Settings::USE_STEP_MAX_global)
        AddStepMax(10.0 * m, 10.0 * m);
}

///////////////////////////////////////////////////////

#include "G4StepLimiter.hh"
#include "G4ProcessManager.hh"

void TGF_PhysicsList::Add_StepMax_for_record_regions() {
    // Step limitation seen as a process
    G4StepLimiter *stepLimiter = new G4StepLimiter();
    ////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();
        pmanager->AddDiscreteProcess(stepLimiter);
        ////pmanager ->AddDiscreteProcess(userCuts);
    }
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "StepMax.hh"

void TGF_PhysicsList::AddStepMax(G4double stepMax_elec, G4double stepMax_phot) {
    // Step limitation seen as a process
    StepMax *stepMaxProcess_elec = new StepMax();
    stepMaxProcess_elec->SetMaxStep(stepMax_elec);

    StepMax *stepMaxProcess_phot = new StepMax();
    stepMaxProcess_phot->SetMaxStep(stepMax_phot);

    auto particleIterator = GetParticleIterator();
    particleIterator->reset();

    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();

        if (stepMaxProcess_elec->IsApplicable(*particle)) {
            if (particle->GetParticleName() == "gamma") {
                pmanager->AddDiscreteProcess(stepMaxProcess_phot);
            } else {
                pmanager->AddDiscreteProcess(stepMaxProcess_elec);
            }
        }
    }
}
