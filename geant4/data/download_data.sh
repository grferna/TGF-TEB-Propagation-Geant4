#!/bin/bash
# This script OBSOLETE because data download and install is done with CMake

rm -rf *.tar.gz
rm -rf ./RealSurface2.1.1
rm -rf ./PhotonEvaporation5.2
rm -rf ./G4EMLOW7.3
rm -rf ./G4ABLA3.1
rm -rf ./RadioactiveDecay5.2
rm -rf ./G4ENSDFSTATE2.2
rm -rf ./G4NDL4.5
rm -rf ./G4NEUTRONXS1.4
rm -rf ./G4SAIDDATA1.1
rm -rf ./G4PII1.3

wget http://cern.ch/geant4-data/datasets/G4EMLOW.7.3.tar.gz
wget http://cern.ch/geant4-data/datasets/G4PhotonEvaporation.5.2.tar.gz
wget http://cern.ch/geant4-data/datasets/G4RadioactiveDecay.5.2.tar.gz
wget http://cern.ch/geant4-data/datasets/G4SAIDDATA.1.1.tar.gz
wget http://cern.ch/geant4-data/datasets/G4NEUTRONXS.1.4.tar.gz
wget http://cern.ch/geant4-data/datasets/G4ABLA.3.1.tar.gz
wget http://cern.ch/geant4-data/datasets/G4PII.1.3.tar.gz
wget http://cern.ch/geant4-data/datasets/G4RealSurface.2.1.1.tar.gz
wget http://cern.ch/geant4-data/datasets/G4ENSDFSTATE.2.2.tar.gz
wget http://cern.ch/geant4-data/datasets/G4NDL.4.5.tar.gz

for filename in *.tar.gz
do
  tar zxf $filename
done

rm -rf *.tar.gz

