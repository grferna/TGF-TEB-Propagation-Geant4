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

#include <Settings.hh>
#include <Analysis.hh>
#include <fstream>
#include <iomanip>

#include "G4Track.hh"
#include "G4SteppingManager.hh"

// singleton pattern initialization
Analysis *Analysis::instance = nullptr;

// constructor
Analysis::Analysis()
{

	const G4double ALT_MAX_RECORDED = *std::max_element(settings->record_altitudes.begin(), settings->record_altitudes.end());

	const G4String output_filename_second_part =
		std::to_string(settings->RANDOM_SEED) + "_" + std::to_string(int(ALT_MAX_RECORDED)) + "_" + std::to_string(int(settings->SOURCE_ALT)) + "_" + std::to_string(int(settings->OPENING_ANGLE)) + "_" + settings->BEAMING_TYPE + "_" +
		std::to_string(int(settings->SOURCE_SIGMA_TIME)) + ".out";

	asciiFileName_phot = "./output/detPhotons_" + output_filename_second_part;
	asciiFileName_lept = "./output/detLeptons_" + output_filename_second_part;

	std::ofstream asciiFile00(asciiFileName2, std::ios::trunc); // to clean the output file
	asciiFile00.close();

	output_lines_lept.clear();
	output_lines_phot.clear();

}

G4int Analysis::get_NB_RECORDED() const
{
	return NB_RECORDED;
}

Analysis::~Analysis() = default;

Analysis *Analysis::getInstance()
{
	if (instance == nullptr)
	{
		instance = new Analysis;
	}

	return instance;
}

void Analysis::save_in_output_buffer(const G4int PDG_NB, const G4double &time, const G4double &energy, const G4double &dist_rad, const G4int ID, const G4double &ecef_x, const G4double &ecef_y, const G4double &ecef_z)
{
	geod_conv::GeodeticConverter::ecef2Geodetic(ecef_x * 1000.0, ecef_y * 1000.0, ecef_z * 1000.0, lat, lon, alt);
	alt = alt / 1000.0;

	std::stringstream buffer;
	buffer << std::scientific << std::setprecision(7); // scientific notation with
	// 7 significant digits
	//   asciiFile << name;
	//   asciiFile << ' ';
	buffer << settings->RANDOM_SEED;
	buffer << ' ';
	buffer << settings->SOURCE_ALT;
	buffer << ' ';
	buffer << settings->OPENING_ANGLE;
	buffer << ' ';
	buffer << settings->TILT_ANGLE;
	buffer << ' ';
	buffer << settings->NB_EVENT;
	buffer << ' ';
	buffer << ID;
	buffer << ' ';
	buffer << PDG_NB;
	buffer << ' ';
	buffer << time;
	buffer << ' ';
	buffer << energy;
	buffer << ' ';
	buffer << alt;
	buffer << ' ';
	buffer << lat;
	buffer << ' ';
	buffer << lon;
	buffer << ' ';

	if (settings->OUTPUT_RadDist)
	{
		buffer << dist_rad;
	}

	if (settings->OUTPUT_ECEF_COORDS)
	{
		buffer << ecef_x;
		buffer << ' ';
		buffer << ecef_y;
		buffer << ' ';
		buffer << ecef_z;
	}

	NB_RECORDED++;

	//    buffer << ' ';
	//    buffer << creator_process;
	buffer << G4endl;

	if (PDG_NB == 22) {
		output_lines_phot.push_back(buffer.str());
	}
	else if (PDG_NB == -11 || PDG_NB == 11) {
		output_lines_lept.push_back(buffer.str());
	}

	write_in_output_file();
}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Analysis::write_in_output_file()
{
	if (output_lines_phot.size() >= output_buffer_size_phot)
	{

		std::ofstream asciiFile;
		asciiFile.open(asciiFileName_phot, std::ios::app);

		if (asciiFile.is_open())
		{
			for (G4String &line : output_lines_phot)
			{
				asciiFile << line;
			}
		}

		asciiFile.close();
		output_lines_phot.clear();

	}

	if (output_lines_lept.size() >= output_buffer_size_lept)
	{

		std::ofstream asciiFile;
		asciiFile.open(asciiFileName_lept, std::ios::app);

		if (asciiFile.is_open())
		{
			for (G4String &line : output_lines_lept)
			{
				asciiFile << line;
			}
		}

		asciiFile.close();
		output_lines_lept.clear();

	}

}

// ....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void Analysis::write_in_output_file_endOfRun()
{
	if (output_lines_phot.size() >= 1)
	{

		std::ofstream asciiFile;
		asciiFile.open(asciiFileName_phot, std::ios::app);

		if (asciiFile.is_open())
		{
			for (G4String &line : output_lines_phot)
			{
				asciiFile << line;
			}
		}

		asciiFile.close();
		output_lines_phot.clear();

	}

	if (output_lines_lept.size() >= 1)
	{

		std::ofstream asciiFile;
		asciiFile.open(asciiFileName_lept, std::ios::app);

		if (asciiFile.is_open())
		{
			for (G4String &line : output_lines_lept)
			{
				asciiFile << line;
			}
		}

		asciiFile.close();
		output_lines_lept.clear();

	}
}
