/*
 * HydroParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>

#include "../include/HydroParameters.h"
#include "../include/Properties.h"

double initialProperTimePoint;
double shearViscosityToEntropyDensity;
double freezeoutTemperatureGeV;
int initializePimunuNavierStokes;

//requires libconfig
/*
void loadHydroParameters(config_t *cfg, const char* configDirectory, void * params) {
	// Read the file
	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "hydro.properties");
	if (!config_read_file(cfg, fname)) {
		fprintf(stderr, "No configuration file  %s found for hydrodynamic parameters - %s.\n", fname, config_error_text(cfg));
		fprintf(stderr, "Using default hydrodynamic configuration parameters.\n");
	}

	getDoubleProperty(cfg, "initialProperTimePoint", &initialProperTimePoint, 0.1);
	getDoubleProperty(cfg, "shearViscosityToEntropyDensity", &shearViscosityToEntropyDensity, 0.0795775);
	getDoubleProperty(cfg, "freezeoutTemperatureGeV", &freezeoutTemperatureGeV, 0.155);

	getIntegerProperty(cfg, "initializePimunuNavierStokes", &initializePimunuNavierStokes, 1);

	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->initialProperTimePoint = initialProperTimePoint;
	hydro->shearViscosityToEntropyDensity = shearViscosityToEntropyDensity;
	hydro->freezeoutTemperatureGeV = freezeoutTemperatureGeV;
	hydro->initializePimunuNavierStokes = initializePimunuNavierStokes;
}
*/

//does not require libconfig
void readHydroParameters(const char* configDirectory, void * params) {
	// Read the file
	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "hydro.properties");
	// std::ifstream is RAII, i.e. no need to call close
	std::ifstream cFile (fname);
	if (cFile.is_open())
	{
		std::string line;

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initialProperTimePoint = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		shearViscosityToEntropyDensity = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		freezeoutTemperatureGeV = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initializePimunuNavierStokes = atoi(line.c_str());

	}
	else std::cerr << "No configuration file  %s found for hydro parameters\n";

	struct HydroParameters * hydro = (struct HydroParameters *) params;
	hydro->initialProperTimePoint = initialProperTimePoint;
	hydro->shearViscosityToEntropyDensity = shearViscosityToEntropyDensity;
	hydro->freezeoutTemperatureGeV = freezeoutTemperatureGeV;
	hydro->initializePimunuNavierStokes = initializePimunuNavierStokes;
}
