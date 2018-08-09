/*
 * InitialConditionParameters.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <string>

#include "../include/InitialConditionParameters.h"
#include "../include/Properties.h"

int initialConditionType;
int numberOfSourceFiles;
int numberOfNucleonsPerNuclei;
double initialEnergyDensity;
double scatteringCrossSectionNN;
double impactParameter;
double fractionOfBinaryCollisions;
// longitudinal energy density profile parameters
double rapidityVariance; // \sigma^{2}_{\eta}
double rapidityMean; // flat region around \ets_s = 0

//uses libconfig
/*
void loadInitialConditionParameters(config_t *cfg, const char* configDirectory, void * params) {
	// Read the file
	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "ic.properties");
	if (!config_read_file(cfg, fname)) {
		fprintf(stderr, "No configuration file  %s found for initial condition parameters - %s.\n", fname, config_error_text(cfg));
		fprintf(stderr, "Using default initial condition configuration parameters.\n");
	}

	getIntegerProperty(cfg, "initialConditionType", &initialConditionType, 2);
	getIntegerProperty(cfg, "numberOfNucleonsPerNuclei", &numberOfNucleonsPerNuclei, 208);
  getIntegerProperty(cfg, "numberOfSourceFiles", &numberOfSourceFiles, 0);

	getDoubleProperty(cfg, "initialEnergyDensity", &initialEnergyDensity, 1.0);
	getDoubleProperty(cfg, "scatteringCrossSectionNN", &scatteringCrossSectionNN, 62);
	getDoubleProperty(cfg, "impactParameter", &impactParameter, 7);
	getDoubleProperty(cfg, "fractionOfBinaryCollisions", &fractionOfBinaryCollisions, 0.5);
	getDoubleProperty(cfg, "rapidityVariance", &rapidityVariance, 0.5);
	getDoubleProperty(cfg, "rapidityMean", &rapidityMean, 0.5);

	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;
	initCond->initialConditionType = initialConditionType;
  initCond->numberOfSourceFiles = numberOfSourceFiles;
	initCond->numberOfNucleonsPerNuclei = numberOfNucleonsPerNuclei;
	initCond->initialEnergyDensity = initialEnergyDensity;
	initCond->scatteringCrossSectionNN = scatteringCrossSectionNN;
	initCond->impactParameter = impactParameter;
	initCond->fractionOfBinaryCollisions = fractionOfBinaryCollisions;
	initCond->rapidityVariance = rapidityVariance;
	initCond->rapidityMean = rapidityMean;
}
*/

//does not require libconfig
void readInitialConditionParameters(const char* configDirectory, void * params) {
	// Read the file
	char fname[255];
	sprintf(fname, "%s/%s", configDirectory, "ic.properties");
	// std::ifstream is RAII, i.e. no need to call close
	std::ifstream cFile (fname);
	if (cFile.is_open())
	{
		std::string line;

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		auto delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initialConditionType = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		numberOfSourceFiles = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		numberOfNucleonsPerNuclei = atoi(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		initialEnergyDensity = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		scatteringCrossSectionNN = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		impactParameter = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		fractionOfBinaryCollisions = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		rapidityVariance = atof(line.c_str());

		getline(cFile, line);
		line.erase(remove_if(line.begin(), line.end(), isspace), line.end());
		delimiterPos = line.find("=");
		line = line.substr(delimiterPos + 1);
		rapidityMean = atof(line.c_str());
	}
	else std::cerr << "No configuration file  %s found for initial condition parameters\n";

	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) params;
	initCond->initialConditionType = initialConditionType;
	initCond->numberOfSourceFiles = numberOfSourceFiles;
	initCond->numberOfNucleonsPerNuclei = numberOfNucleonsPerNuclei;
	initCond->initialEnergyDensity = initialEnergyDensity;
	initCond->scatteringCrossSectionNN = scatteringCrossSectionNN;
	initCond->impactParameter = impactParameter;
	initCond->fractionOfBinaryCollisions = fractionOfBinaryCollisions;
	initCond->rapidityVariance = rapidityVariance;
	initCond->rapidityMean = rapidityMean;

}
