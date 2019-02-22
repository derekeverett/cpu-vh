/*
 * HydroParameters.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef HYDROPARAMETERS_H_
#define HYDROPARAMETERS_H_

//#include <libconfig.h>

#include "DynamicalVariables.h"

struct HydroParameters
{
  double initialProperTimePoint;
  double shearViscosityToEntropyDensity;
  double freezeoutTemperatureGeV;
  int initializePimunuNavierStokes;
  double regulationTemperatureGeV;
  int doFreezeOut;
};

//requires libconfig
//void loadHydroParameters(config_t *cfg, const char* configDirectory, void * params);

//doesnt require libconfig
void readHydroParameters(const char* configDirectory, void * params);

#endif /* HYDROPARAMETERS_H_ */
