/*
 * InitialConditions.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef INITIALCONDITIONS_H_
#define INITIALCONDITIONS_H_

#include "../include/HydroInitialTmunu.h"

void setInitialConditions(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory);
void setICFromEnergyDensityVector(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, HydroInitialTmunu init_tmunu);
void setICFromPreequilVectors(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, HydroInitialTmunu init_tmunu);

#endif /* INITIALCONDITIONS_H_ */
