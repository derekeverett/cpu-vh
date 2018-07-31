#ifndef VORTICITY_H
#define VORTICITY_H

#include "DynamicalVariables.h"

void calculateThermalVorticity(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q,
void * latticeParams, void * hydroParams);

#endif
