/*
 * FluxFunctions.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: Dennis Bazow
 */

#include <stdlib.h>
#include <stdio.h> // for printf

#include "../include/FluxFunctions.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/DynamicalVariables.h"

PRECISION Fx(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return ux * q / ut;
}

PRECISION Fy(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return uy * q / ut;
}

PRECISION Fz(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un) {
	return un * q / ut;
}
