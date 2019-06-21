/*
 * FiniteDifference.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: Dennis Bazow
 */

#include "../include/FiniteDifference.h"
#include "../include/DynamicalVariables.h"

PRECISION finiteDifferenceX(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx) {
	return 0.5 * (var[i+1 + NX * (j + NY * k)]-var[i-1 + NX * (j + NY * k)])/dx;
}

PRECISION finiteDifferenceY(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dx) {
	return 0.5 * (var[i+ NX * (j+1 + NY * k)]-var[i+ NX * (j-1 + NY * k)])/dx;
}

PRECISION finiteDifferenceZ(const PRECISION * const var, int i,int j,int k, int NX, int NY, int NZ, double dz) {
	return 0.5 * (var[i+ NX * (j + NY * (k+1))]-var[i+ NX * (j + NY * (k-1))])/dz;
}
