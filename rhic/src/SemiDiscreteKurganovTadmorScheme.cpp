/*
 * SemiDiscreteKurganovTadmorScheme.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "../include/SemiDiscreteKurganovTadmorScheme.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/DynamicalVariables.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/LocalPropagationSpeed.h"

void flux(const PRECISION * const __restrict__ data, PRECISION * const __restrict__ result,
		PRECISION (* const rightHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const leftHalfCellExtrapolation)(PRECISION qmm, PRECISION qm, PRECISION q, PRECISION qp, PRECISION qpp),
		PRECISION (* const spectralRadius)(PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION (* const fluxFunction)(PRECISION q, PRECISION ut, PRECISION ux, PRECISION uy, PRECISION un),
		PRECISION t, PRECISION ePrev
) {
	// left and right cells
	PRECISION qR[NUMBER_CONSERVED_VARIABLES], qL[NUMBER_CONSERVED_VARIABLES];

	// left and right extrapolated values of the conserved variables
	//int ptr = 0;
	//PRECISION qmm, qm, q, qp, qpp;
  #ifdef SIMD
  #pragma omp simd
  #endif
	for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
    int ptr = 5*n;
		PRECISION qmm 	= *(data+ptr);
		PRECISION qm 	= *(data+ptr+1);
		PRECISION q 		= *(data+ptr+2);
		PRECISION qp 	= *(data+ptr+3);
		PRECISION qpp 	= *(data+ptr+4);
		//ptr+=5;
		qR[n]	= rightHalfCellExtrapolation(qmm, qm, q, qp, qpp);
		qL[n]	= leftHalfCellExtrapolation(qmm, qm, q, qp, qpp);
	}

	// left and right extrapolated values of the primary variables
	PRECISION eR,pR,utR,uxR,uyR,unR;
	getInferredVariables(t,qR,ePrev,&eR,&pR,&utR,&uxR,&uyR,&unR);
	PRECISION eL,pL,utL,uxL,uyL,unL;
	getInferredVariables(t,qL,ePrev,&eL,&pL,&utL,&uxL,&uyL,&unL);

	//PRECISION a,qR_n,qL_n,FqR,FqL,res;
	PRECISION a = localPropagationSpeed(utR,uxR,uyR,unR,utL,uxL,uyL,unL,spectralRadius);
  #ifdef SIMD
  #pragma omp simd
  #endif
	for (unsigned int n = 0; n < NUMBER_CONSERVED_VARIABLES; ++n) {
		PRECISION qR_n = qR[n];
		PRECISION qL_n = qL[n];
		PRECISION FqR = fluxFunction(qR_n, utR, uxR, uyR, unR);
		PRECISION FqL = fluxFunction(qL_n, utL, uxL, uyL, unL);
		PRECISION res = FqR + FqL - a * (qR_n - qL_n);
		res /= 2;
		result[n] = res;
	}
}
