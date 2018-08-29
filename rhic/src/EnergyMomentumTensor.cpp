/*
 * EnergyMomentumTensor.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <math.h> // for math functions
#include <stdio.h>
#include <cstdlib>

#include "../include/EnergyMomentumTensor.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"

#include "../include/FullyDiscreteKurganovTadmorScheme.h" // for const params
#include "../include/EquationOfState.h"

#define MAX_ITERS_NR 10000 //Maximum number of newton-raphson iterations
#define MAX_ITERS_BI 1000000 //Maximum number of newton-raphson iterations
const PRECISION ACC = 1.0e-3; //Maximum relative error in energy density

PRECISION energyDensityFromConservedVariables(PRECISION ePrev, PRECISION M0, PRECISION M, PRECISION Pi) {

	#ifndef CONFORMAL_EOS
	PRECISION e0 = ePrev;	// initial guess for energy density

	//NEWTON-RAPHSON METHOD
	for (int j = 0; j < MAX_ITERS_NR; ++j) {
		PRECISION p = equilibriumPressure(e0);
		PRECISION cs2 = speedOfSoundSquared(e0);
		PRECISION cst2 = p/e0;

		PRECISION A = M0 * (1.0 - cst2) + Pi;
		PRECISION B = M0 * (M0 + Pi) - M;
		PRECISION H = sqrtf( fabsf(A*A + 4*cst2*B) );
		PRECISION D = (A-H) / (2*cst2);

		PRECISION f = e0 + D;
		PRECISION fp = 1.0 - ( (cs2 - cst2) * (B + D*H - ( (cs2 - cst2)*cst2*D*M0)/e0) ) / (cst2*e0*H);
		if (isnan(f/fp)) {printf("f/fp is nan \n"); exit(-1);}
		PRECISION e = e0 - f/fp;
		if (fabsf(e - e0) <=  ACC * fabsf(e)) return e;
		e0 = e;
	}

	//BISECTION METHOD if NEWTON RAPHSON fails
	printf("NR failed to find root. Calling Bisection routine \n");
	PRECISION e_lo = 0.0;
	PRECISION e_mid = ePrev;
	PRECISION e_hi = 2.0 * (ePrev + 1.0);

	for (int j = 0; j < MAX_ITERS_BI; ++j) {
		PRECISION p_lo = equilibriumPressure(e_lo);
		PRECISION p_mid = equilibriumPressure(e_mid);
		PRECISION p_hi = equilibriumPressure(e_hi);

		PRECISION cs2_lo = speedOfSoundSquared(e_lo);
		PRECISION cs2_mid = speedOfSoundSquared(e_mid);
		PRECISION cs2_hi = speedOfSoundSquared(e_hi);

		PRECISION cst2_lo = p_lo/e_lo;
		PRECISION cst2_mid = p_mid/e_mid;
		PRECISION cst2_hi = p_hi/e_hi;

		PRECISION A_lo = M0 * (1.0 - cst2_lo) + Pi;
		PRECISION A_mid = M0 * (1.0 - cst2_mid) + Pi;
		PRECISION A_hi = M0 * (1.0 - cst2_hi) + Pi;

		PRECISION B = M0 * (M0 + Pi) - M;

		PRECISION H_lo = sqrtf( fabsf(A_lo*A_lo + 4*cst2_lo*B) );
		PRECISION H_mid = sqrtf( fabsf(A_mid*A_mid + 4*cst2_mid*B) );
		PRECISION H_hi = sqrtf( fabsf(A_hi*A_hi + 4*cst2_hi*B) );

		PRECISION D_lo = (A_lo-H_lo) / (2*cst2_lo);
		PRECISION D_mid = (A_mid-H_mid) / (2*cst2_mid);
		PRECISION D_hi = (A_hi-H_hi) / (2*cst2_hi);

		PRECISION f_lo = e_lo + D_lo;
		PRECISION f_mid = e_mid + D_mid;
		PRECISION f_hi = e_hi + D_hi;

		if (j == 0 && (f_lo * f_hi > 0.0)) printf("root does not lie within brackets for guess in bisection \n");

		//if we converge within tolerance
		if ( (e_hi - e_lo) < (ACC * e_hi) / 2.0 ) return e_mid;
		else
		{
			//if root lies between lo and mid
			if (f_lo * f_mid < 0.0)
			{
				e_lo = e_lo;
				e_hi = e_mid;
				e_mid = (e_lo + e_hi) / 2.0;
			}
			//if root lies between mid and high
			else if (f_mid * f_hi < 0.0)
			{
				e_lo = e_mid;
				e_hi = e_hi;
				e_mid = (e_lo + e_hi) / 2.0;
			}
		} //else
	} //for (int j = 0; j < MAX_ITERS; ++j)


	//if all root solving routines fail
	printf("Maximum number of iterations exceeded.\tePrev=%.3f,\tM0=%.3f,\t M=%.3f,\t Pi=%.3f\n",ePrev,M0,M,Pi);
	return ePrev;
	#else
	return fabsf( sqrtf( fabsf(4 * M0 * M0 - 3 * M) ) - M0 );
	#endif
}

void getInferredVariables(PRECISION t, const PRECISION * const __restrict__ q, PRECISION ePrev,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
PRECISION * const __restrict__ ut, PRECISION * const __restrict__ ux, PRECISION * const __restrict__ uy, PRECISION * const __restrict__ un
) {
	PRECISION ttt = q[0];
	PRECISION ttx = q[1];
	PRECISION tty = q[2];
	PRECISION ttn = q[3];
#ifdef PIMUNU
	PRECISION pitt = q[4];
	PRECISION pitx = q[5];
	PRECISION pity = q[6];
	PRECISION pitn = q[7];
#else
	PRECISION pitt = 0;
	PRECISION pitx = 0;
	PRECISION pity = 0;
	PRECISION pitn = 0;
#endif
	// \Pi
#ifdef PI
	PRECISION Pi = q[14];
#else
	PRECISION Pi = 0;
#endif

/****************************************************************************\
#ifndef IDEAL
	PRECISION pixx = q[8];
	PRECISION pixy = q[9];
	PRECISION pixn = q[10];
	PRECISION piyy = q[11];
	PRECISION piyn = q[12];
	PRECISION pinn = q[13];
	PRECISION xi0 = (PRECISION)(0.05);
	PRECISION rhomax = (PRECISION)(0.7);
	PRECISION t2 = t*t;
	PRECISION pipi = pitt*pitt-2*pitx*pitx-2*pity*pity+pixx*pixx+2*pixy*pixy+piyy*piyy-2*pitn*pitn*t2+2*pixn*pixn*t2+2*piyn*piyn*t2+pinn*pinn*t2*t2;
	if(isnan(pipi)==1) printf("found pipi Nan\n");
	PRECISION spipi = sqrtf(fabsf(pipi+3*Pi*Pi));
	PRECISION pimumu = pitt - pixx - piyy - pinn*t*t;

	PRECISION pPrev = equilibriumPressure(ePrev);
	PRECISION a1 = spipi/rhomax/sqrtf(ePrev*ePrev+3*pPrev*pPrev);
	PRECISION a2 = pimumu/xi0/rhomax/spipi;
	PRECISION rho = fmaxf(a1,a2);
	PRECISION fac = tanh(rho)/rho;
	if(fabsf(rho)<1.e-7) fac = 1;
	pitt *= fac;
	pitx *= fac;
	pity *= fac;
	pitn *= fac;
#endif
/****************************************************************************/

PRECISION M0 = ttt - pitt;
PRECISION M1 = ttx - pitx;
PRECISION M2 = tty - pity;
PRECISION M3 = ttn - pitn;
PRECISION M = M1 * M1 + M2 * M2 + t * t * M3 * M3;

//Dennis' Scheme of regulating bulk pressure by hand
/*
#ifdef Pi
if ((M0 * M0 - M + M0 * Pi) < 0)
Pi = M / M0 - M0;
#endif
*/

/****************************************************************************/
if (ePrev <= 0.1)
{
	*e = M0 - M / M0;
}
else
{
	*e = energyDensityFromConservedVariables(ePrev, M0, M, Pi);
}
//printf("pitt = %f, pitx = %f, pity = %f, pitn = %f\n", pitt, pitx, pity, pitn);
if (isnan(*e)) {
	printf("\n Found e nan inside getInferredVariables \n");
	printf("M0=%.3f,\t M1=%.3f,\t M2=%.3f,\t M3=%.3f\n", M0, M1, M2, M3);
	printf("ttt=%.3f,\t ttx=%.3f,\t tty=%.3f,\t ttn=%.3f\n", ttt, ttx, tty, ttn);
	printf("pitt=%.3f,\t pitx=%.3f,\t pity=%.3f,\t pitn=%.3f\n", pitt, pitx, pity, pitn);
	exit(-1);
}
	*p = equilibriumPressure(*e);
	if (*e < 1.e-7) {
		*e = 1.e-7;
		*p = 1.e-7;
	}

	PRECISION P = *p + Pi;
	PRECISION E = 1.0/(*e + P);

	*ut = sqrtf(fabsf((M0 + P) * E));
	//*ut = sqrtf( (M0 + P) * E );

	PRECISION E2 = E/(*ut);
	*ux = M1 * E2;
	*uy = M2 * E2;
	*un = M3 * E2;

	//what should we do when solution for flow velocity is too large?
	//MUSIC has revert_grid in this case...
	if ( *ut > 1.0e10 )
	{
		printf("\n found ut = %f in getInferredVariables\n", *ut);
		printf("M0 = %.9f , e = %.9f, p = %.9f, Pi = %.9f, E = %.9f\n", M0, *e, *p, Pi, E);
		printf("***REGULATING FLOW VELOCITY FOR THIS CELL! (ut -> 1.0, ui -> 0.0)*** \n");
		*ut = 1.0;
		*ux = 0.0;
		*uy = 0.0;
		*un = 0.0;
	}
}

void setInferredVariablesKernel(const CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p, FLUID_VELOCITY * const __restrict__ u,
PRECISION t, void * latticeParams
) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int ncx,ncy,ncz;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	ncz = lattice->numComputationalLatticePointsRapidity;

  #pragma omp parallel for collapse(3)
	#ifdef TILE
	#pragma unroll_and_jam
	#endif
	for(int k = 2; k < ncz-2; ++k) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int i = 2; i < ncx-2; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);

				PRECISION q_s[NUMBER_CONSERVED_VARIABLES],_e,_p,ut,ux,uy,un;
				q_s[0] = q->ttt[s];
				q_s[1] = q->ttx[s];
				q_s[2] = q->tty[s];
				q_s[3] = q->ttn[s];
#ifdef PIMUNU
				q_s[4] = q->pitt[s];
				q_s[5] = q->pitx[s];
				q_s[6] = q->pity[s];
				q_s[7] = q->pitn[s];
/****************************************************************************/
				q_s[8] = q->pixx[s];
				q_s[9] = q->pixy[s];
				q_s[10] = q->pixn[s];
				q_s[11] = q->piyy[s];
				q_s[12] = q->piyn[s];
				q_s[13] = q->pinn[s];
/****************************************************************************/
#endif
#ifdef PI
				q_s[14] = q->Pi[s];
#endif
				getInferredVariables(t,q_s,e[s],&_e,&_p,&ut,&ux,&uy,&un);
				e[s] = _e;
				p[s] = _p;
				u->ut[s] = ut;
				u->ux[s] = ux;
				u->uy[s] = uy;
				u->un[s] = un;
			}
		}
	}
}

//===================================================================
// Components of T^{\mu\nu} in (\tau,x,y,\eta_s)-coordinates
//===================================================================
PRECISION Ttt(PRECISION e, PRECISION p, PRECISION ut, PRECISION pitt) {
	return (e+p)*ut*ut-p+pitt;
}

PRECISION Ttx(PRECISION e, PRECISION p, PRECISION ut, PRECISION ux, PRECISION pitx) {
	return (e+p)*ut*ux+pitx;
}

PRECISION Tty(PRECISION e, PRECISION p, PRECISION ut, PRECISION uy, PRECISION pity) {
	return (e+p)*ut*uy+pity;
}

PRECISION Ttn(PRECISION e, PRECISION p, PRECISION ut, PRECISION un, PRECISION pitn) {
	return (e+p)*ut*un+pitn;
}

PRECISION Txx(PRECISION e, PRECISION p, PRECISION ux, PRECISION pixx) {
	return (e+p)*ux*ux+p+pixx;
}

PRECISION Txy(PRECISION e, PRECISION p, PRECISION ux, PRECISION uy, PRECISION pixy) {
	return (e+p)*ux*uy+pixy;
}

PRECISION Txn(PRECISION e, PRECISION p, PRECISION ux, PRECISION un, PRECISION pixn) {
	return (e+p)*ux*un+pixn;
}

PRECISION Tyy(PRECISION e, PRECISION p, PRECISION uy, PRECISION piyy) {
	return (e+p)*uy*uy+p+piyy;
}

PRECISION Tyn(PRECISION e, PRECISION p, PRECISION uy, PRECISION un, PRECISION piyn) {
	return (e+p)*uy*un+piyn;
}

PRECISION Tnn(PRECISION e, PRECISION p, PRECISION un, PRECISION pinn, PRECISION t) {
	return (e+p)*un*un+p/t/t+pinn;
}
