/*
 * MonteCarloGlauberModel.c
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include "../include/GlauberModel.h"
#include "../include/MonteCarloGlauberModel.h"
#include "../include/InitialConditionParameters.h"

#include <gsl/gsl_integration.h>

// generates A samples corresponding to Woods-Saxon nucleus with A nucleons
// results are loaded into x and y arrays (z is not used in current context)
// note: assumes that the the random number generator has already been seeded
void sampleWoodsSaxon(int A, double * const __restrict__ x, double * const __restrict__ y) {
	const double M = 1.01*4.5310551374155095;
	const double rmax = 20;
	int m = 0;
	double r,v,f,u;
	while (m<A) {
		r = rmax*((double) rand())/((double)RAND_MAX);
		u = ((double) rand())/((double)RAND_MAX);
		if (u < r*r*woodsSaxonDistribution(r,A)/M) {
			v = 2*(((double) rand())/((double)RAND_MAX)-0.5);
			f = 2*M_PI*((double) rand())/((double)RAND_MAX);
			x[m] = r*sqrt(1-v*v)*cos(f);
			y[m] = r*sqrt(1-v*v)*sin(f);
			m++;
		}
	}
}

int numberWoundedNucleons(int A, double b, double * const __restrict__ x, double * const __restrict__ y, double snn) {
	double x1[A], y1[A], x2[A], y2[A];
	int l1[A], l2[A];
	sampleWoodsSaxon(A,x1,y1);
	sampleWoodsSaxon(A,x2,y2);
	for (int i=0; i<A; i++) {
		x1[i] -= b/2;
		x2[i] += b/2;
	}
	double dn = sqrt(0.1*snn/M_PI);
	for(int i = 0; i < A; ++i) {
		l1[i] = 0;
		l2[i] = 0;
	}
	for(int i = 0; i < A; ++i) {
		for (int j = 0; j < A; ++j) {
			double dist = pow(x1[i]-x2[j],2);
			dist += pow(y1[i]-y2[j],2);
			dist = sqrt(dist);
			if (dist < dn) {
				l1[i] += 1;
				l2[j] += 1;
			}
		}
	}
	int n = 0;
	for(int i = 0; i < A; ++i) {
		if (l1[i] > 0) {
			x[n] = x1[i];
			y[n] = y1[i];
			n++;
		}
		if (l2[i] > 0) {
			x[n] = x2[i];
			y[n] = y2[i];
			n++;
		}
	}
	return n;
}

void
monteCarloGlauberEnergyDensityTransverseProfile(double * const __restrict__ energyDensityTransverse,
int nx, int ny, double dx, double dy, void * initCondParams
) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int NA = initCond->numberOfNucleonsPerNuclei;
	double b = initCond->impactParameter;
	double e0 = initCond->initialEnergyDensity;
	double snn = initCond->scatteringCrossSectionNN;
	double etaFlat = initCond->rapidityMean;
	double etaVariance = initCond->rapidityVariance;
	//the nucleon width parameter
	double SIG0 = 0.46;
	//double SIG0 = 0.66;

	double xp[2*NA], yp[2*NA];
   srand(1328398221);
	int nNucleons = numberWoundedNucleons(NA,b,xp,yp,snn);
	printf("==> Found %d wounded nucleons.\n", nNucleons);
	for(int i = 0; i < nx; ++i) {
   	for(int j = 0; j < ny; ++j) {
      	energyDensityTransverse[i + nx*j] = 0;
      }
	}
	for(int i = 0; i < nx; ++i) {
   	for(int j = 0; j < ny; ++j) {
      	for (int n = 0; n < nNucleons; ++n) {
         	double x = (i - ((double)nx-1.)/2.)*dx - xp[n];
            double y = (j - ((double)ny-1.)/2.)*dy - yp[n];
						// assumes gaussian bump in density (gaussian energy deposition)
            energyDensityTransverse[i + nx*j] += exp(-x*x/2/SIG0/SIG0-y*y/2/SIG0/SIG0);
         }
		}
	}
}
