#ifndef FREEZEOUT_H_
#define FREEZEOUT_H_

#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"

double linearInterp4D(double x0, double x1, double x2, double x3,
                      double a0000, double a1000, double a0100, double a0010, double a0001,
                      double a1100, double a1010, double a1001,
                      double a0110, double a0101, double a0011,
                      double a1110, double a1101, double a0111, double a1011, double a1111);

void swapAndSetHydroVariables(double ****energy_density_evolution, double *****hydrodynamic_evolution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ);

void setHydroVariables(double ****energy_density_evolution, double *****hydrodynamic_evolution,
                              CONSERVED_VARIABLES * const __restrict__ q, PRECISION * const __restrict__ e,
                              FLUID_VELOCITY * const __restrict__ u, int nx, int ny, int nz, int FOFREQ, int n);

void writeEnergyDensityToHypercube(double ****hyperCube, double ****energy_density_evolution, int it, int ix, int iy, int iz);

double interpolateVariable(double *****hydrodynamic_evolution, int ivar, int it, int ix, int iy, int iz, double tau_frac, double x_frac, double y_frac, double z_frac);

#endif
