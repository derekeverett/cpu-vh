/*
 * EquationOfState.h
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#ifndef EQUATIONOFSTATE_H_
#define EQUATIONOFSTATE_H_

#include "DynamicalVariables.h"

#define CONFORMAL_EOS
//#define WUPERTAL_EOS
//#define HOTQCDHRG_EOS

// ideal gas of massless quarks and gluons
//#define EOS_FACTOR 15.6269 // Nc=3, Nf=3
#define EOS_FACTOR 13.8997 // Nc=3, Nf=2.5

//PRECISION equilibriumPressure(PRECISION e);
//PRECISION speedOfSoundSquared(PRECISION e);
//PRECISION effectiveTemperature(PRECISION e);
//PRECISION equilibriumEnergyDensity(PRECISION T);

class EOS
{
 private:

 public:

    double e_min;
    double e_max;
    double e_spacing;
    int e_length;

    double pressure_tb[100000];
    double temperature_tb[100000];

    PRECISION equilibriumPressure(PRECISION e);
    PRECISION speedOfSoundSquared(PRECISION e);
    PRECISION effectiveTemperature(PRECISION e);
    PRECISION equilibriumEnergyDensity(PRECISION T);

    double interpolate1D(double e, double *table);
    void loadEOSHotQCDHRG();
};

extern EOS eqnOfState;

#endif /* EQUATIONOFSTATE_H_ */
