/*
 * FluxLimiter.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */

#include <math.h>
#include <stdio.h>

#include "../include/FluxLimiter.h"
#include "../include/DynamicalVariables.h"

inline int sign(PRECISION x) {
	if (x<0) return -1;
	else return 1;
}

inline PRECISION minmod(PRECISION x, PRECISION y) {
	return (sign(x)+sign(y))*fmin(fabs(x),fabs(y))/2.0;
}

PRECISION minmod3(PRECISION x, PRECISION y, PRECISION z) {
   return minmod(x,minmod(y,z));
}

PRECISION approximateDerivative(PRECISION x, PRECISION y, PRECISION z) {
	PRECISION l = THETA * (y - x);
	PRECISION c = (z - x) / 2.0;
	PRECISION r = THETA * (z - y);
	return minmod3(l, c, r);
}
