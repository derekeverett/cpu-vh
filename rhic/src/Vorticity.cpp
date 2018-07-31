#include <stdio.h>
#include "../include/Vorticity.h"
#include "../include/LatticeParameters.h"
#include "../include/EquationOfState.h"

void calculateThermalVorticity(PRECISION t, PRECISION dt, CONSERVED_VARIABLES * __restrict__ q, CONSERVED_VARIABLES * __restrict__ Q,
void * latticeParams, void * hydroParams)
{

  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dy = (PRECISION)(lattice->latticeSpacingY);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

  for(int k = 2; k < ncz-2; ++k) {
		for(int j = 2; j < ncy-2; ++j) {
			for(int i = 2; i < ncx-2; ++i) {

				int s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
        int sxm = columnMajorLinearIndex(i-1, j, k, ncx, ncy, ncz);
        int sxp = columnMajorLinearIndex(i+1, j, k, ncx, ncy, ncz);
        int sym = columnMajorLinearIndex(i, j-1, k, ncx, ncy, ncz);
        int syp = columnMajorLinearIndex(i, j+1, k, ncx, ncy, ncz);
        int snm = columnMajorLinearIndex(i, j, k-1, ncx, ncy, ncz);
        int snp = columnMajorLinearIndex(i, j, k+1, ncx, ncy, ncz);

        //get the fluid velocity at neighbors
        PRECISION ut = u->ut[s];
        PRECISION ux = u->ux[s];
        PRECISION uy = u->uy[s];
        PRECISION un = u->un[s];

        PRECISION ut_p = up->ut[s];
        PRECISION ux_p = up->ux[s];
        PRECISION uy_p = up->uy[s];
        PRECISION un_p = up->un[s];

        PRECISION ut_xm = u->ut[sxm];
        PRECISION ux_xm = u->ux[sxm];
        PRECISION uy_xm = u->uy[sxm];
        PRECISION un_xm = u->un[sxm];

        PRECISION ut_xp = u->ut[sxp];
        PRECISION ux_xp = u->ux[sxp];
        PRECISION uy_xp = u->uy[sxp];
        PRECISION un_xp = u->un[sxp];

        PRECISION ut_ym = u->ut[sym];
        PRECISION ux_ym = u->ux[sym];
        PRECISION uy_ym = u->uy[sym];
        PRECISION un_ym = u->un[sym];

        PRECISION ut_yp = u->ut[syp];
        PRECISION ux_yp = u->ux[syp];
        PRECISION uy_yp = u->uy[syp];
        PRECISION un_yp = u->un[syp];

        PRECISION ut_nm = u->ut[snm];
        PRECISION ux_nm = u->ux[snm];
        PRECISION uy_nm = u->uy[snm];
        PRECISION un_nm = u->un[snm];

        PRECISION ut_np = u->ut[snp];
        PRECISION ux_np = u->ux[snp];
        PRECISION uy_np = u->uy[snp];
        PRECISION un_np = u->un[snp];

        //get the Temperature at neighbors
        PRECISION T = effectiveTemperature(e[s]);
        /*
        PRECISION T_p = effectiveTemperature(e[s]);
        PRECISION T_xm = effectiveTemperature(e[sxm]);
        PRECISION T_xp = effectiveTemperature(e[sxp]);
        PRECISION T_ym = effectiveTemperature(e[sym]);
        PRECISION T_yp = effectiveTemperature(e[syp]);
        PRECISION T_nm = effectiveTemperature(e[snm]);
        PRECISION T_np = effectiveTemperature(e[snp]);
        */
        //set the thermal velocity beta^\mu at neighbors

        PRECISION betat = ut / T;
        PRECISION betax = ux / T;
        PRECISION betay = uy / T;
        PRECISION betan = un / T;

        PRECISION betat_p = ut_p / T;
        PRECISION betax_p = ux_p / T;
        PRECISION betay_p = uy_p / T;
        PRECISION betan_p = un_p / T;

        PRECISION betat_xm = ut_xm / T;
        PRECISION betax_xm = ux_xm / T;
        PRECISION betay_xm = uy_xm / T;
        PRECISION betan_xm = un_xm / T;

        PRECISION betat_xp = ut_xp / T;
        PRECISION betax_xp = ux_xp / T;
        PRECISION betay_xp = uy_xp / T;
        PRECISION betan_xp = un_xp / T;

        PRECISION betat_ym = ut_ym / T;
        PRECISION betax_ym = ux_ym / T;
        PRECISION betay_ym = uy_ym / T;
        PRECISION betan_ym = un_ym / T;

        PRECISION betat_yp = ut_yp / T;
        PRECISION betax_yp = ux_yp / T;
        PRECISION betay_yp = uy_yp / T;
        PRECISION betan_yp = un_yp / T;

        PRECISION betat_nm = ut_nm / T;
        PRECISION betax_nm = ux_nm / T;
        PRECISION betay_nm = uy_nm / T;
        PRECISION betan_nm = un_nm / T;

        PRECISION betat_np = ut_np / T;
        PRECISION betax_np = ux_np / T;
        PRECISION betay_np = uy_np / T;
        PRECISION betan_np = un_np / T;

        //calculate time derivatives with backwards differences
        PRECISION ddt_betat = (betat_p - betat) / dt;
        PRECISION ddt_betax = (betax_p - betax) / dt;
        PRECISION ddt_betay = (betay_p - betay) / dt;
        PRECISION ddt_betan = (betan_p - betan) / dt;

        //calculate spatial derivatives using central differences
        PRECISION ddx_betat = (betat_xp - betat_xm) / (2.0 * dx);
        PRECISION ddy_betat = (betat_yp - betat_ym) / (2.0 * dy);
        PRECISION ddn_betat = (betat_np - betat_nm) / (2.0 * dz);

        PRECISION ddx_betax = (betax_xp - betax_xm) / (2.0 * dx);
        PRECISION ddy_betax = (betax_yp - betax_ym) / (2.0 * dy);
        PRECISION ddn_betax = (betax_np - betax_nm) / (2.0 * dz);

        PRECISION ddx_betay = (betay_xp - betay_xm) / (2.0 * dx);
        PRECISION ddy_betay = (betay_yp - betay_ym) / (2.0 * dy);
        PRECISION ddn_betay = (betay_np - betay_nm) / (2.0 * dz);

        PRECISION ddx_betan = (betan_xp - betan_xm) / (2.0 * dx);
        PRECISION ddy_betan = (betan_yp - betan_ym) / (2.0 * dy);
        PRECISION ddn_betan = (betan_np - betan_nm) / (2.0 * dz);

        //set the thermal vorticity tensor
        PRECISION wtx = (-1.0/2.0) *  ( ddt_betax - ddx_betat );
        PRECISION wty = (-1.0/2.0) *  ( ddt_betay - ddy_betat );
        PRECISION wtn = (-1.0/2.0) *  ( ddt_betan - ddn_betat );
        PRECISION wxy = (-1.0/2.0) *  ( ddx_betay - ddy_betax );
        PRECISION wxn = (-1.0/2.0) *  ( ddx_betan - ddn_betax );
        PRECISION wyn = (-1.0/2.0) *  ( ddy_betan - ddn_betay );

        //save values of vorticity tensor here to later write to FO surface
        
      } // for(int i = 2; i < ncx-2; ++i)
    } // for(int j = 2; j < ncy-2; ++j)
  } //for(int k = 2; k < ncz-2; ++k)

}
