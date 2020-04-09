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
        //save the covariant components of thermal vorticity tensor \omega_{\mu\nu}
        PRECISION ut = (u->ut[s]) * -1.0;
        PRECISION ux = (u->ux[s]) * -1.0;
        PRECISION uy = (u->uy[s]) * -1.0;;
        PRECISION un = (u->un[s]) * -1.0;;

        PRECISION ut_p = (up->ut[s]) * -1.0;;
        PRECISION ux_p = (up->ux[s]) * -1.0;;
        PRECISION uy_p = (up->uy[s]) * -1.0;;
        PRECISION un_p = (up->un[s]) * -1.0;;

        PRECISION ut_xm = (u->ut[sxm]) * -1.0;;
        PRECISION ux_xm = (u->ux[sxm]) * -1.0;;
        PRECISION uy_xm = (u->uy[sxm]) * -1.0;;
        PRECISION un_xm = (u->un[sxm]) * -1.0;;

        PRECISION ut_xp = (u->ut[sxp]) * -1.0;;
        PRECISION ux_xp = (u->ux[sxp]) * -1.0;;
        PRECISION uy_xp = (u->uy[sxp]) * -1.0;;
        PRECISION un_xp = (u->un[sxp]) * -1.0;;

        PRECISION ut_ym = (u->ut[sym]) * -1.0;;
        PRECISION ux_ym = (u->ux[sym]) * -1.0;;
        PRECISION uy_ym = (u->uy[sym]) * -1.0;;
        PRECISION un_ym = (u->un[sym]) * -1.0;;

        PRECISION ut_yp = (u->ut[syp]) * -1.0;;
        PRECISION ux_yp = (u->ux[syp]) * -1.0;;
        PRECISION uy_yp = (u->uy[syp]) * -1.0;;
        PRECISION un_yp = (u->un[syp]) * -1.0;;

        PRECISION ut_nm = (u->ut[snm]) * -1.0;;
        PRECISION ux_nm = (u->ux[snm]) * -1.0;;
        PRECISION uy_nm = (u->uy[snm]) * -1.0;;
        PRECISION un_nm = (u->un[snm]) * -1.0;;

        PRECISION ut_np = (u->ut[snp]) * -1.0;;
        PRECISION ux_np = (u->ux[snp]) * -1.0;;
        PRECISION uy_np = (u->uy[snp]) * -1.0;;
        PRECISION un_np = (u->un[snp]) * -1.0;;

        //get the Temperature at neighbors
        PRECISION T = eqnOfState.effectiveTemperature(e[s]);
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

        beta_mu->beta_t[s] = betat;
        beta_mu->beta_x[s] = betax;
        beta_mu->beta_y[s] = betay;
        beta_mu->beta_n[s] = betan;  

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
        PRECISION ddt_betat = (betat - betat_p) / dt;
        PRECISION ddt_betax = (betax - betax_p) / dt;
        PRECISION ddt_betay = (betay - betay_p) / dt;
        PRECISION ddt_betan = (betan - betan_p) / dt;

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
        wmunu->wtx[s] = (-1.0/2.0) *  ( ddt_betax - ddx_betat );
        wmunu->wty[s] = (-1.0/2.0) *  ( ddt_betay - ddy_betat );
        wmunu->wtn[s] = (-1.0/2.0) *  ( ddt_betan - ddn_betat );
        wmunu->wxy[s] = (-1.0/2.0) *  ( ddx_betay - ddy_betax );
        wmunu->wxn[s] = (-1.0/2.0) *  ( ddx_betan - ddn_betax );
        wmunu->wyn[s] = (-1.0/2.0) *  ( ddy_betan - ddn_betay );

      } // for(int i = 2; i < ncx-2; ++i)
    } // for(int j = 2; j < ncy-2; ++j)
  } //for(int k = 2; k < ncz-2; ++k)

}
