/*
 * InitialConditions.c
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP

#include "../include/InitialConditions.h"
#include "../include/HydroInitialTmunu.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/GlauberModel.h"
#include "../include/MonteCarloGlauberModel.h"
#include "../include/HydroParameters.h"
#include "../include/EquationOfState.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>
#define REGULATE 0 // 1 to regulate flow in dilute regions
#define GAMMA_MAX 10.0

#ifdef _OPENMP
#include <omp.h>
#endif

#define THETA_FUNCTION(X) ((double)X < (double)0 ? (double)0 : (double)1)

//*********************************************************************************************************\
//* Read in all initial profiles from a single or seperate file
//*********************************************************************************************************/

//this reads all hydro variables from a single file; this way we do not need to fetch the coordinates many times
//note that the file must contain values for all dissipative currents, even if they are zero !!!
void setInitialTmunuFromFile(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    float x, y, z, e_in, p_in, ut_in, ux_in, uy_in, un_in;
    //#ifdef PIMUNU
    float pitt_in, pitx_in, pity_in, pitn_in, pixx_in, pixy_in, pixn_in, piyy_in, piyn_in, pinn_in;
    //#endif
    //#ifdef PI
    float Pi_in;
    //#endif
    FILE *fileIn;
    char fname[255];

    sprintf(fname, "%s/%s", rootDirectory, "/input/Tmunu.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open Tmunu.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", &x, &y, &z, &e_in, &p_in, &ut_in, &ux_in, &uy_in, &un_in, &pitt_in, &pitx_in, &pity_in, &pitn_in, &pixx_in, &pixy_in, &pixn_in, &piyy_in, &piyn_in, &pinn_in, &Pi_in);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    e[s] =  (PRECISION) e_in;
                    //ep[s] = (PRECISION) e_in; //set previous step to same value
                    p[s] = p_in;
                    u->ut[s] = ut_in;
                    u->ux[s] = ux_in;
                    u->uy[s] = uy_in;
                    u->un[s] = un_in;
                    up->ut[s] = ut_in; //set previous step to same value
                    up->ux[s] = ux_in; //...
                    up->uy[s] = uy_in;
                    up->un[s] = un_in;
#ifdef PIMUNU
                    q->pitt[s] = pitt_in;
                    q->pitx[s] = pitx_in;
                    q->pity[s] = pity_in;
                    q->pitn[s] = pitn_in;
                    q->pixx[s] = pixx_in;
                    q->pixy[s] = pixy_in;
                    q->pixn[s] = pixn_in;
                    q->piyy[s] = piyy_in;
                    q->piyn[s] = piyn_in;
                    q->pinn[s] = pinn_in;
#endif
#ifdef PI
                    q->Pi[s] = Pi_in;
#endif
                }
            }
        }
    }
    fclose(fileIn);
}

//this function reads a separate file for every hydrodynamic variable
void setInitialTmunuFromFiles(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    float x, y, z, value;
    FILE *fileIn;
    char fname[255];

    //energy density and pressure
    sprintf(fname, "%s/%s", rootDirectory, "/input/e.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open e.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    e[s] =  (PRECISION) value;
                    p[s] = equilibriumPressure( e[s] );
                    //ep[s] = (PRECISION) value;
                    //printf("e [ %d ] = %f\n", s, e[s]);
                }
            }
        }
    }
    fclose(fileIn);

    //pressure
    /*
    sprintf(fname, "%s/%s", rootDirectory, "/input/p.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open p.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    p[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
    */

    //ut
    sprintf(fname, "%s/%s", rootDirectory, "/input/ut.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ut.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    u->ut[s] =  (PRECISION) value;
                    up->ut[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //ux
    sprintf(fname, "%s/%s", rootDirectory, "/input/ux.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open ux.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    u->ux[s] =  (PRECISION) value;
                    up->ux[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //uy
    sprintf(fname, "%s/%s", rootDirectory, "/input/uy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open uy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    u->uy[s] =  (PRECISION) value;
                    up->uy[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //un
    sprintf(fname, "%s/%s", rootDirectory, "/input/un.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open un.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    u->un[s] =  (PRECISION) value;
                    up->un[s] = (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

#ifdef PIMUNU
    //pitt
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitt.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitt.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pitt[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pitx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pitx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pity
    sprintf(fname, "%s/%s", rootDirectory, "/input/pity.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pity.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pity[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pitn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pitn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pitn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pitn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pixx
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixx.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixx.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pixx[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pixy
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pixy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pixn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pixn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pixn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pixn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //piyy
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyy.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyy.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->piyy[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //piyn
    sprintf(fname, "%s/%s", rootDirectory, "/input/piyn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open piyn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->piyn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);

    //pinn
    sprintf(fname, "%s/%s", rootDirectory, "/input/pinn.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open pinn.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->pinn[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
#ifdef PI
    //bulk
    sprintf(fname, "%s/%s", rootDirectory, "/input/bulk.dat");
    fileIn = fopen(fname, "r");
    if (fileIn == NULL)
    {
        printf("Couldn't open bulk.dat!\n");
    }
    else
    {
        for(int i = 2; i < nx+2; ++i) {
            for(int j = 2; j < ny+2; ++j) {
                for(int k = 2; k < nz+2; ++k) {
                    fscanf(fileIn, "%f %f %f %f\n", &x, &y, &z, &value);
                    int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                    q->Pi[s] =  (PRECISION) value;
                }
            }
        }
    }
    fclose(fileIn);
#endif
}

/*********************************************************************************************************\
 * Set initial flow profile
 *		- u^\mu = (1, 0, 0, 0)
 * 	- No transverse flow (ux = uy = 0)
 *		- Longitudinal scaling flow (u_z = z/t, i.e. un = 0)
/*********************************************************************************************************/
void setFluidVelocityInitialCondition(void * latticeParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double t0 = hydro->initialProperTimePoint;

  #pragma omp parallel for collapse(3)
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				PRECISION ux = 0.;
				PRECISION uy = 0.;
				PRECISION un = 0.;
				u->ux[s] = 0.;
				u->uy[s] = 0.;
				u->un[s] = 0.;
				u->ut[s] = sqrt(1.0+ux*ux+uy*uy+t0*t0*un*un);

        up->ux[s] = 0.;
				up->uy[s] = 0.;
				up->un[s] = 0.;
				up->ut[s] = sqrt(1.0+ux*ux+uy*uy+t0*t0*un*un);

        uS->ux[s] = 0.;
				uS->uy[s] = 0.;
				uS->un[s] = 0.;
				uS->ut[s] = sqrt(1.0+ux*ux+uy*uy+t0*t0*un*un);
			}
		}
	}
}

/*********************************************************************************************************\
 * Set initial shear-stress tensor \pi^\mu\nu
 *		- Navier-Stokes value, i.e. \pi^\mu\nu = 2 * (\epsilon + P) / T * \eta/S * \sigma^\mu\nu
 * 	- No initial pressure anisotropies (\pi^\mu\nu = 0)
/*********************************************************************************************************/
void setPimunuNavierStokesInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	PRECISION dx = (PRECISION)(lattice->latticeSpacingX);
	PRECISION dz = (PRECISION)(lattice->latticeSpacingRapidity);

	PRECISION etabar = (PRECISION)(hydro->shearViscosityToEntropyDensity);
	PRECISION t = hydro->initialProperTimePoint;

	PRECISION e0 = initCond->initialEnergyDensity;

  #pragma omp parallel for collapse(3)
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
//				double T = pow(e[s]/e0, 0.25);
				PRECISION T = effectiveTemperature(e[s]);
				if (T == 0) T = 1.e-3;
				//PRECISION pinn = -2/(3*t*t*t)*etabar*(e[s]+p[s])/T; //wrong by factor of 2
				PRECISION pinn = -4.0/(3.0*t*t*t)*etabar*(e[s] + p[s]) / T;
#ifdef PIMUNU
				q->pitt[s] = 0;
				q->pitx[s] = 0;
				q->pity[s] = 0;
				q->pitn[s] = 0;
				q->pixx[s] = -t*t*pinn/2;
				q->pixy[s] = 0;
				q->pixn[s] = 0;
				q->piyy[s] = -t*t*pinn/2;
				q->piyn[s] = 0;
				q->pinn[s] = pinn;
#endif
#ifdef PI
#define A_1 -13.77
#define A_2 27.55
#define A_3 13.45

#define LAMBDA_1 0.9
#define LAMBDA_2 0.25
#define LAMBDA_3 0.9
#define LAMBDA_4 0.22

#define SIGMA_1 0.025
#define SIGMA_2 0.13
#define SIGMA_3 0.0025
#define SIGMA_4 0.022
				PRECISION x = T/1.01355;
				PRECISION zetabar = A_1*x*x + A_2*x - A_3;
				if(x > 1.05)
					zetabar = LAMBDA_1*exp(-(x-1)/SIGMA_1) + LAMBDA_2*exp(-(x-1)/SIGMA_2)+0.001;
				else if(x < 0.995)
					zetabar = LAMBDA_3*exp((x-1)/SIGMA_3)+ LAMBDA_4*exp((x-1)/SIGMA_4)+0.03;
				q->Pi[s] = -zetabar*(e[s]+p[s])/T/t;
#endif
			}
		}
	}
}

void setPimunuInitialCondition(void * latticeParams, void * initCondParams, void * hydroParams) {
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
	int initializePimunuNavierStokes = hydro->initializePimunuNavierStokes;
	if (initializePimunuNavierStokes==1) {
		printf("Initialize \\pi^\\mu\\nu to its asymptotic Navier-Stokes value.\n");
#ifdef PI
		printf("Initialize \\Pi to its asymptotic Navier-Stokes value.\n");
#endif
		setPimunuNavierStokesInitialCondition(latticeParams, initCondParams, hydroParams);
		return;
	}
	else {
		printf("Initialize \\pi^\\mu\\nu to zero.\n");
		struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
		int nx = lattice->numLatticePointsX;
		int ny = lattice->numLatticePointsY;
		int nz = lattice->numLatticePointsRapidity;

    #pragma omp parallel for collapse(3)
		for(int i = 2; i < nx+2; ++i) {
			for(int j = 2; j < ny+2; ++j) {
				for(int k = 2; k < nz+2; ++k) {
					int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
#ifdef PIMUNU
			  		q->pitt[s] = 0;
			  		q->pitx[s] = 0;
			  		q->pity[s] = 0;
			  		q->pitn[s] = 0;
			  		q->pixx[s] = 0;
			  		q->pixy[s] = 0;
			  		q->pixn[s] = 0;
			  		q->piyy[s] = 0;
			  		q->piyn[s] = 0;
			  		q->pinn[s] = 0;
#endif
#ifdef PI
			  		q->Pi[s] = 0;
#endif
				}
			}
		}
		return;
	}
}


//*********************************************************************************************************\
//* Initial conditions for hydro with dynamical sources
//*********************************************************************************************************/
void setICfromSource(void * latticeParams, void * initCondParams, void * hydroParams, const char * rootDirectory){

    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
    struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

    int sourceType = initCond->sourceType;
    
    double initialEnergyDensity = initCond->initialEnergyDensity;

    int ncx = lattice->numComputationalLatticePointsX;
    int ncy = lattice->numComputationalLatticePointsY;
    int ncz = lattice->numComputationalLatticePointsRapidity;

    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    double dx = lattice->latticeSpacingX;
    double dy = lattice->latticeSpacingY;
    double dz = lattice->latticeSpacingRapidity;
    double t0 = hydro->initialProperTimePoint;
    
    int DIM_X = lattice->numLatticePointsX;
    int DIM_Y = lattice->numLatticePointsY;
    int DIM_ETA = lattice->numLatticePointsRapidity;
    int DIM = ncx*ncy*ncz;
    
    double DX = lattice->latticeSpacingX;
    double DY = lattice->latticeSpacingY;
    double DETA = lattice->latticeSpacingRapidity;
    double TAU = hydro->initialProperTimePoint;

    //==================================================
    // Initialize to vacuum
    //==================================================
    
    switch(sourceType){
        case 0: {
            printf("initialized to be vaccum \n");
            
            double ed = initialEnergyDensity;
            double pd = equilibriumPressure(ed);
            
            //--------------------------------------------------
            // Initialize energy, baryon and pressure density, also
            // shear, bulk, flow velocity and baryon diffusion current
            //--------------------------------------------------
            
            printf("Initialize \\pi^\\mu\\nu to zero.\n");
            printf("Initialize \\nb^\\mu to zero.\n");
            
            for(int i = 2; i < nx+2; ++i) {
                for(int j = 2; j < ny+2; ++j) {
                    for(int k = 2; k < nz+2; ++k) {
                        int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                        
                        e[s] = (PRECISION) ed;
                        p[s] = pd;
                        
                        PRECISION ux = 0;
                        PRECISION uy = 0;
                        PRECISION un = 0;
                        
                        u->ux[s] = 0;
                        u->uy[s] = 0;
                        u->un[s] = 0;
                        u->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
                        
                        up->ux[s] = 0;
                        up->uy[s] = 0;
                        up->un[s] = 0;
                        up->ut[s] = sqrt(1+ux*ux+uy*uy+t0*t0*un*un);
#ifdef PIMUNU
                        q->pitt[s] = 0;
                        q->pitx[s] = 0;
                        q->pity[s] = 0;
                        q->pitn[s] = 0;
                        q->pixx[s] = 0;
                        q->pixy[s] = 0;
                        q->pixn[s] = 0;
                        q->piyy[s] = 0;
                        q->piyn[s] = 0;
                        q->pinn[s] = 0;
#endif
#ifdef PI
                        q->Pi[s] = 0;
#endif
                    }
                }
            }
            return;
        }
            
        case 1: {
            
            printf("initialized with distributions in the Landau frame \n");
            
            //--------------------------------------------------
            // Solve eigenvalue problem in the Landau Frame
            //--------------------------------------------------
            
            double ttt_in, ttx_in, tty_in, ttn_in, txx_in, txy_in, txn_in, tyy_in, tyn_in, tnn_in;
            double stressTensor[10][DIM];
            
            FILE *fileIn;
            char fname[255];
            
            sprintf(fname, "%s/%s", rootDirectory, "input/Tmunu.dat");
            
            fileIn = fopen(fname, "r");
            if (fileIn == NULL)
            {
                printf("Couldn't open initialTmunu.dat!\n");
            }
            else
            {
                for(int i = 2; i < DIM_X+2; ++i) {
                    for(int j = 2; j < DIM_Y+2; ++j) {
                        for(int k = 2; k < DIM_ETA+2; ++k) {
                            fscanf(fileIn, "%*s%*s%*s%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &ttt_in, &ttx_in, &tty_in, &ttn_in, &txx_in, &txy_in, &txn_in, &tyy_in, &tyn_in, &tnn_in);
                            int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4, DIM_ETA+4);
                            stressTensor[0][is] = ttt_in;
                            stressTensor[1][is] = ttx_in;
                            stressTensor[2][is] = tty_in;
                            stressTensor[3][is] = ttn_in;
                            stressTensor[4][is] = txx_in;
                            stressTensor[5][is] = txy_in;
                            stressTensor[6][is] = txn_in;
                            stressTensor[7][is] = tyy_in;
                            stressTensor[8][is] = tyn_in;
                            stressTensor[9][is] = tnn_in;
                        }
                    }
                }
            }
            
            float tolerance = 1.0e-3; //set energy density to tolerance if it is less than tolerance and if REGULATE is true
            
            for (int is = 0; is < DIM; is++)
            {
                gsl_matrix *Tmunu; //T^(mu,nu) with two contravariant indices; we need to lower an index
                //using the metric to find the eigenvectors of T^(mu)_(nu) with one contravariant and one contravariant index
                Tmunu = gsl_matrix_alloc(4,4);
                gsl_matrix *gmunu;
                gmunu = gsl_matrix_alloc(4,4);
                gsl_matrix_complex *eigen_vectors;
                eigen_vectors = gsl_matrix_complex_alloc(4,4);
                gsl_vector_complex *eigen_values;
                eigen_values = gsl_vector_complex_alloc(4);
                //set the values of the energy momentum tensor
                gsl_matrix_set(Tmunu, 0, 0, stressTensor[0][is]); //tau,tau
                gsl_matrix_set(Tmunu, 0, 1, stressTensor[1][is]); //tau,x
                gsl_matrix_set(Tmunu, 0, 2, stressTensor[2][is]); //tau,y
                gsl_matrix_set(Tmunu, 0, 3, stressTensor[3][is]); //tau,eta
                gsl_matrix_set(Tmunu, 1, 1, stressTensor[4][is]); //x,x
                gsl_matrix_set(Tmunu, 1, 2, stressTensor[5][is]); //x,y
                gsl_matrix_set(Tmunu, 1, 3, stressTensor[6][is]); //x,eta
                gsl_matrix_set(Tmunu, 2, 2, stressTensor[7][is]); //y,y
                gsl_matrix_set(Tmunu, 2, 3, stressTensor[8][is]); //y,eta
                gsl_matrix_set(Tmunu, 3, 3, stressTensor[9][is]); //eta,eta
                gsl_matrix_set(Tmunu, 1, 0, stressTensor[1][is]); //x,tau
                gsl_matrix_set(Tmunu, 2, 0, stressTensor[2][is]); //y,tau
                gsl_matrix_set(Tmunu, 3, 0, stressTensor[3][is]); //eta,tau
                gsl_matrix_set(Tmunu, 2, 1, stressTensor[5][is]); //y,x
                gsl_matrix_set(Tmunu, 3, 1, stressTensor[6][is]); //eta,x
                gsl_matrix_set(Tmunu, 3, 2, stressTensor[8][is]); //eta,y
                
                //set the values of the "metric"; not really the metric, but the numerical constants
                //which are multiplied by the elements of T^(mu,nu) to get the values of T^(mu)_(nu)
                //note factors of TAU appropriate for milne coordinates g_(mu.nu) = diag(1,-1,-1,-TAU^2)
                gsl_matrix_set(gmunu, 0, 0, 1.0); //tau,tau
                gsl_matrix_set(gmunu, 0, 1, -1.0); //tau,x
                gsl_matrix_set(gmunu, 0, 2, -1.0); //tau,y
                gsl_matrix_set(gmunu, 0, 3, -1.0*TAU*TAU); //tau,eta
                gsl_matrix_set(gmunu, 1, 0, 1.0); //x,tau
                gsl_matrix_set(gmunu, 1, 1, -1.0); //x,x
                gsl_matrix_set(gmunu, 1, 2, -1.0); //x,y
                gsl_matrix_set(gmunu, 1, 3, -1.0*TAU*TAU); //x,eta
                gsl_matrix_set(gmunu, 2, 0, 1.0); //y,tau
                gsl_matrix_set(gmunu, 2, 1, -1.0); //y,x
                gsl_matrix_set(gmunu, 2, 2, -1.0); //y,y
                gsl_matrix_set(gmunu, 2, 3, -1.0*TAU*TAU); //y,eta
                gsl_matrix_set(gmunu, 3, 0, 1.0); //eta,tau
                gsl_matrix_set(gmunu, 3, 1, -1.0); //eta,x
                gsl_matrix_set(gmunu, 3, 2, -1.0); //eta,y
                gsl_matrix_set(gmunu, 3, 3, -1.0*TAU*TAU); //eta,eta
                //lower one index of the stress tensor; save it to the same matrix to save memory
                gsl_matrix_mul_elements(Tmunu, gmunu); //result stored in Tmunu !this multiplies element-wise, not ordinary matrix multiplication!
                gsl_eigen_nonsymmv_workspace *eigen_workspace;
                eigen_workspace = gsl_eigen_nonsymmv_alloc(4);
                gsl_eigen_nonsymmv(Tmunu, eigen_values, eigen_vectors, eigen_workspace);
                gsl_eigen_nonsymmv_free(eigen_workspace);
                
                //***does this have a solution for energy density and flow at every point?
                int eigenvalue_exists = 0;
                for (int i = 0; i < 4; i++)
                {
                    gsl_complex eigenvalue = gsl_vector_complex_get(eigen_values, i);
                    
                    if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
                    {
                        gsl_complex v0 = gsl_matrix_complex_get(eigen_vectors, 0 , i);
                        gsl_complex v1 = gsl_matrix_complex_get(eigen_vectors, 1 , i);
                        gsl_complex v2 = gsl_matrix_complex_get(eigen_vectors, 2 , i);
                        gsl_complex v3 = gsl_matrix_complex_get(eigen_vectors, 3 , i);
                        
                        if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
                        {
                            double minkowskiLength = GSL_REAL(v0)*GSL_REAL(v0) - (GSL_REAL(v1)*GSL_REAL(v1) + GSL_REAL(v2)*GSL_REAL(v2) + TAU*TAU*GSL_REAL(v3)*GSL_REAL(v3));
                            double factor = 1.0 / sqrt(minkowskiLength);
                            
                            if (GSL_REAL(v0) < 0) factor=-factor;
                            
                            //ignore eigenvectors with gamma >~ 60
                            if ( (GSL_REAL(v0) * factor) < GAMMA_MAX)
                            {
                                eigenvalue_exists = 1;
                                e[is] = GSL_REAL(eigenvalue);
                                u->ut[is] = GSL_REAL(v0) * factor;
                                u->ux[is] = GSL_REAL(v1) * factor;
                                u->uy[is] = GSL_REAL(v2) * factor;
                                u->un[is] = GSL_REAL(v3) * factor;
                            }
                            
                        } // if (GSL_IMAG(v0) == 0 && (2.0 * GSL_REAL(v0) * GSL_REAL(v0) - 1.0 - (GSL_REAL(v3) * GSL_REAL(v3) * (TAU * TAU - 1.0) )) > 0) //choose timelike eigenvector
                    } // if (GSL_REAL(eigenvalue) > 0.0 && GSL_IMAG(eigenvalue) == 0) //choose eigenvalue
                } //for (int i = 0; i < 4; ...)
                
                if (eigenvalue_exists == 0)
                {
                    //in dilute regions where we can't find a timelike eigenvector, set e = 0, u^t = 1, u^x=u^y=u^n=0
                    e[is] = 0.0;
                    u->ut[is] = 1.0;
                    u->ux[is] = 0.0;
                    u->uy[is] = 0.0;
                    u->un[is] = 0.0;
                }
            } // for (int is; is < DIM; ...)
            
            
            //--------------------------------------------------
            // Initialize flow velocity, energy, pressure and baryon densities
            //--------------------------------------------------
            
            //now regulate the flow velocity by a smoothing procedure. Flow can be too large in dilute regions, cause hydro to crash...
            //this method doesnt yield a smooth profile
            
            //try scaling the flow velocity by a smooth profile which goes to zero after some finite radius
            if (REGULATE)
            {
                printf("Regulating flow velocity profile in dilute regions \n");
                for(int i = 2; i < DIM_X+2; ++i) {
                    double x = (i-2 - (DIM_X-1)/2.)*DX;
                    for(int j = 2; j < DIM_Y+2; ++j) {
                        double y = (j-2 - (DIM_Y-1)/2.)*DY;
                        for(int k = 2; k < DIM_ETA+2; ++k) {
                            
                            int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4, DIM_ETA+4);
                            
                            double r = sqrt(x*x + y*y);
                            //printf("r=%f\n",r);
                            
                            float R_WIDTH = 0.6;
                            float R_FLAT = 4.5;
                            float arg = (-1.0) * (r - R_FLAT) * (r - R_FLAT) / (2.0 * R_WIDTH * R_WIDTH);
                            arg = arg * THETA_FUNCTION(r - R_FLAT);
                            
                            u->ux[is] = u->ux[is] * exp(arg);
                            u->uy[is] = u->uy[is] * exp(arg);
                            u->un[is] = 0.0;
                            
                            u->ut[is] = sqrt( 1 + u->ux[is]*u->ux[is] + u->uy[is]*u->uy[is] + TAU*TAU*u->un[is]*u->un[is]);
                        }
                    }
                }
            }
            
            for(int i = 2; i < DIM_X+2; ++i) {
                for(int j = 2; j < DIM_Y+2; ++j) {
                    for(int k = 2; k < DIM_ETA+2; ++k) {
                        int is = columnMajorLinearIndex(i, j, k, DIM_X+4, DIM_Y+4, DIM_ETA+4);
                        
                        u->ux[is] = 0.0;
                        u->uy[is] = 0.0;
                        u->un[is] = 0.0;
                        u->ut[is] = sqrt( 1 + u->ux[is]*u->ux[is] + u->uy[is]*u->uy[is] + TAU*TAU*u->un[is]*u->un[is]);
                        
                        e[is] = e[is] + 1.e-3;
                        up->ut[is] = u->ut[is];
                        up->ux[is] = u->ux[is];
                        up->uy[is] = u->uy[is];
                        up->un[is] = u->un[is];
                        p[is] = equilibriumPressure(e[is]);
                    }
                }
            }
            fclose(fileIn);
            
            //--------------------------------------------------
            // Initialize shear, bulk
            //--------------------------------------------------
            
            printf("Initialize \\pi^\\mu\\nu to zero.\n");
            
            int nx = lattice->numLatticePointsX;
            int ny = lattice->numLatticePointsY;
            int nz = lattice->numLatticePointsRapidity;
            for(int i = 2; i < nx+2; ++i) {
                for(int j = 2; j < ny+2; ++j) {
                    for(int k = 2; k < nz+2; ++k) {
                        int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
#ifdef PIMUNU
                        q->pitt[s] = 0;
                        q->pitx[s] = 0;
                        q->pity[s] = 0;
                        q->pitn[s] = 0;
                        q->pixx[s] = 0;
                        q->pixy[s] = 0;
                        q->pixn[s] = 0;
                        q->piyy[s] = 0;
                        q->piyn[s] = 0;
                        q->pinn[s] = 0;
#endif
#ifdef PI
                        q->Pi[s] = 0;
#endif
                    }
                }
            }
            return;
        }
    }
}

/*********************************************************************************************************\
 * Constant initial energy density distribution
/*********************************************************************************************************/
void setConstantEnergyDensityInitialCondition(void * latticeParams, void * initCondParams) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	double initialEnergyDensity = initCond->initialEnergyDensity;

	double T0 = 3.05;
	double ed = equilibriumEnergyDensity(T0);

	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

  #pragma omp parallel for collapse(3)
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the sound propagation test
 *		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/*********************************************************************************************************\
void setSoundPropagationInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	double initialEnergyDensity = initCond->initialEnergyDensity;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double T0 = 3.05; // -> 0.6 GeV
	double e0 = initialEnergyDensity*pow(T0, 4);

	double cs = 0.57735;
	double de = e0/100.;
	double PI = 3.141592653589793;
	double p0=e0/3;
	double lambda = (nx-1)*dx/2.;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		double vx = cs*de/(e0+p0)*sin(2*PI*x/lambda);
		double ed = e0 + de*sin(2*PI*x/lambda);
		// periodic boundary conditions
		if (i==2) {
			vx = cs*de/(e0+p0)*sin(2*PI*abs(x)/lambda);
			ed = e0 + de*sin(2*PI*abs(x)/lambda);
		}

		double u0 = 1/sqrt(1-vx*vx);

		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4);
				e[s] = ed;
				p[s] = Pressure(e[s]);
				u->ux[s] = u0*vx;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = u0;
				// initialize \pi^\mu\nu to zero
        		q->pitt[s] = 0;
        		q->pitx[s] = 0;
        		q->pity[s] = 0;
        		q->pitn[s] = 0;
        		q->pixx[s] = 0;
        		q->pixy[s] = 0;
        		q->pixn[s] = 0;
        		q->piyy[s] = 0;
        		q->piyn[s] = 0;
        		q->pinn[s] = 0;
			}
		}
	}
}

/*********************************************************************************************************\
 * Longitudinal initial energy density distribution
/*********************************************************************************************************/
void longitudinalEnergyDensityDistribution(double * const __restrict__ eL, void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nz = lattice->numLatticePointsRapidity;

	double dz = lattice->latticeSpacingRapidity;

	double etaFlat = initCond->rapidityMean;
	double etaVariance = initCond->rapidityVariance;

  #pragma omp parallel for
	for(int k = 0; k < nz; ++k) {
		double eta = (k - (nz-1)/2)*dz;
		double etaScaled = fabs(eta) - etaFlat/2;
		double arg = -etaScaled * etaScaled / etaVariance / 2 * THETA_FUNCTION(etaScaled);
		eL[k] = exp(arg);
	}
}

/*********************************************************************************************************\
 * Continuous optical glauber Glauber initial energy density distribution
/*********************************************************************************************************/
void setGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
	double T0 = 2.05;
//	e0 *= pow(T0,4);
	e0 = (double) equilibriumEnergyDensity(T0);

	double eT[nx*ny], eL[nz];
	energyDensityTransverseProfileAA(eT, nx, ny, dx, dy, initCondParams);
	longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);

  #pragma omp parallel for collapse(3)
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
        double energyDensityTransverse = e0 * eT[i-2+(j-2)*nx];
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				double energyDensityLongitudinal = eL[k-2];
				double ed = (energyDensityTransverse * energyDensityLongitudinal) + 1.e-3;
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Monte carlo Glauber initial energy density distribution
/*********************************************************************************************************/
void setMCGlauberInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;
//	double T0 = 3.05;
	double T0 = 2.03;
//	e0 *= pow(T0,4);
	e0 = (double) equilibriumEnergyDensity(T0);

	double eT[nx*ny], eL[nz];
	monteCarloGlauberEnergyDensityTransverseProfile(eT, nx, ny, dx, dy, initCondParams);
	longitudinalEnergyDensityDistribution(eL, latticeParams, initCondParams);

  #pragma omp parallel for collapse(3)
	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			for(int k = 2; k < nz+2; ++k) {
        double energyDensityTransverse = e0 * eT[i-2 + nx*(j-2)];
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				double energyDensityLongitudinal = eL[k-2];
				double ed = (energyDensityTransverse * energyDensityLongitudinal) + 1.e-3;
				e[s] = (PRECISION) ed;
				p[s] = equilibriumPressure(e[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the Gubser ideal hydro test
 *		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/*********************************************************************************************************/
void setIdealGubserInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			double T = 1.9048812623618392/pow(1 + pow(1 - pow(x,2) - pow(y,2),2) + 2*(1 + pow(x,2) + pow(y,2)),0.3333333333333333);
			double r = sqrt(x*x+y*y);
			double phi = atanh(2*1*r/(1+1+x*x+y*y));

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);

				e[s] = (PRECISION) (e0 * pow(T,4));
				p[s] = e[s]/3;
				u->ux[s] = (PRECISION) (sinh(phi)*x/r);
				u->uy[s] = (PRECISION) (sinh(phi)*y/r);
				u->un[s] = 0;
				u->ut[s] = sqrt(1 + u->ux[s]*u->ux[s] + u->uy[s]*u->uy[s]);
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the Gubser viscous hydro test
 *		- set energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny
/*********************************************************************************************************/
void setISGubserInitialCondition(void * latticeParams, const char *rootDirectory) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double x,y,ed,u1,u2,pitt,pitx,pity,pixx,pixy,piyy,pinn;

	FILE *file;
	char fname[255];
	sprintf(fname, "%s/%s", rootDirectory, "/rhic/rhic-trunk/src/test/resources/gubser/viscous/gubserIC.dat");
	file = fopen(fname, "r");

	double pitn=0;
	double pixn=0;
	double piyn=0;

	for(int i = 2; i < nx+2; ++i) {
		for(int j = 2; j < ny+2; ++j) {
			int status = fscanf(file,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		    		&x,&y,&ed,&u1,&u2,&pixx,&piyy,&pixy,&pitt,&pitx,&pity,&pinn);
			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);

				e[s] = (PRECISION) ed;
				p[s] = e[s]/3;
				u->ux[s] = u1;
				u->uy[s] = u2;
				u->un[s] = 0;
				u->ut[s] = sqrt(1 + u1*u1 + u2*u2);
#ifdef PIMUNU
        		q->pitt[s] = (PRECISION) pitt;
        		q->pitx[s] = (PRECISION) pitx;
        		q->pity[s] = (PRECISION) pity;
        		q->pitn[s] = (PRECISION) pitn;
        		q->pixx[s] = (PRECISION) pixx;
        		q->pixy[s] = (PRECISION) pixy;
        		q->pixn[s] = (PRECISION) pixn;
        		q->piyy[s] = (PRECISION) piyy;
        		q->piyn[s] = (PRECISION) piyn;
        		q->pinn[s] = (PRECISION) pinn;
#endif
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the relativistic Sod shock-tube test
 *		- set energy density, pressure, fluid velocity u^\mu
/*********************************************************************************************************/
void setSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				if(x > 0) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);
//				if(y > 0) 	e[s] = (PRECISION) (0.00778147);
//				else 			e[s] = (PRECISION) (0.124503);
//				if(x > 0) 	e[s] = (PRECISION) (1.0);
//				else 			e[s] = (PRECISION) (100.0);
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the relativistic Sod shock-tube test
 *		- set energy density, pressure, fluid velocity u^\mu
/*********************************************************************************************************/
void set2dSodShockTubeInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				if(y > x) 	e[s] = (PRECISION) (0.00778147);
//				if(atan(y/x)>0.7853981634) 	e[s] = (PRECISION) (0.00778147);
				else 			e[s] = (PRECISION) (0.124503);
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the implosion in a box test
 *		- set energy density, pressure, fluid velocity u^\mu
/*********************************************************************************************************/
void setImplosionBoxInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double e0 = initCond->initialEnergyDensity;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
//				e[s] = (PRECISION) (0.00778147);
//				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.124503);

//				e[s] = (PRECISION) (0.124503);
//				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.00778147);
				e[s] = (PRECISION) (0.00778147);
				if (sqrt(x*x+y*y)<=0.15) e[s] = (PRECISION) (0.124503);

//				e[s] = (PRECISION) (1.0);
/*
				if (x < 1) {
					if (y < (1-x))
						e[s] = (PRECISION) (0.00778147);
				}
				else e[s] = (PRECISION) (0.124503);
//*/
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the implosion in a box test
 *		- set energy density, pressure, fluid velocity u^\mu
/*********************************************************************************************************/
void setRayleighTaylorInstibilityInitialCondition(void * latticeParams, void * initCondParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;

	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	double Lx = ( (nx-1)/2.)*dx;
	double Ly = ( (ny-1)/2.)*dy;

	for(int i = 2; i < nx+2; ++i) {
		double x = (i-2 - (nx-1)/2.)*dx;
		for(int j = 2; j < ny+2; ++j) {
			double y = (j-2 - (ny-1)/2.)*dy;

			for(int k = 2; k < nz+2; ++k) {
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);

				double gasGamma = 1.4;
				double gravity = 0.1;

				double rhoTop = 2.0;
				double rhoBot = 1.0;
				double pr0 = 0.01;
				double pert = 0.01;

				double yloc = 0.5 + pert*cos(M_PI*x);
				double pr;
				if(y > yloc) 	pr = rhoTop*gravity*(1-y);
				else 				pr = rhoTop*gravity*(1-yloc) + rhoBot*gravity*(yloc-y);
				e[s] = pr;
				p[s] = e[s]/3;
				u->ux[s] = 0;
				u->uy[s] = 0;
				u->un[s] = 0;
				u->ut[s] = 1;
			}
		}
	}
}

/*********************************************************************************************************\
 * Initial conditions for the implosion in a box test
 *		- set energy density, pressure, fluid velocity u^\mu
/*********************************************************************************************************/
void setGaussianPulseInitialCondition(void * latticeParams, void * initCondParams) {
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

  int nx = lattice->numLatticePointsX;
  int ny = lattice->numLatticePointsY;
  int nz = lattice->numLatticePointsRapidity;

  double dx = lattice->latticeSpacingX;
  double dy = lattice->latticeSpacingY;
  double dz = lattice->latticeSpacingRapidity;

  double Lx = ( (nx-1)/2.)*dx;
  double Ly = ( (ny-1)/2.)*dy;

  for(int i = 2; i < nx+2; ++i) {
    double x = (i-2 - (nx-1)/2.)*dx;
    for(int j = 2; j < ny+2; ++j) {
      double y = (j-2 - (ny-1)/2.)*dy;

      for(int k = 2; k < nz+2; ++k) {
        int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);

        double xc = 0.5;
        double yc = 0.5;
        double beta = 50.0;
        double pr = 1.0 + 1e-1*exp(-beta*((x-xc)*(x-xc)+(y-yc)*(y-yc)));

        e[s] = (PRECISION) (pr/(1.4-1.));
        p[s] = e[s]/3.0;
        u->ux[s] = 0;
        u->uy[s] = (1+cos(2*M_PI*x/Lx))*(1+cos(2*M_PI*y/Ly));
        //u->uy[s] = 0;
        u->un[s] = 0;
        u->ut[s] = sqrt(1+u->uy[s]*u->uy[s]);
      }
    }
  }
}

/*********************************************************************************************************\
 * Initial conditions to use.
 *	Set the energy density, pressure, fluid velocity u^\mu, and \pi^\mu\ny.
 * 	0 - constant energy density
 *		1 - Isreal-Stewart hydrodynamic Gubser flow test
 *		2 - Continous optical Glauber
 *		3 - Ideal hydrodynamic Gubser flow test
 *		4 - Monte carlo Glauber
 *		5 - Relativistic Sod shock-tube test
/*********************************************************************************************************/
void setInitialConditions(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory) {
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	int initialConditionType = initCond->initialConditionType;
	printf("Setting initial conditions: ");
	switch (initialConditionType) {
		case 0: {
			printf("constant energy density.\n");
			setConstantEnergyDensityInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 1: {
			printf("Isreal-Stewart hydrodynamic Gubser flow test.\n");
			setISGubserInitialCondition(latticeParams, rootDirectory);
			return;
		}
		case 2: {
			printf("Continous optical Glauber.\n");
			setGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 3: {
			printf("Ideal hydrodynamic Gubser flow test.\n");
			setIdealGubserInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 4: {
			printf("Monte carlo Glauber.\n");
			setMCGlauberInitialCondition(latticeParams, initCondParams);
			setFluidVelocityInitialCondition(latticeParams, hydroParams);
			setPimunuInitialCondition(latticeParams, initCondParams, hydroParams);
			return;
		}
		case 5: {
			printf("Relativistic Sod shock-tube test.\n");
			setSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 6: {
			printf("Implosion in a box test.\n");
			setImplosionBoxInitialCondition(latticeParams, initCondParams);
			return;
		}

		case 7: {
			printf("Rayleigh-Taylor instability test.\n");
			setRayleighTaylorInstibilityInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 8: {
			printf("Implosion in a box test.\n");
			setGaussianPulseInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 9: {
			printf("Relativistic 2d Sod shock-tube test.\n");
			set2dSodShockTubeInitialCondition(latticeParams, initCondParams);
			return;
		}
		case 10: {
            printf("Reading initial T ^mu nu from input/e.dat , input/p.dat , etc... \n");
            setInitialTmunuFromFiles(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }
        case 11: {
            printf("Reading initial T ^mu nu from /input/Tmunu.dat \n");
            setInitialTmunuFromFile(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
		}
        case 12:{
            printf("Hydro with dynamical sources...\n");
            setICfromSource(latticeParams, initCondParams, hydroParams, rootDirectory);
            return;
        }

		default: {
			printf("Initial condition type not defined. Exiting ...\n");
			exit(-1);
		}
	}
}

void setICFromEnergyDensityVector(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, HydroInitialTmunu init_tmunu)
{
  printf("Initialize from initial energy density vector ...\n");
  const double hbarc = 0.197326938;
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
  struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
  int nx = lattice->numLatticePointsX;
  int ny = lattice->numLatticePointsY;
  int nz = lattice->numLatticePointsRapidity;

  for (int i = 2; i < nx+2; ++i) {
    for (int j = 2; j < ny+2; ++j) {
      for (int k = 2; k < nz+2; ++k) {
        //s runs over hydro grid with ghost cells
        int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
        //sm runs over preequilibrium grid without ghost cells
        int im = i-2;
        int jm = j-2;
        int km = k-2;
        int sm = columnMajorLinearIndex(im, jm, km, nx, ny, nz);

        e[s] =  (PRECISION) (init_tmunu.e_in[sm] / hbarc) + (PRECISION)(1.0e-3);
        p[s] = equilibriumPressure( e[s] );
        u->ut[s] = 1.0;
        u->ux[s] = 0.0;
        u->uy[s] = 0.0;
        u->un[s] = 0.0;
        up->ut[s] = 1.0; //set previous step to same value
        up->ux[s] = 0.0; //...
        up->uy[s] = 0.0;
        up->un[s] = 0.0;
        #ifdef PIMUNU
        q->pitt[s] = 0.0;
        q->pitx[s] = 0.0;
        q->pity[s] = 0.0;
        q->pitn[s] = 0.0;
        q->pixx[s] = 0.0;
        q->pixy[s] = 0.0;
        q->pixn[s] = 0.0;
        q->piyy[s] = 0.0;
        q->piyn[s] = 0.0;
        q->pinn[s] = 0.0;
        #endif
        #ifdef PI
        q->Pi[s] = 0.0;
        #endif
      }
    } // for (int j = 2; j < ny+2; ++j)
  } //for (int i = 2; i < nx+2; ++i)
}

void setICFromPreequilVectors(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, HydroInitialTmunu init_tmunu)
{
  printf("Initialize from preequilibrium vectors ...\n");
  const double hbarc = 0.197326938;
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
  struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;
  int nx = lattice->numLatticePointsX;
  int ny = lattice->numLatticePointsY;
  int nz = lattice->numLatticePointsRapidity;

  for (int i = 2; i < nx+2; ++i) {
    for (int j = 2; j < ny+2; ++j) {
      for (int k = 2; k < nz+2; ++k) {
        //s runs over hydro grid with ghost cells
        int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
        //sm runs over preequilibrium grid without ghost cells
        int im = i-2;
        int jm = j-2;
        int km = k-2;
        int sm = columnMajorLinearIndex(im, jm, km, nx, ny, nz);
        e[s] =  (PRECISION) (init_tmunu.e_in[sm] / hbarc) + (PRECISION)(1.0e-3);
        p[s] = equilibriumPressure( e[s] );

        PRECISION ux, uy, un, ut;
        if ( fabs(init_tmunu.ux_in[sm]) < 1.0e-7 ) ux = 0.0;
        else ux = (PRECISION) init_tmunu.ux_in[sm];
        if ( fabs(init_tmunu.uy_in[sm]) < 1.0e-7 ) uy = 0.0;
        else uy = (PRECISION) init_tmunu.uy_in[sm];
        if ( fabs(init_tmunu.un_in[sm]) < 1.0e-7 ) un = 0.0;
        else un = (PRECISION) init_tmunu.un_in[sm];

        ut = 1.0 + ux*ux + uy*uy + un*un;

        u->ux[s] = ux;
        u->uy[s] = uy;
        u->un[s] = un;
        u->ut[s] = ut;

        up->ut[s] = ut; //set previous step to same value
        up->ux[s] = ux; //...
        up->uy[s] = uy;
        up->un[s] = un;

        #ifdef PIMUNU
        q->pitt[s] = (PRECISION) (init_tmunu.pitt_in[sm] / hbarc);
        q->pitx[s] = (PRECISION) (init_tmunu.pitx_in[sm] / hbarc);
        q->pity[s] = (PRECISION) (init_tmunu.pity_in[sm] / hbarc);
        q->pitn[s] = (PRECISION) (init_tmunu.pitn_in[sm] / hbarc);
        q->pixx[s] = (PRECISION) (init_tmunu.pixx_in[sm] / hbarc);
        q->pixy[s] = (PRECISION) (init_tmunu.pixy_in[sm] / hbarc);
        q->pixn[s] = (PRECISION) (init_tmunu.pixn_in[sm] / hbarc);
        q->piyy[s] = (PRECISION) (init_tmunu.piyy_in[sm] / hbarc);
        q->piyn[s] = (PRECISION) (init_tmunu.piyn_in[sm] / hbarc);
        q->pinn[s] = (PRECISION) (init_tmunu.pinn_in[sm] / hbarc);
        #endif
        #ifdef PI
        q->Pi[s] = (PRECISION) (init_tmunu.Pi_in[sm] / hbarc);
        #endif
      }
    } // for (int j = 2; j < ny+2; ++j)
  } //for (int i = 2; i < nx+2; ++i)
}
