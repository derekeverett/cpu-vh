/*
 * DynamicalVariables.cpp
 *
 *  Created on: Oct 22, 2015
 *      Author: bazow
 */
#include <stdlib.h>

#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h" // for ghost cells

#include <omp.h>

//#define ROW_MAJ

CONSERVED_VARIABLES *q,*Q,*qS;

FLUID_VELOCITY *u,*up,*uS;

PRECISION *e, *p;

int columnMajorLinearIndex(int i, int j, int k, int nx, int ny, int nz) {
	#ifdef ROW_MAJ
	return (ny * nz) * i + (nz) * j + k;
	#else
	return i + nx * (j + ny * k);
	#endif
}

void allocateHostMemory(int len) {
	size_t bytes = sizeof(PRECISION);

	//=======================================================
	// Primary variables
	//=======================================================
	e = (PRECISION *)calloc(len, bytes);
	p = (PRECISION *)calloc(len,bytes);
	// fluid velocity at current time step
	u = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	u->ut = (PRECISION *)calloc(len,bytes);
	u->ux = (PRECISION *)calloc(len,bytes);
	u->uy = (PRECISION *)calloc(len,bytes);
	u->un = (PRECISION *)calloc(len,bytes);
	// fluid velocity at previous time step
	up = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	up->ut = (PRECISION *)calloc(len,bytes);
	up->ux = (PRECISION *)calloc(len,bytes);
	up->uy = (PRECISION *)calloc(len,bytes);
	up->un = (PRECISION *)calloc(len,bytes);
	// fluid velocity at intermediate time step
	uS = (FLUID_VELOCITY *)calloc(1, sizeof(FLUID_VELOCITY));
	uS->ut = (PRECISION *)calloc(len,bytes);
	uS->ux = (PRECISION *)calloc(len,bytes);
	uS->uy = (PRECISION *)calloc(len,bytes);
	uS->un = (PRECISION *)calloc(len,bytes);

	//=======================================================
	// Conserved variables
	//=======================================================
	q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	q->ttt = (PRECISION *)calloc(len, bytes);
	q->ttx = (PRECISION *)calloc(len, bytes);
	q->tty = (PRECISION *)calloc(len, bytes);
	q->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	q->pitt = (PRECISION *)calloc(len, bytes);
	q->pitx = (PRECISION *)calloc(len, bytes);
	q->pity = (PRECISION *)calloc(len, bytes);
	q->pitn = (PRECISION *)calloc(len, bytes);
	q->pixx = (PRECISION *)calloc(len, bytes);
	q->pixy = (PRECISION *)calloc(len, bytes);
	q->pixn = (PRECISION *)calloc(len, bytes);
	q->piyy = (PRECISION *)calloc(len, bytes);
	q->piyn = (PRECISION *)calloc(len, bytes);
	q->pinn = (PRECISION *)calloc(len, bytes);
#endif
	// allocate space for \Pi
#ifdef PI
	q->Pi = (PRECISION *)calloc(len, bytes);
#endif
	// upated variables at the n+1 time step
	Q = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	Q->ttt = (PRECISION *)calloc(len, bytes);
	Q->ttx = (PRECISION *)calloc(len, bytes);
	Q->tty = (PRECISION *)calloc(len, bytes);
	Q->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	Q->pitt = (PRECISION *)calloc(len, bytes);
	Q->pitx = (PRECISION *)calloc(len, bytes);
	Q->pity = (PRECISION *)calloc(len, bytes);
	Q->pitn = (PRECISION *)calloc(len, bytes);
	Q->pixx = (PRECISION *)calloc(len, bytes);
	Q->pixy = (PRECISION *)calloc(len, bytes);
	Q->pixn = (PRECISION *)calloc(len, bytes);
	Q->piyy = (PRECISION *)calloc(len, bytes);
	Q->piyn = (PRECISION *)calloc(len, bytes);
	Q->pinn = (PRECISION *)calloc(len, bytes);
#endif
	// allocate space for \Pi
#ifdef PI
	Q->Pi = (PRECISION *)calloc(len, bytes);
#endif
	// updated variables at the intermediate time step
	qS = (CONSERVED_VARIABLES *)calloc(1, sizeof(CONSERVED_VARIABLES));
	qS->ttt = (PRECISION *)calloc(len, bytes);
	qS->ttx = (PRECISION *)calloc(len, bytes);
	qS->tty = (PRECISION *)calloc(len, bytes);
	qS->ttn = (PRECISION *)calloc(len, bytes);
#ifdef PIMUNU
	qS->pitt = (PRECISION *)calloc(len, bytes);
	qS->pitx = (PRECISION *)calloc(len, bytes);
	qS->pity = (PRECISION *)calloc(len, bytes);
	qS->pitn = (PRECISION *)calloc(len, bytes);
	qS->pixx = (PRECISION *)calloc(len, bytes);
	qS->pixy = (PRECISION *)calloc(len, bytes);
	qS->pixn = (PRECISION *)calloc(len, bytes);
	qS->piyy = (PRECISION *)calloc(len, bytes);
	qS->piyn = (PRECISION *)calloc(len, bytes);
	qS->pinn = (PRECISION *)calloc(len, bytes);
#endif
	// allocate space for \Pi
#ifdef PI
	qS->Pi = (PRECISION *)calloc(len, bytes);
#endif
}

void setConservedVariables(double t, void * latticeParams) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	#pragma omp parallel for collapse(3)
	#ifdef TILE
	#pragma unroll_and_jam
	#endif
	//#endif
	for (int k = N_GHOST_CELLS_M; k < nz+N_GHOST_CELLS_M; ++k) {
		for (int j = N_GHOST_CELLS_M; j < ny+N_GHOST_CELLS_M; ++j) {
			for (int i = N_GHOST_CELLS_M; i < nx+N_GHOST_CELLS_M; ++i) {
				int s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);

				PRECISION ux_s = u->ux[s];
				PRECISION uy_s = u->uy[s];
				PRECISION un_s = u->un[s];
				PRECISION ut_s = u->ut[s];
				PRECISION e_s = e[s];
				PRECISION p_s = p[s];

				PRECISION pitt_s = 0;
				PRECISION pitx_s = 0;
				PRECISION pity_s = 0;
				PRECISION pitn_s = 0;
#ifdef PIMUNU
				pitt_s = q->pitt[s];
				pitx_s = q->pitx[s];
				pity_s = q->pity[s];
				pitn_s = q->pitn[s];
#endif
				PRECISION Pi_s = 0;
#ifdef PI
				Pi_s = q->Pi[s];
#endif

				q->ttt[s] = Ttt(e_s, p_s+Pi_s, ut_s, pitt_s);
				q->ttx[s] = Ttx(e_s, p_s+Pi_s, ut_s, ux_s, pitx_s);
				q->tty[s] = Tty(e_s, p_s+Pi_s, ut_s, uy_s, pity_s);
				q->ttn[s] = Ttn(e_s, p_s+Pi_s, ut_s, un_s, pitn_s);
			}
		}
	}
}

void setGhostCells(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams
) {
	setGhostCellsKernelI(q,e,p,u,latticeParams);
	setGhostCellsKernelJ(q,e,p,u,latticeParams);
	setGhostCellsKernelK(q,e,p,u,latticeParams);
}

void setGhostCellVars(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u,
int s, int sBC) {
	e[s] = e[sBC];
	p[s] = p[sBC];
	u->ut[s] = u->ut[sBC];
	u->ux[s] = u->ux[sBC];
	u->uy[s] = u->uy[sBC];
	u->un[s] = u->un[sBC];
	q->ttt[s] = q->ttt[sBC];
	q->ttx[s] = q->ttx[sBC];
	q->tty[s] = q->tty[sBC];
	q->ttn[s] = q->ttn[sBC];
	// set \pi^\mu\nu ghost cells if evolved
#ifdef PIMUNU
	q->pitt[s] = q->pitt[sBC];
	q->pitx[s] = q->pitx[sBC];
	q->pity[s] = q->pity[sBC];
	q->pitn[s] = q->pitn[sBC];
	q->pixx[s] = q->pixx[sBC];
	q->pixy[s] = q->pixy[sBC];
	q->pixn[s] = q->pixn[sBC];
	q->piyy[s] = q->piyy[sBC];
	q->piyn[s] = q->piyn[sBC];
	q->pinn[s] = q->pinn[sBC];
#endif
	// set \Pi ghost cells if evolved
#ifdef PI
	q->Pi[s] = q->Pi[sBC];
#endif
}

void setGhostCellsKernelI(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams
) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nx,ncx,ncy,ncz;
	nx = lattice->numLatticePointsX;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	ncz = lattice->numComputationalLatticePointsRapidity;

	int iBC,s,sBC;
	#pragma omp parallel for
	for(int j = 2; j < ncy; ++j) {
		for(int k = 2; k < ncz; ++k) {
			iBC = 2;
			//#pragma omp parallel for
			for (int i = 0; i <= 1; ++i) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
			iBC = nx + 1;
			//#pragma omp parallel for
			for (int i = nx + 2; i <= nx + 3; ++i) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(iBC, j, k, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
		}
	}
}

void setGhostCellsKernelJ(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams
) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int ny,ncx,ncy,ncz;
	ny = lattice->numLatticePointsY;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	ncz = lattice->numComputationalLatticePointsRapidity;

	int jBC,s,sBC;
	#pragma omp parallel for
	for(int i = 2; i < ncx; ++i) {
		for(int k = 2; k < ncz; ++k) {
			jBC = 2;
			//#pragma omp parallel for
			for (int j = 0; j <= 1; ++j) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
			jBC = ny + 1;
			//#pragma omp parallel for
			for (int j = ny + 2; j <= ny + 3; ++j) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(i, jBC, k, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
		}
	}
}

void setGhostCellsKernelK(CONSERVED_VARIABLES * const __restrict__ q,
PRECISION * const __restrict__ e, PRECISION * const __restrict__ p,
FLUID_VELOCITY * const __restrict__ u, void * latticeParams
) {
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;

	int nz,ncx,ncy;
	nz = lattice->numLatticePointsRapidity;
	ncx = lattice->numComputationalLatticePointsX;
	ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;

	int kBC,s,sBC;
	#pragma omp parallel for
	for(int i = 2; i < ncx; ++i) {
		for(int j = 2; j < ncy; ++j) {
			kBC = 2;
			//#pragma omp parallel for
			for (int k = 0; k <= 1; ++k) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
			kBC = nz + 1;
			//#pragma omp parallel for
			for (int k = nz + 2; k <= nz + 3; ++k) {
				s = columnMajorLinearIndex(i, j, k, ncx, ncy, ncz);
				sBC = columnMajorLinearIndex(i, j, kBC, ncx, ncy, ncz);
				setGhostCellVars(q,e,p,u,s,sBC);
			}
		}
	}
}

void swap(CONSERVED_VARIABLES **arr1, CONSERVED_VARIABLES **arr2) {
	CONSERVED_VARIABLES *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}

void setCurrentConservedVariables() {
	swap(&q, &Q);
}

void swapFluidVelocity(FLUID_VELOCITY **arr1, FLUID_VELOCITY **arr2) {
	FLUID_VELOCITY *tmp = *arr1;
	*arr1 = *arr2;
	*arr2 = tmp;
}

void freeHostMemory() {
	free(e);
	free(p);
	free(u->ut);
	free(u->ux);
	free(u->uy);
	free(u->un);

	free(q->ttt);
	free(q->ttx);
	free(q->tty);
	free(q->ttn);
	// free \pi^\mu\nu
#ifdef PIMUNU
	free(q->pitt);
	free(q->pitx);
	free(q->pity);
	free(q->pitn);
	free(q->pixx);
	free(q->pixy);
	free(q->pixn);
	free(q->piyy);
	free(q->piyn);
	free(q->pinn);
#endif
	// free \Pi
#ifdef PI
	free(q->Pi);
#endif
	free(q);
}
