/*
 * DynamicalSources.cpp
 *
 *  Created on: Nov 25, 2017
 *  Author: Lipei
 */

#include <math.h> // for math functions
#include <stdio.h> // for printf
#include <stdlib.h> //TEMP
#include <iostream>
#include <istream>
#include <fstream>
#include <cassert>
#include <string>

#include <iomanip>//by Lipei

#include "../include/DynamicalVariables.h"
#include "../include/DynamicalSources.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"

#include <H5Cpp.h>
#include <H5File.h>

using namespace std;
//*********************************************************************************************************\
//* Initialize the dynamical source terms
//*********************************************************************************************************/

void readInSource(int n, void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory)
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
	struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	// For each time step, the total number of cells in 3D space
	int nElements = nx * ny * nz;

	double t0 = hydro->initialProperTimePoint;
	double dt = lattice->latticeSpacingProperTime;

	double t = t0 + (n-1)* dt;
	double time;

	//FILE *sourcefile;
	char fname[255];
	sprintf(fname, "%s/%s%d.h5", rootDirectory, "input/DynamicalSources/Sources",n);
	//sourcefile = fopen(fname, "r");

	H5::H5File file( fname, H5F_ACC_RDONLY );
	H5::DataSet dataset = file.openDataSet( "data" );
	H5::DataSpace dataspace = dataset.getSpace();
  int rank = dataspace.getSimpleExtentNdims();
  hsize_t dims_out[4];
  int ndims = dataspace.getSimpleExtentDims( dims_out, NULL);
	//H5::DataSpace mspace(RANK, dims);

	//read the source terms (Sb, St, Sx, Sy, Sn) into a single array
	float *Sall;
	Sall = (float *)calloc( 5*nElements, sizeof(float) );
	dataset.read( Sall, H5::PredType::NATIVE_FLOAT, dataspace);

	//split this array into corresponding source terms
	for(int i = 2; i < nx+2; ++i){
		for(int j = 2; j < ny+2; ++j){
			for(int k = 2; k < nz+2; ++k){
				//this index runs over grid with ghost cells
				int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
				//this index runs over grid without ghost cells
				int s_m = columnMajorLinearIndex(i-2, j-2, k-2, nx, ny, nz);
				Source->sourceb[s]=(PRECISION)Sall[s_m];
				Source->sourcet[s]=(PRECISION)Sall[nElements + s_m];
				Source->sourcex[s]=(PRECISION)Sall[2*nElements + s_m];
				Source->sourcey[s]=(PRECISION)Sall[3*nElements + s_m];
				if (nz > 1) Source->sourcen[s]=(PRECISION)Sall[4*nElements + s_m];
			} //for(int k = 2; k < nz+2; ++k)
		} // for(int j = 2; j < ny+2; ++j)
	} // for(int i = 2; i < nx+2; ++i)

	/*
  cout << "rank " << rank << ", dimensions " << (unsigned long)(dims_out[0])
	<< " x " << (unsigned long)(dims_out[1]) << " x " << (unsigned long)(dims_out[2])
	<< " x " << (unsigned long)(dims_out[3]) << endl;
	*/

	/*
	if(sourcefile==NULL){
		printf("The source files could not be opened...\n");
		exit(-1);
	}
	else
	{
		fseek(sourcefile,0L,SEEK_SET);

		//for(int i=0; i<(nElements+1)*(n-1); i++) fscanf(sourcefile,"%*[^\n]%*c");//Skip the title line and all the cells read in by previous steps, (nElements+1) lines

		//fscanf(sourcefile,"%*s%le%*c", &time);
		//printf("time=%lf\n",time);
		//if(time-t>1.e-20) printf("The dynamical source at a wrong time step is being read in. tSource=%lf, tCode=%lf\n", time, t);
		//if(time==t) printf("The dynamical source starts to be read in at %lf.\n", time);

		for(int i = 2; i < nx+2; ++i){
			for(int j = 2; j < ny+2; ++j){
				for(int k = 2; k < nz+2; ++k){
					int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
					fscanf(sourcefile,"%*s%*s%*s %le %le %le %le %le", & Source->sourcet[s], & Source->sourcex[s], & Source->sourcey[s], & Source->sourcen[s], & Source->sourceb[s]);
					//printf("%le\t %le\t %le\t %le\t %le\n", Source->sourcet[s], Source->sourcex[s], Source->sourcey[s], Source->sourcen[s], Source->sourceb[s]);
					//Source->sourcet[s]=10*Source->sourcet[s];
					//Source->sourcex[s]=10*Source->sourcex[s];
					//Source->sourcex[s]=0;
					//Source->sourcey[s]=0;
					//Source->sourceb[s]=10*Source->sourceb[s];
				}
			}
		}
	}
	fclose(sourcefile);
	*/
}

void noSource(void * latticeParams, void * initCondParams)
{
    struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
    struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

    int nx = lattice->numLatticePointsX;
    int ny = lattice->numLatticePointsY;
    int nz = lattice->numLatticePointsRapidity;

    for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                Source->sourcet[s] = 0.0;
                Source->sourcex[s] = 0.0;
                Source->sourcey[s] = 0.0;
                Source->sourcen[s] = 0.0;
                Source->sourceb[s] = 0.0;
            }//k
        }//j
    }//i
}

//*********************************************************************************************************\
//* Dynamical source terms from the jet traversing the medium
//*********************************************************************************************************/

void setDynamicalSources(void * latticeParams, void * initCondParams, double *dp_dtau, double *pos) //dp_dtau is the jet energy loss, pos is the jet position
{
	struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
	struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;

	int nx = lattice->numLatticePointsX;
	int ny = lattice->numLatticePointsY;
	int nz = lattice->numLatticePointsRapidity;
	int ncx = lattice->numComputationalLatticePointsX;
	int ncy = lattice->numComputationalLatticePointsY;
	int ncz = lattice->numComputationalLatticePointsRapidity;
	double dx = lattice->latticeSpacingX;
	double dy = lattice->latticeSpacingY;
	double dz = lattice->latticeSpacingRapidity;

	PRECISION dummy = 0.0;

    double xmin = (-1.0) * ((double)(nx - 1) / 2.0) * dx;
	double ymin = (-1.0) * ((double)(ny - 1) / 2.0) * dy;
	double zmin = (-1.0) * ((double)(nz - 1) / 2.0) * dz;

	//construct an array of the gaussian smeared jet position
	double smearedPosition[ncx * ncy * ncz];
	double width = 0.1; //width of gaussian smearing

	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                double x = (double)i * dx + xmin;
                double y = (double)j * dy + ymin;
                double z = (double)k * dz + zmin;
                smearedPosition[s] = exp((-1.0)*(pos[1] - x) * (pos[1] - x) / width) * exp((-1.0)*(pos[2] - y) * (pos[2] - y) / width) * exp((-1.0)*(pos[3] - z) * (pos[3] - z) / width);
            }//k
        }//j
    }//i

    //now multiply the smeared position by energy loss corresponding to vector components
	for(int i = 2; i < nx+2; ++i){
        for(int j = 2; j < ny+2; ++j){
            for(int k = 2; k < nz+2; ++k){
                int s = columnMajorLinearIndex(i, j, k, nx+4, ny+4, nz+4);
                Source->sourcet[s] = -dp_dtau[0] * smearedPosition[s];
                Source->sourcex[s] = -dp_dtau[1] * smearedPosition[s];
                Source->sourcey[s] = -dp_dtau[2] * smearedPosition[s];
                Source->sourcen[s] = -dp_dtau[3] * smearedPosition[s];
                Source->sourceb[s] = dummy;
            }//k
        }//j
    }//i
}
