/*
* HydroPlugin.c
*
*  Created on: Oct 23, 2015
*      Author: bazow
*/

#include <stdlib.h>
#include <stdio.h> // for printf

// for timing
#include <ctime>
#include <iostream>

//for cornelius and writing freezeout file
#include <fstream>
#include "cornelius.h"
#include "FreezeOut.h"
#include "Memory.h"

#include "edu/osu/rhic/harness/hydro/HydroPlugin.h"
#include "edu/osu/rhic/trunk/hydro/DynamicalVariables.h"
#include "edu/osu/rhic/harness/lattice/LatticeParameters.h"
#include "edu/osu/rhic/harness/ic/InitialConditionParameters.h"
#include "edu/osu/rhic/harness/hydro/HydroParameters.h"
#include "edu/osu/rhic/harness/io/FileIO.h"
#include "edu/osu/rhic/trunk/ic/InitialConditions.h"
#include "edu/osu/rhic/trunk/hydro/FullyDiscreteKurganovTadmorScheme.h"
#include "edu/osu/rhic/trunk/hydro/EnergyMomentumTensor.h"
#include "edu/osu/rhic/trunk/eos/EquationOfState.h"

#define FREQ 10 //write output to file every FREQ timesteps
#define FOFREQ 5 //call freezeout surface finder every FOFREQ timesteps
#define FOTEST 1 //if true, freezeout surface file is written with proper times rounded (down) to step size

void outputDynamicalQuantities(double t, const char *outputDir, void * latticeParams)
{
  output(e, t, outputDir, "e", latticeParams);
  output(u->ux, t, outputDir, "ux", latticeParams);
  output(u->uy, t, outputDir, "uy", latticeParams);
  //	output(u->un, t, outputDir, "un", latticeParams);
  output(u->ut, t, outputDir, "ut", latticeParams);
  //	output(q->ttt, t, outputDir, "ttt", latticeParams);
  //	output(q->ttn, t, outputDir, "ttn", latticeParams);
  #ifdef PIMUNU
  //output(q->pixx, t, outputDir, "pixx", latticeParams);
  //output(q->pixy, t, outputDir, "pixy", latticeParams);
  //output(q->pixn, t, outputDir, "pixn", latticeParams);
  //output(q->piyy, t, outputDir, "piyy", latticeParams);
  //output(q->piyn, t, outputDir, "piyn", latticeParams);
  //output(q->pinn, t, outputDir, "pinn", latticeParams);
  #endif
  #ifdef PI
  //output(q->Pi, t, outputDir, "Pi", latticeParams);
  #endif

}

void run(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir)
{
  struct LatticeParameters * lattice = (struct LatticeParameters *) latticeParams;
  struct InitialConditionParameters * initCond = (struct InitialConditionParameters *) initCondParams;
  struct HydroParameters * hydro = (struct HydroParameters *) hydroParams;

  /************************************************************************************	\
  * System configuration
  /************************************************************************************/
  int nt = lattice->numProperTimePoints;
  int nx = lattice->numLatticePointsX;
  int ny = lattice->numLatticePointsY;
  int nz = lattice->numLatticePointsRapidity;

  if ( (nx < 2) || (ny < 2) || (nz < 2) ) printf("!simulation is not in 3+1d; freezeout finder will not work - fix this!\n");
  int ncx = lattice->numComputationalLatticePointsX;
  int ncy = lattice->numComputationalLatticePointsY;
  int ncz = lattice->numComputationalLatticePointsRapidity;
  int nElements = ncx * ncy * ncz;

  double t0 = hydro->initialProperTimePoint;
  double dt = lattice->latticeSpacingProperTime;
  double dx = lattice->latticeSpacingX;
  double dy = lattice->latticeSpacingY;
  double dz = lattice->latticeSpacingRapidity;
  double e0 = initCond->initialEnergyDensity;

  double freezeoutTemperatureGeV = hydro->freezeoutTemperatureGeV;
  const double hbarc = 0.197326938;
  const double freezeoutTemperature = freezeoutTemperatureGeV/hbarc;
  //const double freezeoutEnergyDensity = e0*pow(freezeoutTemperature,4);
  const double freezeoutEnergyDensity = equilibriumEnergyDensity(freezeoutTemperature);
  printf("Grid size = %d x %d x %d\n", nx, ny, nz);
  printf("spatial resolution = (%.3f, %.3f, %.3f)\n", lattice->latticeSpacingX, lattice->latticeSpacingY, lattice->latticeSpacingRapidity);
  printf("freezeout temperature = %.3f [fm^-1] (eF = %.3f [fm^-4])\n", freezeoutTemperature, freezeoutEnergyDensity);

  // allocate memory
  allocateHostMemory(nElements);

  //initialize cornelius for freezeout surface finding
  //see example_4d() in example_cornelius
  //this works only for full 3+1 d simulation? need to find a way to generalize to n+1 d
  int dim;
  double *lattice_spacing;
  if ((nx > 1) && (ny > 1) && (nz > 1))
  {
    dim = 4;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
    lattice_spacing[3] = dz;
  }
  else if ((nx > 1) && (ny > 1) && (nz < 2))
  {
    dim = 3;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
  }
  else if ((nx > 1) && (ny < 2) && (nz < 2))
  {
    dim = 2;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt;
    lattice_spacing[1] = dx;
  }

  Cornelius cor;
  cor.init(dim, freezeoutEnergyDensity, lattice_spacing);

  double ****energy_density_evoution;
  energy_density_evoution = calloc4dArray(energy_density_evoution, FOFREQ+1, nx, ny, nz);

  //make an array to store all the hydrodynamic variables for FOFREQ time steps
  //to be written to file once the freezeout surface is determined by the critical energy density
  int n_hydro_vars = 16; //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, the temperature and pressure are calclated with EoS
  double *****hydrodynamic_evoution;
  hydrodynamic_evoution = calloc5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny, nz);

  //find a way to make this more compact/readable, and more readily parallelizable
  double ****hyperCube;
  hyperCube = calloc4dArray(hyperCube, 2, 2, 2, 2);
  //open the freezeout surface file
  ofstream freezeoutSurfaceFile;
  freezeoutSurfaceFile.open("output/freezeoutSurface.dat"); //find a way to do this using command line args

  /************************************************************************************	\
  * Fluid dynamic initialization
  /************************************************************************************/
  double t = t0;
  // generate initial conditions
  setInitialConditions(latticeParams, initCondParams, hydroParams, rootDirectory);
  // Calculate conserved quantities
  setConservedVariables(t, latticeParams);
  // impose boundary conditions with ghost cells
  setGhostCells(q,e,p,u,latticeParams);

  /************************************************************************************	\
  * Evolve the system in time
  /************************************************************************************/
  int ictr = (nx % 2 == 0) ? ncx/2 : (ncx-1)/2;
  int jctr = (ny % 2 == 0) ? ncy/2 : (ncy-1)/2;
  int kctr = (nz % 2 == 0) ? ncz/2 : (ncz-1)/2;
  int sctr = columnMajorLinearIndex(ictr, jctr, kctr, ncx, ncy);

  std::clock_t t1,t2;

  double totalTime = 0;
  int nsteps = 0;

  int accumulator1 = 0;
  int accumulator2 = 0;
  // evolve in time
  for (int n = 1; n <= nt+1; ++n)
  {
    // copy variables back to host and write to disk
    if ((n-1) % FREQ == 0) {
      printf("n = %d:%d (t = %.3f),\t (e, p) = (%.3f, %.3f) [fm^-4],\t (T = %.3f [GeV]),\t",
      n - 1, nt, t, e[sctr], p[sctr], effectiveTemperature(e[sctr])*hbarc);
      outputDynamicalQuantities(t, outputDir, latticeParams);
      // end hydrodynamic simulation if the temperature is below the freezeout temperature
      //if(e[sctr] < freezeoutEnergyDensity) {
      //printf("\nReached freezeout temperature at the center.\n");
      //break;
      //}
    }

    //append the energy density and all hydro variables to storage arrays
    int nFO = n % FOFREQ;

    if(nFO == 0) //swap in the old values so that freezeout volume elements have overlap between calls to finder
    {
      swapAndSetHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ);
    }
    else //update the values of the rest of the array with current time step
    {
      setHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ, n);
    }

    //the n=1 values are written to the it = 2 index of array, so don't start until here
    int start;
    if (n <= FOFREQ) start = 2;
    else start = 0;
    if (nFO == FOFREQ - 1) //call the freezeout finder should this be put before the values are set?
    {

      //besides writing centroid and normal to file, write all the hydro variables
      for (int it = start; it < FOFREQ; it++) //note* avoiding boundary problems (reading outside array)
      {
        for (int ix = 0; ix < nx-1; ix++)
        {
          for (int iy = 0; iy < ny-1; iy++)
          {
            for (int iz = 0; iz < nz-1; iz++)
            {
              //write the values of energy density to all corners of the hyperCube
              writeEnergyDensityToHypercube(hyperCube, energy_density_evoution, it, ix, iy, iz);

              //the freezeout surface file is written in the same format that is
              //written by MUSIC hydro code (see readFreezeOutSurface() in freeze.cpp)

              //use cornelius to find the centroid and normal vector of each hyperCube
              cor.find_surface_4d(hyperCube);
              //write centroid and normal of each surface element to file
              for (int i = 0; i < cor.get_Nelements(); i++)
              {
                double temp = 0.0; //temporary variable
                //first write the position of the centroid of surface element
                double cell_tau = t0 + ((double)(n - FOFREQ + it)) * dt; //check if this is the correct time!
                double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);
                double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
                double cell_z = (double)iz * dz  - (((double)(nz-1)) / 2.0 * dz);

                if (FOTEST) {freezeoutSurfaceFile << cell_tau << " ";}
                else {freezeoutSurfaceFile << cor.get_centroid_elem(i,0) + cell_tau << " ";}
                freezeoutSurfaceFile << cor.get_centroid_elem(i,1) + cell_x << " ";
                freezeoutSurfaceFile << cor.get_centroid_elem(i,2) + cell_y << " ";
                freezeoutSurfaceFile << cor.get_centroid_elem(i,3) + cell_z << " ";
                //then the surface normal element; note jacobian factors of +/- tau for milne coordinates
                freezeoutSurfaceFile << t * cor.get_normal_elem(i,0) << " ";
                freezeoutSurfaceFile << (-1.0) * t * cor.get_normal_elem(i,1) << " ";
                freezeoutSurfaceFile << (-1.0) * t * cor.get_normal_elem(i,2) << " ";
                freezeoutSurfaceFile << (-1.0) * t * cor.get_normal_elem(i,3) << " ";
                //write all the necessary hydro dynamic variables by first performing linear interpolation from values at
                //corners of hypercube
                double tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
                double x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
                double y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];
                double z_frac = cor.get_centroid_elem(i,3) / lattice_spacing[3];

                //first write the flow velocity
                for (int ivar = 0; ivar < 4; ivar++)
                {
                  temp = interpolateVariable(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                  freezeoutSurfaceFile << temp << " ";
                }

                //write the energy density
                temp = interpolateVariable(hydrodynamic_evoution, 4, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                freezeoutSurfaceFile << hbarc * temp << " "; //note factors of hbarc to give units (GeV/fm^3)
                //the temperature !this needs to be checked
                freezeoutSurfaceFile << hbarc * effectiveTemperature(temp) << " ";
                //the baryon chemical potential, writing zero for now
                freezeoutSurfaceFile << 0.0 << " ";
                //  (e + P) / T , the entropy density for zero chem. potentials !check this, note we could be a divide by zero problem if T=0!
                double e_plus_P_over_T = (temp + equilibriumPressure(temp)) / effectiveTemperature(temp);
                freezeoutSurfaceFile << e_plus_P_over_T << " ";
                //write ten components of pi_(mu,nu) shear viscous tensor
                for (int ivar = 5; ivar < 15; ivar++)
                {
                  temp = interpolateVariable(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                  freezeoutSurfaceFile << hbarc * temp << " ";
                }
                //write the bulk pressure Pi, and start a new line
                temp = interpolateVariable(hydrodynamic_evoution, 15, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                freezeoutSurfaceFile << hbarc * temp << endl;
              }
            }
          }
        }
      }
    }

    //if all cells are below freezeout temperature end hydro
    accumulator1 = 0;
    for (int ix = 2; ix < nx+2; ix++)
    {
      for (int iy = 2; iy < ny+2; iy++)
      {
        for (int iz = 2; iz < nz+2; iz++)
        {
          int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4);
          if (e[s] > freezeoutEnergyDensity) accumulator1 = accumulator1 + 1;
        }
      }
    }
    if (accumulator1 == 0) accumulator2 += 1;
    if (accumulator2 >= FOFREQ+1) //only break once freezeout finder has had a chance to search/write to file
    {
      printf("\nAll cells have dropped below freezeout energy density\n");
      break;
    }

    t1 = std::clock();
    rungeKutta2(t, dt, q, Q, latticeParams, hydroParams);
    t2 = std::clock();
    double delta_time = (t2 - t1) / (double)(CLOCKS_PER_SEC / 1000);
    if ((n-1) % FREQ == 0) printf("(Elapsed time: %.3f ms)\n",delta_time);
    totalTime+=delta_time;
    ++nsteps;

    setCurrentConservedVariables();

    t = t0 + n * dt;
  }
  printf("Average time/step: %.3f ms\n",totalTime/((double)nsteps));

  freezeoutSurfaceFile.close();
  /************************************************************************************	\
  * Deallocate host memory
  /************************************************************************************/
  freeHostMemory();

  //Deallocate memory used for freezeout finding
  free4dArray(energy_density_evoution, FOFREQ+1, nx, ny);
  free5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny);
  delete [] lattice_spacing;

  free4dArray(hyperCube, 2, 2, 2);
}
