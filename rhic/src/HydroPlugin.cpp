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

//for FO surface
#include <vector>

//for cornelius and writing freezeout file
#include <fstream>
//#include "cornelius-c++-1.3/cornelius.h"
#include "FreezeOut.cpp"
#include "Memory.cpp"

#include "../include/HydroPlugin.h"
#include "../include/DynamicalVariables.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroInitialTmunu.h"
#include "../include/HydroParameters.h"
#include "../include/FileIO.h"
#include "../include/InitialConditions.h"
#include "../include/FullyDiscreteKurganovTadmorScheme.h"
#include "../include/EnergyMomentumTensor.h"
#include "../include/EquationOfState.h"
#include "../include/DynamicalSources.h"
#include "../include/FreezeOut.h"

#include "../include/Vorticity.h" //for polarization studies

#define FREQ 10 //write output to file every FREQ timesteps
#define FOFREQ 10 //call freezeout surface finder every FOFREQ timesteps
#define FOTEST 0 //if true, freezeout surface file is written with proper times rounded (down) to step size
#define FOFORMAT 0 // 0 : write f.o. surface to ASCII file ;  1 : write to binary file

void outputDynamicalQuantities(double t, const char *outputDir, void * latticeParams)
{
  output(e, t, outputDir, "e", latticeParams);
  output(u->ux, t, outputDir, "ux", latticeParams);
  output(u->uy, t, outputDir, "uy", latticeParams);
  output(u->un, t, outputDir, "un", latticeParams);
  output(u->ut, t, outputDir, "ut", latticeParams);
  output(q->ttt, t, outputDir, "ttt", latticeParams);
  output(q->ttn, t, outputDir, "ttn", latticeParams);
  #ifdef PIMUNU
  output(q->pixx, t, outputDir, "pixx", latticeParams);
  output(q->pixy, t, outputDir, "pixy", latticeParams);
  output(q->pixn, t, outputDir, "pixn", latticeParams);
  //output(q->piyy, t, outputDir, "piyy", latticeParams);
  //output(q->piyn, t, outputDir, "piyn", latticeParams);
  output(q->pinn, t, outputDir, "pinn", latticeParams);
  output(regFacShear, t, outputDir, "regFacShear", latticeParams);
  #endif
  #ifdef PI
  output(q->Pi, t, outputDir, "Pi", latticeParams);
  output(regFacBulk, t, outputDir, "regFacBulk", latticeParams);
  #endif
  #ifdef THERMAL_VORTICITY
  output(wmunu->wtx, t, outputDir, "wtx", latticeParams);
  output(wmunu->wty, t, outputDir, "wty", latticeParams);
  output(wmunu->wtn, t, outputDir, "wtn", latticeParams);
  output(wmunu->wxy, t, outputDir, "wxy", latticeParams);
  output(wmunu->wxn, t, outputDir, "wxn", latticeParams);
  output(wmunu->wyn, t, outputDir, "wyn", latticeParams);
  #endif
}

void run(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir, HydroInitialTmunu init_tmunu)
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
  int coarseProperTimeFOFactor = lattice->coarseProperTimeFOFactor;

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

  int initialConditionType = initCond->initialConditionType;
  int numberOfSourceFiles = initCond->numberOfSourceFiles;

  double freezeoutTemperatureGeV = hydro->freezeoutTemperatureGeV;
  const double hbarc = 0.197326938;
  const double freezeoutTemperature = freezeoutTemperatureGeV/hbarc;
  const double freezeoutEnergyDensity = equilibriumEnergyDensity(freezeoutTemperature);
  printf("Grid size = %d x %d x %d\n", nx, ny, nz);
  printf("spatial resolution = (%.3f, %.3f, %.3f)\n", lattice->latticeSpacingX, lattice->latticeSpacingY, lattice->latticeSpacingRapidity);
  printf("freezeout temperature = %.3f [fm^-1] (eF = %.3f [fm^-4])\n", freezeoutTemperature, freezeoutEnergyDensity);
  printf("Coarse graining FO surface by factor %d in proper time\n", coarseProperTimeFOFactor);
  #ifndef IDEAL
  printf("Regulating viscous currents using : ");
  #ifdef REG_SCHEME_1
  printf("Regulation Scheme 1\n");
  #endif
  #ifdef REG_SCHEME_2
  printf("Regulation Scheme 2\n");
  #endif
  #ifdef REG_SCHEME_3
  printf("Regulation Scheme 3\n");
  #endif
  #endif
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
    lattice_spacing[0] = dt * coarseProperTimeFOFactor;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
    lattice_spacing[3] = dz;
  }
  else if ((nx > 1) && (ny > 1) && (nz < 2))
  {
    dim = 3;
    lattice_spacing = new double[dim];
    lattice_spacing[0] = dt * coarseProperTimeFOFactor;
    lattice_spacing[1] = dx;
    lattice_spacing[2] = dy;
  }
  else printf("simulation is not in 3+1D or 2+1D; freezeout finder will not work!\n");


  //Cornelius cor;
  //cor.init(dim, freezeoutEnergyDensity, lattice_spacing);

  //declare a vector of FO_Element to hold FO cell info
  std::vector<FO_Element> fo_surf;

  double ****energy_density_evoution;
  energy_density_evoution = calloc4dArray(energy_density_evoution, FOFREQ+1, nx, ny, nz);

  //make an array to store all the hydrodynamic variables for FOFREQ time steps
  //to be written to file once the freezeout surface is determined by the critical energy density
  #ifdef THERMAL_VORTICITY
  int n_hydro_vars = 16; //u1, u2, u3, e, pi11, pi12, pi13, pi22, pi23, Pi, wtx, wty, wtn, wxy, wxn, wyn (the temperature and pressure are calclated with EoS)
  #else
  int n_hydro_vars = 10; //u1, u2, u3, e, pi11, pi12, pi13, pi22, pi23, Pi (the temperature and pressure are calclated with EoS)
  #endif

  //rather than declaring a multidimensional array, declare a vector of FO_elements to which we save the hydro info?
  //but we need to interpolate between cells, and this is easily done with a multidimensional array ...?
  double *****hydrodynamic_evoution;
  hydrodynamic_evoution = calloc5dArray(hydrodynamic_evoution, n_hydro_vars, FOFREQ+1, nx, ny, nz);

  //for 3+1D simulations
  double ****hyperCube4D;
  hyperCube4D = calloc4dArray(hyperCube4D, 2, 2, 2, 2);
  //for 2+1D simulations
  double ***hyperCube3D;
  hyperCube3D = calloc3dArray(hyperCube3D, 2, 2, 2);

  //open the freezeout surface file
  ofstream freezeoutSurfaceFile;
  std::string surf_fname = std::string(outputDir) + std::string("/surface.dat");
  freezeoutSurfaceFile.open(surf_fname);

  /************************************************************************************	\
  * Fluid dynamic initialization
  /************************************************************************************/
  double t = t0;
  // generate initial conditions from internal routines or reading from files
  if (initialConditionType < 13) setInitialConditions(latticeParams, initCondParams, hydroParams, rootDirectory);
  //set initial conditions from initial energy density vector
  else if (initialConditionType == 13) setICFromEnergyDensityVector(latticeParams, initCondParams, hydroParams, rootDirectory, init_tmunu);
  //set initial conditions from preequilibrium vectors
  else if (initialConditionType == 14) setICFromPreequilVectors(latticeParams, initCondParams, hydroParams, rootDirectory, init_tmunu);
  else {printf("NOT A VALID INITIAL CONDITION! \n"); exit(-1);}

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
  int sctr = columnMajorLinearIndex(ictr, jctr, kctr, ncx, ncy, ncz);

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
    }

    //************************************************************************************\
    // Freeze-out finder (Derek)
    // the freezeout surface file is written in the format which can
    // be read by iS3D : https://github.com/derekeverett/iS3D
    //************************************************************************************/

    //append the energy density and all hydro variables to storage arrays
    int nFO = n % FOFREQ;

    //for vorticity and polzn studies
    #ifdef THERMAL_VORTICITY
    //swap in the old values so that freezeout volume elements have overlap between calls to finder
    if (nFO == 0) swapAndSetHydroVariables_Vorticity(energy_density_evoution, hydrodynamic_evoution, q, e, u, wmunu, nx, ny, nz, FOFREQ);
    //update the values of the rest of the array with current time step
    else setHydroVariables_Vorticity(energy_density_evoution, hydrodynamic_evoution, q, e, u, wmunu, nx, ny, nz, FOFREQ, n);
    #else
    //swap in the old values so that freezeout volume elements have overlap between calls to finder
    if (nFO == 0) swapAndSetHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ);
    //update the values of the rest of the array with current time step
    else setHydroVariables(energy_density_evoution, hydrodynamic_evoution, q, e, u, nx, ny, nz, FOFREQ, n);
    #endif

    //the n=1 values are written to the it = 2 index of array, so don't start until here
    int start;
    if (n <= FOFREQ) start = 2;
    //if (n <= FOFREQ) start = 1;
    else start = 0;
    //call the freezeout finder
    if (nFO == FOFREQ - 1) callFOFinder(dim, start, nx, ny, nz, n, t0, dt, t, dx, dy, dz, lattice_spacing, freezeoutEnergyDensity,
                                        hyperCube4D, hyperCube3D, energy_density_evoution, hydrodynamic_evoution,
                                        freezeoutSurfaceFile, fo_surf);


    //if all cells are below freezeout temperature end hydro
    accumulator1 = 0;
    accumulator1 = checkForCellsAboveTc(nx, ny, nz, freezeoutEnergyDensity, e);
    if (accumulator1 == 0) accumulator2 += 1;
    if (accumulator2 >= FOFREQ + 1) //only break once freezeout finder has had a chance to search/write to file
    {
      //make sure we dont break if we are running dynamical sources
      if (initialConditionType != 12)
      {
        printf("\nAll cells have dropped below freezeout energy density\n");
        break;
      }
      else if ( (initialConditionType == 12) && (n > numberOfSourceFiles) )
      {
        printf("\nAll cells have dropped below freezeout energy density and source terms finished \n");
        break;
      }
    }

    t1 = std::clock();

    // Read in source terms from particles
    if (initialConditionType == 12) {
        if (n <= numberOfSourceFiles) readInSource(n, latticeParams, initCondParams, hydroParams, rootDirectory);
        else if (n == numberOfSourceFiles + 1) noSource(latticeParams, initCondParams);
    }

    rungeKutta2(t, dt, q, Q, latticeParams, hydroParams);

    t2 = std::clock();
    double delta_time = (t2 - t1) / (double)(CLOCKS_PER_SEC / 1000);
    if ((n-1) % FREQ == 0) printf("(Elapsed time: %.3f ms)\n",delta_time);
    totalTime+=delta_time;
    ++nsteps;

    setCurrentConservedVariables();

    //calculate the thermal vorticity tensor for use in polarization studies
    #ifdef THERMAL_VORTICITY
    calculateThermalVorticity(t, dt, q, Q, latticeParams, hydroParams);
    #endif

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

  free4dArray(hyperCube4D, 2, 2, 2);
}
