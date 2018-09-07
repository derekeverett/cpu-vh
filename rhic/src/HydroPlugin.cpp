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
#include "cornelius-c++-1.3/cornelius.cpp"
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
  else printf("simulation is not in 3+1D or 2+1D; freezeout finder will not work!\n");

  Cornelius cor;
  cor.init(dim, freezeoutEnergyDensity, lattice_spacing);

  //declare an array of FO_Element to hold FO cell info
  std::vector<FO_Element> fo_surf;

  double ****energy_density_evoution;
  energy_density_evoution = calloc4dArray(energy_density_evoution, FOFREQ+1, nx, ny, nz);

  //make an array to store all the hydrodynamic variables for FOFREQ time steps
  //to be written to file once the freezeout surface is determined by the critical energy density
  #ifdef THERMAL_VORTICITY
  int n_hydro_vars = 22; //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, wtx, wty, wtn, wxy, wxn, wyn (the temperature and pressure are calclated with EoS)
  #else
  int n_hydro_vars = 16; //u0, u1, u2, u3, e, pi00, pi01, pi02, pi03, pi11, pi12, pi13, pi22, pi23, pi33, Pi, the temperature and pressure are calclated with EoS
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
  if (FOFORMAT == 0) freezeoutSurfaceFile.open(surf_fname);
  else freezeoutSurfaceFile.open(surf_fname, ios::binary);

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
    if (nFO == FOFREQ - 1) //call the freezeout finder should this be put before the values are set?
    {
      //besides writing centroid and normal to file, write all the hydro variables
      int dimZ;
      if (dim == 4) dimZ = nz-1; //enter the loop over iz, and avoid problems at boundary
      else if (dim == 3) dimZ = 1; //we need to enter the 'loop' over iz rather than skipping it
      for (int it = start; it < FOFREQ; it++) //note* avoiding boundary problems (reading outside array)
      {
        for (int ix = 0; ix < nx-1; ix++)
        {
          for (int iy = 0; iy < ny-1; iy++)
          {
            for (int iz = 0; iz < dimZ; iz++)
            {
              //write the values of energy density to all corners of the hyperCube
              if (dim == 4) writeEnergyDensityToHypercube4D(hyperCube4D, energy_density_evoution, it, ix, iy, iz);
              else if (dim == 3) writeEnergyDensityToHypercube3D(hyperCube3D, energy_density_evoution, it, ix, iy);

              //use cornelius to find the centroid and normal vector of each hyperCube
              if (dim == 4) cor.find_surface_4d(hyperCube4D);
              else if (dim == 3) cor.find_surface_3d(hyperCube3D);
              //write centroid and normal of each surface element to file
              for (int i = 0; i < cor.get_Nelements(); i++)
              {
                //declare a new fo cell to hold info, later push back to vector
                FO_Element fo_cell;
                double temp = 0.0; //temporary variable
                //first write the position of the centroid of surface element
                double cell_tau;
                if (n <= FOFREQ) cell_tau = t0 + ( (double)(n - FOFREQ + (it-1) ) )* dt; //check if this is the correct time!
                else cell_tau = t0 + ((double)(n - FOFREQ + it)) * dt; //check if this is the correct time!
                double cell_x = (double)ix * dx  - (((double)(nx-1)) / 2.0 * dx);
                double cell_y = (double)iy * dy  - (((double)(ny-1)) / 2.0 * dy);
                double cell_z = (double)iz * dz  - (((double)(nz-1)) / 2.0 * dz);

                double tau_frac = cor.get_centroid_elem(i,0) / lattice_spacing[0];
                double x_frac = cor.get_centroid_elem(i,1) / lattice_spacing[1];
                double y_frac = cor.get_centroid_elem(i,2) / lattice_spacing[2];
                double z_frac;
                if (dim == 4) z_frac = cor.get_centroid_elem(i,3) / lattice_spacing[3];
                else z_frac = 0.0;

                if (FOFORMAT == 0) //write ASCII file
                {
                  if (FOTEST)
                  {
                    freezeoutSurfaceFile << cell_tau << " ";
                    fo_cell.tau = cell_tau;
                  }
                  else
                  {
                    freezeoutSurfaceFile << cor.get_centroid_elem(i,0) + cell_tau << " ";
                    fo_cell.tau = (cor.get_centroid_elem(i,0) + cell_tau);
                  }
                  freezeoutSurfaceFile << cor.get_centroid_elem(i,1) + cell_x << " ";
                  fo_cell.x = cor.get_centroid_elem(i,1) + cell_x;
                  freezeoutSurfaceFile << cor.get_centroid_elem(i,2) + cell_y << " ";
                  fo_cell.y = cor.get_centroid_elem(i,2) + cell_y;
                  if (dim == 4)
                  {
                    freezeoutSurfaceFile << cor.get_centroid_elem(i,3) + cell_z << " ";
                    fo_cell.eta = cor.get_centroid_elem(i,3) + cell_z;
                  }
                  else
                  {
                    freezeoutSurfaceFile << cell_z << " ";
                    fo_cell.eta = cell_z;
                  }
                  //then the (covariant?) surface normal element; check jacobian factors of tau for milne coordinates!
                  //acording to cornelius user guide, corenelius returns covariant components of normal vector without jacobian factors
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,0) << " ";
                  fo_cell.dat = t * cor.get_normal_elem(i,0);
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,1) << " ";
                  fo_cell.dax = t * cor.get_normal_elem(i,1);
                  freezeoutSurfaceFile << t * cor.get_normal_elem(i,2) << " ";
                  fo_cell.day = t * cor.get_normal_elem(i,2);
                  if (dim == 4)
                  {
                    freezeoutSurfaceFile << t * cor.get_normal_elem(i,3) << " ";
                    fo_cell.dan = t * cor.get_normal_elem(i,3);
                  }
                  else
                  {
                    freezeoutSurfaceFile << 0.0 << " ";
                    fo_cell.dan = 0.0;
                  }
                  //write all the necessary hydro dynamic variables by first performing linear interpolation from values at
                  //corners of hypercube

                  if (dim == 4) // for 3+1D
                  {
                    //first write the contravariant flow velocity
                    for (int ivar = 0; ivar < dim; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the energy density
                    temp = interpolateVariable4D(hydrodynamic_evoution, 4, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile << temp << " "; //note : iSpectra reads in file in fm^x units e.g. energy density should be written in fm^-4
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile << effectiveTemperature(temp) << " ";
                    //the thermal pressure
                    freezeoutSurfaceFile << equilibriumPressure(temp) << " ";
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the bulk pressure Pi
                    temp = interpolateVariable4D(hydrodynamic_evoution, 15, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                    freezeoutSurfaceFile << temp;
                    //write 6 components of w^\mu\nu thermal vorticity tensor
                    #ifdef THERMAL_VORTICITY
                    for (int ivar = 16; ivar < 22; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    #endif
                    freezeoutSurfaceFile << endl;
                  }

                  else //for 2+1D
                  {
                    //first write the flow velocity
                    for (int ivar = 0; ivar < 4; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the energy density
                    temp = interpolateVariable3D(hydrodynamic_evoution, 4, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile << temp << " "; //note units of fm^-4 appropriate for iSpectra reading
                    //the temperature !this needs to be checked
                    freezeoutSurfaceFile << effectiveTemperature(temp) << " ";
                    //the thermal pressure
                    freezeoutSurfaceFile << equilibriumPressure(temp) << " ";
                    //write ten components of pi_(mu,nu) shear viscous tensor
                    for (int ivar = 5; ivar < 15; ivar++)
                    {
                      temp = interpolateVariable3D(hydrodynamic_evoution, ivar, it, ix, iy, tau_frac, x_frac, y_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    //write the bulk pressure Pi
                    temp = interpolateVariable3D(hydrodynamic_evoution, 15, it, ix, iy, tau_frac, x_frac, y_frac);
                    freezeoutSurfaceFile << temp;
                    #ifdef THERMAL_VORTICITY
                    for (int ivar = 16; ivar < 22; ivar++)
                    {
                      temp = interpolateVariable4D(hydrodynamic_evoution, ivar, it, ix, iy, iz, tau_frac, x_frac, y_frac, z_frac);
                      freezeoutSurfaceFile << temp << " ";
                    }
                    #endif
                    freezeoutSurfaceFile << endl;
                  }
                } // if (FOFORMAT == 0)

                //add the fo cell to fo surface
                fo_surf.push_back(fo_cell);

              } //for (int i = 0; i < cor.get_Nelements(); i++)
            } // for (int iz = 0; iz < dimZ; iz++)
          } // for (int iy = 0; iy < ny-1; iy++)
        } // for (int ix = 0; ix < nx-1; ix++)
      } //for (int it = start; it < FOFREQ; it++)
    } // if (nFO == FOFREQ - 1)

    //if all cells are below freezeout temperature end hydro
    accumulator1 = 0;
    for (int ix = 2; ix < nx+2; ix++)
    {
      for (int iy = 2; iy < ny+2; iy++)
      {
        for (int iz = 2; iz < nz+2; iz++)
        {
          int s = columnMajorLinearIndex(ix, iy, iz, nx+4, ny+4, nz+4);
          if (e[s] > freezeoutEnergyDensity) accumulator1 = accumulator1 + 1;
        }
      }
    }
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
