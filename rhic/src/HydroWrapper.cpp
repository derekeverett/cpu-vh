//This is a c++ wrapper for the OSU hydro codes (cpu-vh, gpu-vh ...)
//Written by Derek Everett 2018

#include "../include/HydroWrapper.h"

const char *version = "";
const char *address = "";

HYDRO::HYDRO()
{
}

HYDRO::~HYDRO()
{
}

void HYDRO::initialize_from_vector(std::vector<double> e)
{
  initial_energy_density = e;
}

void HYDRO::initialize_from_vectors(std::vector<double>& e, //e
                        //std::vector<double>& p, //p
                        std::vector<double>& ut, //ut
                        std::vector<double>& ux, //ux
                        std::vector<double>& uy, //uy
                        std::vector<double>& un, //un
                        std::vector<double>& pitt, //pitt
                        std::vector<double>& pitx, //pitx
                        std::vector<double>& pity, //pity
                        std::vector<double>& pitn, //pitn
                        std::vector<double>& pixx, //pixx
                        std::vector<double>& pixy, //pixy
                        std::vector<double>& pixn, //pixn
                        std::vector<double>& piyy, //piyy
                        std::vector<double>& piyn, //piyn
                        std::vector<double>& pinn, //pinn
                        std::vector<double>& Pi) //Pi
{
  initial_energy_density = e;
  //initial_pressure = p;
  initial_ut = ut;
  initial_ux = ux;
  initial_uy = uy;
  initial_un = un;
  initial_pitt = pitt;
  initial_pitx = pitx;
  initial_pity = pity;
  initial_pitn = pitn;
  initial_pixx = pixx;
  initial_pixy = pixy;
  initial_pixn = pixn;
  initial_piyy = piyy;
  initial_piyn = piyn;
  initial_pinn = pinn;
  initial_Pi = Pi;
}

//this function does not accept command line arguments
int HYDRO::run_hydro_no_cli()
{
  struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);

	// Print argument values
	printf("configDirectory = %s\n", cli.configDirectory);
	printf("outputDirectory = %s\n", cli.outputDirectory);
	if (cli.runHydro) printf("runHydro = True\n");
	else printf("runHydro = False\n");
	if (cli.runTest) printf("runTest = True\n");
	else printf("runTest = False\n");

  //clock
  double sec = 0.0;
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
	//=========================================
	// Run hydro
	//=========================================
	if (cli.runHydro) {
		run(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory);
		printf("Done hydro.\n");
	}
  //clock
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  printf("Hydro took %f seconds\n", sec);
	return 0;
}

//this function accepts command line arguments
int HYDRO::run_hydro(int argc, char **argv)
{
  struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	loadCommandLineArguments(argc, argv, &cli, version, address);

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);

	// Print argument values
	printf("configDirectory = %s\n", cli.configDirectory);
	printf("outputDirectory = %s\n", cli.outputDirectory);
	if (cli.runHydro) printf("runHydro = True\n");
	else printf("runHydro = False\n");
	if (cli.runTest) printf("runTest = True\n");
	else printf("runTest = False\n");

	//=========================================
	// Set parameters from configuration files
	//=========================================

  //this parser requires libconfig
  /*
	config_t latticeConfig, initCondConfig, hydroConfig;
	// Set lattice parameters from configuration file
	config_init(&latticeConfig);
	loadLatticeParameters(&latticeConfig, cli.configDirectory, &latticeParams);
	config_destroy(&latticeConfig);
	// Set initial condition parameters from configuration file
	config_init(&initCondConfig);
	loadInitialConditionParameters(&initCondConfig, cli.configDirectory, &initCondParams);
	config_destroy(&initCondConfig);
	// Set hydrodynamic parameters from configuration file
	config_init(&hydroConfig);
	loadHydroParameters(&hydroConfig, cli.configDirectory, &hydroParams);
	config_destroy (&hydroConfig);
  */

  //this (simpler) parser doesnt require libconfig
  readLatticeParameters(cli.configDirectory, &latticeParams);
  readInitialConditionParameters(cli.configDirectory, &initCondParams);
  readHydroParameters(cli.configDirectory, &hydroParams);

  //clock
  double sec = 0.0;
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
	//=========================================
	// Run hydro
	//=========================================
	if (cli.runHydro) {
		run(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory);
		printf("Done hydro.\n");
	}
  //clock
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  printf("Hydro took %f seconds\n", sec);
	// TODO: Probably should free host memory here since the freezeout plugin will need
	// to access the energy density, pressure, and fluid velocity.
	return 0;
}
