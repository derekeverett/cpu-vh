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
	if (cli.runHydro)
		printf("runHydro = True\n");
	else
		printf("runHydro = False\n");
	if (cli.runTest)
		printf("runTest = True\n");
	else
		printf("runTest = False\n");

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
