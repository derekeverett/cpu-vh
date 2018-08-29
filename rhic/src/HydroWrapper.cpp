//This is a c++ wrapper for the OSU hydro codes (cpu-vh, gpu-vh ...)
//Written by Derek Everett 2018

#include "../include/HydroWrapper.h"
#include "../include/HydroInitialTmunu.h"
#include "../include/RuntimeParameters.h"

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

void HYDRO::initialize_from_vectors(std::vector<double> e, //e
                        //std::vector<double> p, //p
                        std::vector<double> ut, //ut
                        std::vector<double> ux, //ux
                        std::vector<double> uy, //uy
                        std::vector<double> un, //un
                        std::vector<double> pitt, //pitt
                        std::vector<double> pitx, //pitx
                        std::vector<double> pity, //pity
                        std::vector<double> pitn, //pitn
                        std::vector<double> pixx, //pixx
                        std::vector<double> pixy, //pixy
                        std::vector<double> pixn, //pixn
                        std::vector<double> piyy, //piyy
                        std::vector<double> piyn, //piyn
                        std::vector<double> pinn, //pinn
                        std::vector<double> Pi) //Pi
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
//useful for running in JETSCAPE or other c++ frameworks
int HYDRO::run_hydro_no_cli()
{
  //set properties of initialization
  HydroInitialTmunu init_tmunu;
  init_tmunu.e_in = initial_energy_density;

  init_tmunu.ut_in = initial_ut;
  init_tmunu.ux_in = initial_ux;
  init_tmunu.uy_in = initial_uy;
  init_tmunu.un_in = initial_un;

  #ifdef PIMUNU
  init_tmunu.pitt_in = initial_pitt;
  init_tmunu.pitx_in = initial_pitx;
  init_tmunu.pity_in = initial_pity;
  init_tmunu.pitn_in = initial_pitn;
  init_tmunu.pixx_in = initial_pixx;
  init_tmunu.pixy_in = initial_pixy;
  init_tmunu.pixn_in = initial_pixn;
  init_tmunu.piyy_in = initial_piyy;
  init_tmunu.piyn_in = initial_piyn;
  init_tmunu.pinn_in = initial_pinn;
  #endif
  
  #ifdef PI
  init_tmunu.Pi_in = initial_Pi;
  #endif

  RuntimeParameters run_params;
  run_params.configDirectory = "rhic-conf";
  run_params.outputDirectory = "cpu-vh_output";
  run_params.runHydro = true;

  //struct CommandLineArguments cli;
	struct LatticeParameters latticeParams;
	struct InitialConditionParameters initCondParams;
	struct HydroParameters hydroParams;

	char *rootDirectory = NULL;
	size_t size;
	rootDirectory = getcwd(rootDirectory,size);

	// Print argument values
	printf("configDirectory = %s\n", run_params.configDirectory);
	printf("outputDirectory = %s\n", run_params.outputDirectory);
	//if (run_params.runTest) printf("runTest = True\n");
	//else printf("runTest = False\n");

  //this (simpler) parser doesnt require libconfig
  readLatticeParameters(run_params.configDirectory, &latticeParams);
  readInitialConditionParameters(run_params.configDirectory, &initCondParams);
  readHydroParameters(run_params.configDirectory, &hydroParams);


  //clock
  double sec = 0.0;
  #ifdef _OPENMP
  sec = omp_get_wtime();
  #endif
	//=========================================
	// Run hydro
	//=========================================
	if (run_params.runHydro) {
		run(&latticeParams, &initCondParams, &hydroParams, rootDirectory, run_params.outputDirectory, init_tmunu);
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
  HydroInitialTmunu init_tmunu;

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
		run(&latticeParams, &initCondParams, &hydroParams, rootDirectory, cli.outputDirectory, init_tmunu);
		printf("Done hydro.\n");
	}
  //clock
  #ifdef _OPENMP
  sec = omp_get_wtime() - sec;
  #endif
  printf("Hydro took %f seconds\n", sec);

	return 0;
}
