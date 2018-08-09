//This is a c++ wrapper for the OSU hydro codes (cpu-vh, gpu-vh ...)
//Written by Derek Everett 2018


#ifndef _HYDROWRAPPER_SRC
#define _HYDROWRAPPER_SRC

#include <stdlib.h>
#include <stdio.h> // for printf
#include <sys/time.h> // for timing
#include <unistd.h>  // for current working directory
//#include <libconfig.h>
#include <vector>

#include "../include/CommandLineArguments.h"
#include "../include/LatticeParameters.h"
#include "../include/InitialConditionParameters.h"
#include "../include/HydroParameters.h"
#include "../include/HydroPlugin.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

class HYDRO {
 private:

 public:
    HYDRO();
    ~HYDRO();

    std::vector<double> initial_energy_density;
    std::vector<double> initial_pressure;
    std::vector<double> initial_ut;
    std::vector<double> initial_ux;
    std::vector<double> initial_uy;
    std::vector<double> initial_un;
    std::vector<double> initial_pitt;
    std::vector<double> initial_pitx;
    std::vector<double> initial_pity;
    std::vector<double> initial_pitn;
    std::vector<double> initial_pixx;
    std::vector<double> initial_pixy;
    std::vector<double> initial_pixn;
    std::vector<double> initial_piyy;
    std::vector<double> initial_piyn;
    std::vector<double> initial_pinn;
    std::vector<double> initial_Pi;

    //run the hydro until freezeout
    int run_hydro(int, char **);

    //run the hydro for one time step
    //int run_hydro_time_step();

    //save the freezeout surface info
    //void save_fo_history();

    //support to initilialize the energy density from a vector - useful for JETSCAPE
    //note units of argument should be GeV / fm^3
    //then we convert to fm^(-4)
    //void initialize_ed_from_vector(std::vector<double>);


    //support to initilialize all components of T^\mu\nu from vectors - useful for JETSCAPE
    /*
    void initialize_from_vectors(std::vector<double>&, //e
                            std::vector<double>&, //p
                            std::vector<double>&, //ut
                            std::vector<double>&, //ux
                            std::vector<double>&, //uy
                            std::vector<double>&, //un
                            std::vector<double>&, //pitt
                            std::vector<double>&, //pitx
                            std::vector<double>&, //pity
                            std::vector<double>&, //pitn
                            std::vector<double>&, //pixx
                            std::vector<double>&, //pixy
                            std::vector<double>&, //pixn
                            std::vector<double>&, //piyy
                            std::vector<double>&, //piyn
                            std::vector<double>&, //pinn
                            std::vector<double>&); //Pi
    */
};

#endif
