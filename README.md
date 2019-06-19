The following papers should be cited when referring to cpu-vh:
1) D. Bazow, U. W. Heinz and M. Strickland, Comput. Phys. Commun. 225, 92 (2018) doi:10.1016/j.cpc.2017.01.015
[arXiv:1608.06577 [physics.comp-ph]]
2) (L. Du) in preparation

# cpu-vh
cpu-vh is a code designed for the hydrodynamic simulation of heavy ion collisions.
Please see arXiv:1608.06577 for a description of the physics, as well as the KT algorithm which 
is used to solve the hydrodynamic equations of motion. 
This code is the CPU version of the algorithm described in the paper, and has been further optimized.

The up to date GPU version of this code can be found at https://github.com/derekeverett/gpu-vh
The two codes are designed to be as similar as possible, with the flexibility of running on 
heterogeneous computing platforms.

# Installation
To compile with cmake:
   mkdir build & cd build
   cmake ..
   make
   make install

There should now exist an executable cpu-vh in the parent directory. 

# Usage

To run the code:
   sh run.sh <NUMTHREADS>

where <NUMTHREADS> is the number of cores on which to run the code in parallel.  
Alternatively one can export the enviroment variable `OMP_NUM_THREADS`.

All parameters can be set in the files inside the 'rhic-conf' directory.

'lattice.properties' contains parameters which determine the hydrodynamic grid size, 
number of points, the time step size and number of time steps.
NOTE* The number of points in x, y and eta need to be odd for the grid to be centered!

'ic.properties' contains the parameters for the initial conditions for hydrodynamics. 

'hydro.properties' contains parameters related to the viscous hydrodynamic evolution 
(for example the shear viscosity). 

The freezeout surface file 'surface.dat' is written to the directory `output`; 
this directory must exist at runtime running cpu-vh.
The freezeout file is written in a format readable by the Cooper Frye and 
Sampling Code iS3D : \url{https://github.com/derekeverett/iS3D}.

If reading initial conditions from file, the input files are read from the `input` directory. 


