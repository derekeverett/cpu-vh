#ifndef HYDRO_INITIAL_TMUNU_H
#define HYDRO_INITIAL_TMUNU_H

#include <stdlib.h>
#include <vector>

#include "DynamicalVariables.h"

using namespace std;

class HydroInitialTmunu {
private:

public:
  std::vector<PRECISION> e_in;
  std::vector<PRECISION> ut_in;
  std::vector<PRECISION> ux_in;
  std::vector<PRECISION> uy_in;
  std::vector<PRECISION> un_in;
  std::vector<PRECISION> pitt_in;
  std::vector<PRECISION> pitx_in;
  std::vector<PRECISION> pity_in;
  std::vector<PRECISION> pitn_in;
  std::vector<PRECISION> pixx_in;
  std::vector<PRECISION> pixy_in;
  std::vector<PRECISION> pixn_in;
  std::vector<PRECISION> piyy_in;
  std::vector<PRECISION> piyn_in;
  std::vector<PRECISION> pinn_in;
  std::vector<PRECISION> Pi_in;
};
#endif
