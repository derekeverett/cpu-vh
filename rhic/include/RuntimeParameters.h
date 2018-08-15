#ifndef RUNTIME_PARAMETERS_H
#define RUNTIME_PARAMETERS_H

#include <stdlib.h>


using namespace std;

class RuntimeParameters {
private:

public:
  //bool runTest;
  bool runHydro;
  char *configDirectory;              /* The -v flag */
  char *outputDirectory;            /* Argument for -o */

  //other options? 
};
#endif
