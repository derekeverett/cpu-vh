/*
 * HydroPlugin.h
 *
 *  Created on: Oct 23, 2015
 *      Author: bazow
 */

#ifndef HYDROPLUGIN_H_
#define HYDROPLUGIN_H_

#include "HydroInitialTmunu.h"

void run(void * latticeParams, void * initCondParams, void * hydroParams, const char *rootDirectory, const char *outputDir, HydroInitialTmunu init_tmunu);

#endif /* HYDROPLUGIN_H_ */
