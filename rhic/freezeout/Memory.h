#ifndef MEMORY_H_
#define MEMORY_H_

double ** calloc2dArray(double **array, int dim1, int dim2);

double *** calloc3dArray(double ***array, int dim1, int dim2, int dim3);

double **** calloc4dArray(double ****array, int dim1, int dim2, int dim3, int dim4);

double ***** calloc5dArray(double *****array, int dim1, int dim2, int dim3, int dim4, int dim5);

void free2dArray(double **array, int dim1);

void free3dArray(double ***array, int dim1, int dim2);

void free4dArray(double ****array, int dim1, int dim2, int dim3);

void free5dArray(double *****array, int dim1, int dim2, int dim3, int dim4);

#endif
