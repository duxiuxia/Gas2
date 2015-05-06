#ifndef READCDF_H_INCLUDED
#define READCDF_H_INCLUDED
#include <netcdfcpp.h>
#include "common.h"
//void readCDF(const char * sourcefile,int & nTotalScans,int & nTotalPoints,double* & arrayRT, int* &arrayCount,int* &arrayInt,int* &arrayMz);
void readCDF(const char * sourcefile,int & nTotalScans,int & nTotalPoints,int* &arrayCount,int* &arrayInt,int* &arrayMz,double &delaytime,double &scanInterval);

#endif // READCDF_H_INCLUDED
