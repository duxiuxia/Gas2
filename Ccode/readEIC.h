#ifndef READEIC_H_INCLUDED
#define READEIC_H_INCLUDED

#include "common.h"

void initEIC(int *array,int n);
void getEIC (EICtype* EIC,int nTotalScans,const int* arrayCount,const int* arrayInt,const int* arrayMz);
void printEIC(EICtype & EIC,int nTotalScans);
int saveEIC(EICtype & EIC,unsigned nTotalScans,std::string outputDir,std::string filename);
int saveEIC1(EICtype & EIC,unsigned nTotalScans,std::string outputDir,std::string filename);

int loadEIC(EICtype & EIC,unsigned & nTotalScans ,unsigned startMz,unsigned endMz,std::string outputDir,std::string filename);
int loadEIC1(EICtype & EIC,unsigned & nTotalScans,std::string outputDir,std::string filename);

void dispose_EIC(EICtype& EIC);

#endif // READEIC_H_INCLUDED
