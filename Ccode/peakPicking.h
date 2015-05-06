#ifndef PEAKPICKING_H_INCLUDED
#define PEAKPICKING_H_INCLUDED
/*
 *  peakPicking.h
 *  Janus
 *
 *  Created by mike on 2/3/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */
#include "common.h"


struct pkDetectPoint {
	int val;
	int counter;
};


void getFeatureList(std::vector<FEATURE_POINT>* vFeature,const std::vector<int> & vValleyIndex, const std::vector<int>& vPeakIndex,const std::vector<int>& vecInt);
void mergeLowerAjacentPk(int i,std::vector<FEATURE_POINT> & vFeature);
void DenoiseByArea(std::vector<FEATURE_POINT>& vFeature,int nLocalWindowSize,int EICPkRatio);
void DenoiseByHeight(std::vector<FEATURE_POINT>& vFeature,unsigned nLocalWindowSize,unsigned EICPkRatio);
void peakBoundary(const int* arrayValleyInd,const int nValley,const int *arrayPkInd,const int nPeak,double * lbound,double *rbound,int *area,int *isRemove,const double *arrayRT,const int *arrayInt,int nEIClength);
void peakDetect(const int sign,const std::vector<int>  x,std::vector<int>& resultInd,unsigned nPkSpan);
void saveFeatureList(FeatureType & featurelist,std::string outputDir,std::string filename);
void loadFeatureList(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename);

bool compPkDetectPoints(pkDetectPoint p1,pkDetectPoint p2);


#endif // READCDF_H_INCLUDED
