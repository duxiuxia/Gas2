/*
 *  deconvolution.h
 *  Janus
 *
 *  Created by mike on 2/9/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */
#ifndef DECONVOLUTION_H_INCLUDED
#define DECONVOLUTION_H_INCLUDED

#include "common.h"
#include "kmedoid.h"
void deconvolution(std::vector<FEATURE_POINT> & featurelist, EICtype & EIC,unsigned maxScanShift,unsigned FragPkRatio,unsigned nFragments,double maxIntraDist);
void saveFeatureListDecon(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename);
void convertFeatureListDecon(std::vector<FEATURE_POINT> & featurelist, EICtype &EIC,double scanInterval,double delaytime,std::string filename);
bool LoadFeatureListDecon(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename);


#endif // READCDF_H_INCLUDED
