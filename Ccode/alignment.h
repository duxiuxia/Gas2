/*
 *  alignment.h
 *  Janus
 *
 *  Created by mike on 2/24/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */

#ifndef ALIGNMENT_H_INCLUDED
#define ALIGNMENT_H_INCLUDED
#include "common.h"
void alignment(std::vector<ALIGNED_COMPOUNDS> &,std::vector<std::list<COMPOUND> > & ,int ,double ,double ,bool ,double);
void getCompoundList(std::vector<FEATURE_POINT> &featurelist,std::list<COMPOUND> * compounds, unsigned FileID);
void saveAlignedFeatureList(std::vector<ALIGNED_COMPOUNDS> &alignedCompList,std::string outputDir,std::string projectName);

#endif // READCDF_H_INCLUDED
