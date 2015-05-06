#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED
/*
 *  common.h
 *  Janus
 *
 *  Created by mike on 2/3/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */
#include <time.h>
#include <string>
#include <vector>
#include <list>
#include <map>
#include <stdio.h>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <sstream>
#include <math.h>
#include <cstring>
#include <iomanip>
//#include "mxml.h"

//#include <valarray>
std::string stoupper(std::string s);
struct PARAMS {
	std::string JobName;
	std::string Owner;
	unsigned JobID;
	std::string WorkDir;
	std::string DataDir;
	std::string outputDir;
	std::vector<std::string> DataFiles;
	unsigned LocalWindowSize;
	unsigned EICPkRatio;
	unsigned MaxScanShift;
	unsigned FragPkRatio;
	unsigned Fragments;
	unsigned MaxIntraDist;
	bool TwoPhase;
	bool isSplitEIC;
	bool isSmoothing;
	bool isBaseline;
	unsigned nNode;
	bool isPeakpicking;
	bool isDeconvolution;
	bool isAlignment;
	//unsigned SampleRate;
	unsigned MaxCompoundETShift;
	double MinSpecSimilarity1;
	double MinSpecSimilarity2;
};
struct FEATURE_POINT {
	int lboundInd;
	int rboundInd;
	int pkInd;
	int Intensity;
	double weightedInt;
	int area;
	unsigned mz;
	bool isRm_ComBoundary;
	bool isRm_SmArea;
	bool isRm_SmHeight;
	//int detectCode;
	int windowID;
	int compoundID;
        // Added by wenchao zhang for three columns in Decon.csv file for ADAP2
        bool isModel;
        double factor;
        bool   isUnique;
};
struct COMPOUND {
	std::vector<FEATURE_POINT> fragments;
	unsigned pkInd;
	unsigned FileID;
	unsigned ID;
	//double specAngle;
	double specSimilarity;
	double RTdiff;//RT difference between current component and reference
	double score;
};


struct EIC_POINT {
	double RT;
	int Intensity;
};

void initFeature(FEATURE_POINT* Feature);
typedef std::map<unsigned,int*> EICtype;
typedef std::map<unsigned,std::vector<FEATURE_POINT> > FeatureType;

struct ALIGNED_COMPOUNDS {
	int alignedClustID;
	COMPOUND refComp;
	std::list<COMPOUND>::iterator refCompIt;
	//double avgSpecAngle;
	double avgSpecSimilarity;
	std::vector<COMPOUND> matchedComps;
	std::vector<std::list<COMPOUND>::iterator > matchedCompsIt;
};
void parseFileName (const std::string & filePath,std::string &filename);
bool compFeaturePointHeight(FEATURE_POINT p1,FEATURE_POINT p2);

bool compFeaturePointArea(FEATURE_POINT p1,FEATURE_POINT p2);
bool compFeaturePointRT(FEATURE_POINT p,FEATURE_POINT val);
void saveFirstScanRT(std::vector<double> scanIntervalVec,std::vector<double> delaytime,std::vector<std::string> DataFiles,std::string outputDir,std::string projectName);
void LoadFirstScanRT(std::vector<double> &scanIntervalVec,std::vector<double> &delaytime,std::string outputDir,std::string projectName);

#define Version_ADAP_2 1 //Use ADAP2 Version  
#endif // READCDF_H_INCLUDED
