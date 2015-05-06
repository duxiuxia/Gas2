/*
 *  common.cpp
 *  Janus
 *
 *  Created by mike on 2/3/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */

#include "common.h"

void initFeature(FEATURE_POINT* Feature)
{
	Feature->isRm_SmArea=0;
	Feature->isRm_SmHeight=0;
	Feature->isRm_ComBoundary=0;
	Feature->lboundInd=-1;
	Feature->rboundInd=-1;
	Feature->area=0;
	Feature->Intensity=0;
	Feature->pkInd=-1;
	Feature->mz=0;
	Feature->compoundID=0;
	Feature->windowID=0;
        //Added by wenchao zhang, add 3 column in Decon.csv file for ADAP2
        Feature->isModel=0;
        Feature->factor =0.0;
        Feature->isUnique=0;    
	
}

void parseFileName (const std::string & filePath,std::string & filename)
{
    int nstart=filePath.find_last_of("/")+1;
    int nend=filePath.find_last_of(".");
    int nsize=nend-nstart;
    filename.append(filePath.substr(nstart,nsize));
}
bool compFeaturePointHeight(FEATURE_POINT p1,FEATURE_POINT p2)
{
	return p1.Intensity<p2.Intensity;
}

bool compFeaturePointArea(FEATURE_POINT p1,FEATURE_POINT p2)
{
	return p1.area<p2.area;
}
bool compFeaturePointRT(FEATURE_POINT p,FEATURE_POINT val)
{
	return p.pkInd<val.pkInd;
}
//bool compFeaturePointRT_low(FEATURE_POINT p,FEATURE_POINT val)
//{
//	return p.pkInd<val.pkInd;
//}

void saveFirstScanRT(std::vector<double> scanIntervalVec,std::vector<double> delaytime,std::vector<std::string> DataFiles,std::string outputDir,std::string projectName)
{
	
	projectName.insert(0,outputDir);
	projectName.append("_RT.csv");
	std::ofstream outputfile;
	outputfile.open(projectName.c_str(), std::ios::out);
	outputfile<<"FileID,FileName,firstRT,ScanInterval"<<std::endl;//
	
	FeatureType::iterator it;
	for(unsigned i=0;i<DataFiles.size();i++)
	{
		//std::cout<<delaytime.at(i)<<std::endl;
		std::string filename;
		parseFileName (DataFiles.at(i),filename);
		outputfile<<(i+1)<<","<<filename<<","<<delaytime.at(i)/60<<","<<scanIntervalVec.at(i)/60<<std::endl;
	}
	outputfile.close();
	
}
void LoadFirstScanRT(std::vector<double> &scanIntervalVec,std::vector<double> &delaytime,std::string outputDir,std::string projectName)
{

	projectName.insert(0,outputDir);
	projectName.append("_RT.csv");

	char dimSym='\n';
	/*read file*/
	std::cout << "read:"<<projectName<<std::endl;
	std::ifstream inputfile;
	inputfile.open(projectName.c_str(), std::ios::in);
	std::string csvLine;
	std::getline(inputfile,csvLine,dimSym);//skip the first header line
	//unsigned i=1;
	while(std::getline(inputfile,csvLine,dimSym))
	{
		std::stringstream ss;
		std::string csvElement;
		std::istringstream csvStream(csvLine);//put current line into string stream
		double curRT,curInterval;
		std::getline(csvStream, csvElement, ',');//read first col
		std::getline(csvStream, csvElement, ',');//read second col

		std::getline(csvStream, csvElement, ',');//read third col, which is the RT


		ss.str(csvElement);
		ss>>curRT;//read the current element 
		ss.clear();				
		delaytime.push_back(curRT);

		std::getline(csvStream, csvElement, ',');//read 4th col, which is the scanInterval
		ss.str(csvElement);
		ss>>curInterval;//read the current element
		ss.clear();
		scanIntervalVec.push_back(curInterval);
	}
	inputfile.close();
	
	
	
}

std::string stoupper(std::string s)
{
	
	
	std::string::iterator i = s.begin();
	std::string::iterator end = s.end();
	
	
	while (i != end) {
		*i = std::toupper((unsigned char)*i);
		++i;
	}
	return s;
	
}

