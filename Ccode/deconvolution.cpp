/*
 *  deconvolution.cpp
 *  Janus
 *
 *  Created by mike on 2/9/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */

#include "deconvolution.h"
bool findZeroCompID(FEATURE_POINT p)
{
	return p.compoundID==0;
	
}
//std::pair< unsigned,silhouetteType> 
double addSi(double initVal,std::pair< unsigned,silhouetteType> p)
{
	 //initVal.second.si+= p.second.si;
	initVal+= p.second.si;
	return initVal;
}
bool compUndetectedFeaturePointRT(FEATURE_POINT val,FEATURE_POINT p)
{
	return val.pkInd<p.pkInd&&p.compoundID==0;
}
bool compUndetectedFeaturePointArea(FEATURE_POINT val,FEATURE_POINT p)
{
	return val.area<p.area&&p.compoundID==0;
}
bool compUndetectedFeaturePointHeight(FEATURE_POINT val,FEATURE_POINT p)
{
	return val.Intensity<p.Intensity&&p.compoundID==0;
}
bool compSi(ClusterSolution val,ClusterSolution p)
{
	//std::pair< unsigned,silhouetteType> s1,s2;
	//s1.second.si=0;
	//s2.second.si=0;
	double init=0;
	double s1=accumulate(val.clustInfo.begin(), val.clustInfo.end(), init,addSi);
	double s2=accumulate(p.clustInfo.begin(), p.clustInfo.end(), init,addSi);
	
	//return (s1.second.si/val.clustInfo.size())<(s2.second.si/p.clustInfo.size());
	return (s1/val.clustInfo.size())<(s2/p.clustInfo.size());
}

void assignCompIDforCurWindow(FEATURE_POINT * p)
{
	p->compoundID=-2;
}




//run clustering and save the cluster results
void clust(std::vector<FEATURE_POINT*> curWindow, EICtype & EIC,std::vector<ClusterSolution>& clusterResults,double maxIntraDist)
{
	unsigned nPoints=curWindow.size();
	
	double** distmatrix=distancematrix(curWindow,EIC);
	unsigned npass=500;
	unsigned *clusterid=new unsigned[nPoints];
	double error;
	int ifound;
	
	double curIntraDist=90;//init the intraDist with maximum value
	for (unsigned k=1;k<=3&&curIntraDist>maxIntraDist;k++)//avoid splitting when the intraDist is less than threshold
	{
		kmedoids(k,nPoints,distmatrix,npass,clusterid, &error, &ifound);
		
		//fill the result into struct
		ClusterSolution curResult;
		for(unsigned i=0;i<nPoints;i++)
		{
			clustPointType curPoint;
			curPoint.clustAssignment=clusterid[i];
			curPoint.neigbourDist=DBL_MAX;//init the neighbour distance as the maximum possible value 
			curResult.point.push_back(curPoint);
		}
		
		silhouette(curResult,distmatrix);//calcluate the silhouette
		clusterResults.push_back(curResult);
		//update current intra_dist
		curIntraDist=curResult.clustInfo.begin()->second.ai;
	}

	//garbage collection
	for (unsigned i = 0; i < nPoints; i++)
		delete [] distmatrix[i];
	delete [] distmatrix;
}

void deconvolution(std::vector<FEATURE_POINT> & featurelist, EICtype & EIC,unsigned maxScanShift,unsigned FragPkRatio,unsigned nFragments,double maxIntraDist)
{
	/*sort the peak list by RT so that component window can be quickly detected by binary search*/
	stable_sort(featurelist.begin(),featurelist.end(),compFeaturePointRT);

	int nWindowCount=0;//component window
	/*search for the earliest peak with component ID as 0, which indicates it has not been processed yet
	 the lower_it store the lower bound of the current component window	 */
	std::vector<FEATURE_POINT>::iterator lower_it=find_if(featurelist.begin(), featurelist.end(),findZeroCompID);
	while(lower_it!=featurelist.end())//continue to process until the end of peak list
	{
		
		nWindowCount++;// increment of component window ID for each new window
		FEATURE_POINT val;
		int curInd=lower_it->pkInd;//start index of current RT window
		int endInd=curInd+maxScanShift;// the end index of current RT window
		val.pkInd=endInd;
		
		/*get upper bound the unprocessed peaks within current window range*/
		std::vector<FEATURE_POINT>::iterator upper_it=upper_bound(lower_it, featurelist.end(), val,compUndetectedFeaturePointRT);
		/*search the smallest peak within the current window*/
		std::vector<FEATURE_POINT>::iterator minFP=min_element(lower_it, upper_it,compUndetectedFeaturePointHeight);
		/*denoised the current window by SNR*/
		int minPkHeight=minFP->Intensity;
		std::vector<FEATURE_POINT*> curWindow;
		for(std::vector<FEATURE_POINT>::iterator it=lower_it;it!=upper_it;it++)//transpass the local window to set flag
		{
			if(it->Intensity/minPkHeight>(int)FragPkRatio)
			{
				if(it->compoundID==0)
					curWindow.push_back(&(*it));
			}
			else
				it->compoundID=-1;//filter by pk area
		}
		/*if the peak numbers is larger than threshold after denoising, then perform clustering to group peaks into different clusters*/
		if(curWindow.size()>=nFragments)
		{
			//the vector to store differenet clustering solution resulted by diffent k value from k-medoid
			std::vector<ClusterSolution> clusterResults;

			//run clustering and save the cluster results
			clust(curWindow,EIC,clusterResults,maxIntraDist);
			
			//get the best solution with the maximum si value
			std::vector<ClusterSolution>::iterator bestClustResult=max_element(clusterResults.begin(), clusterResults.end(),compSi);	
			
			//get clusting information from the best solution
			std::map<unsigned,silhouetteType>* clustInfo= &(bestClustResult->clustInfo);
			//get points in each cluster
			std::vector<clustPointType>* clustPoints=&(bestClustResult->point);
			/*set compoundID for each group*/
			for(unsigned i=0;i<clustPoints->size();i++)
			{
				unsigned curClustID=clustPoints->at(i).clustAssignment;
				unsigned nCount=(*clustInfo)[curClustID].ncount;
				
				if(nCount>=nFragments)
				{
					curWindow.at(i)->compoundID=curClustID+1;//increase the original clustID by 1 to avoid the the ClustID with zero value
					curWindow.at(i)->windowID=nWindowCount;//set the windowID
				}
				else
					curWindow.at(i)->compoundID=-3;//not enough fragments after clustering
			}

		}
		else//set compounID=-2 if not enough fragments after denoising
		{
			for_each(curWindow.begin(), curWindow.end(), assignCompIDforCurWindow);	
		}
		//move forward to the next undetected feature
		lower_it=find_if(featurelist.begin(), featurelist.end(),findZeroCompID);
		

	}

}

/*save the deconvolued peak list into csv*/
void saveFeatureListDecon(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename)
{
	std::cout << " writing deconvoluted peak list file" <<filename<< std::endl;

	filename.insert(0,"decon/");
	filename.insert(0,outputDir);
	filename.append("Decon.csv");
	
	std::ofstream outputfile;
	outputfile.open(filename.c_str(), std::ios::out);
	outputfile<<"mz,lboundInd,rboundInd,pkInd,Intensity,area,windowID,compoundID"<<std::endl;//isRm_SmArea,isRm_ComBoundary"<<std::endl;
	std::vector<FEATURE_POINT>::iterator it;
	for(it=featurelist.begin();it!=featurelist.end();++it)
	{

		//if(it->compoundID>0)
				outputfile<<it->mz<<","<<it->lboundInd<<","<<it->rboundInd<<","<<it->pkInd<<","<<it->Intensity<<","<<it->area<<","<<it->windowID<<","<<it->compoundID<<std::endl;
		
	}
	
	outputfile.close();
	
}

//load feature points from all files into single featurelist
bool LoadFeatureListDecon(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename)
{
	
	std::cout<<"loading deconvoluted peak list:"<<filename<<std::endl;
	// time_t time1=time(NULL);
	filename.insert(0,"decon/");
	filename.insert(0,outputDir);
	filename.append("Decon.csv");
	
    /*read file*/
   	
	std::ifstream inputfile;
	inputfile.open(filename.c_str(), std::ios::in);
        if(!inputfile.is_open())
        {
           std::cout << "Load feature fail! File" << filename <<"isn't exist" << std::endl;
           return false;
        }
	std::string csvLine;
	std::getline(inputfile,csvLine);//skip the first header line
	while(std::getline(inputfile,csvLine))
	{
		std::stringstream ss;
		std::string csvElement;
		std::istringstream csvStream(csvLine);//put current line into string stream
		FEATURE_POINT fpoint;
		initFeature(&fpoint);
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.mz;//read the current element 
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.lboundInd;
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.rboundInd;
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.pkInd;
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.Intensity;
		ss.clear();		
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.area;
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.windowID;
		ss.clear();

		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.compoundID;
		ss.clear();

                if(Version_ADAP_2)
                {
                   std::getline(csvStream, csvElement, ',');
		   ss.str(csvElement);
		   ss>>fpoint.isModel;
		   ss.clear();

                   std::getline(csvStream, csvElement, ',');
		   ss.str(csvElement);
		   ss>>fpoint.isUnique;
		   ss.clear();

                   std::getline(csvStream, csvElement, ',');
		   ss.str(csvElement);
		   ss>>fpoint.factor;
		   ss.clear();                   
                }
		
		if(fpoint.compoundID>0)
			featurelist.push_back(fpoint);
	}
	inputfile.close();
	return true;
}




