/*
 *  resolveEIC.cpp
 *  Janus
 *
 *  Created by mike on 4/1/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */

#include "resolveEIC.h"

//void resolveEIC(std::vector<FEATURE_POINT> & featurelist, EICtype & EIC)
//{
//	stable_sort(featurelist.begin(),featurelist.end(),compFeaturePointRT);
//	
//	int nWindowCount=0;
//	std::vector<FEATURE_POINT>::iterator lower_it=find_if(featurelist.begin(), featurelist.end(),findZeroCompID);
//	//int ntotal=featurelist.size();
//	while(lower_it!=featurelist.end())
//	{
//		//		double Xvalue =(featurelist.size()-distance(lower_it,featurelist.end()))*100/featurelist.size();
//		//std::cout<<'\r'<<Xvalue<<'%';
//		//std::cout<<".";
//		
//		nWindowCount++;
//		FEATURE_POINT val;
//		int curInd=lower_it->pkInd;
//		int endInd=curInd+maxScanShift;
//		val.pkInd=endInd;
//		std::vector<FEATURE_POINT>::iterator upper_it=upper_bound(lower_it, featurelist.end(), val,compUndetectedFeaturePointRT);
//		
//		std::vector<FEATURE_POINT>::iterator minFP=min_element(lower_it, upper_it,compUndetectedFeaturePointHeight);
//		int minPkHeight=minFP->Intensity;
//		std::vector<FEATURE_POINT*> curWindow;
//		for(std::vector<FEATURE_POINT>::iterator it=lower_it;it!=upper_it;it++)//transpass the local window to set flag
//		{
//			if(it->Intensity/minPkHeight>FragPkRatio)
//			{
//				if(it->compoundID==0)
//					curWindow.push_back(&(*it));
//			}
//			else
//				it->compoundID=-1;//filter by pk area
//		}
//		
//		if(curWindow.size()>=nFragments)
//		{
//			//std::cout<<"window:"<<nWindowCount<<std::endl;
//			std::vector<ClusterSolution> clusterResults;
//			
//			//run clustering and save the cluster results
//			clust(curWindow,EIC,clusterResults);
//			
//			//get the best solution with the maximum si value
//			std::vector<ClusterSolution>::iterator bestClustResult=max_element(clusterResults.begin(), clusterResults.end(),compSi);	
//			
//			//get clusting information from the best solution
//			std::map<unsigned,silhouetteType>* clustInfo= &(bestClustResult->clustInfo);
//			//get points in each cluster
//			std::vector<clustPointType>* clustPoints=&(bestClustResult->point);
//			/*set compoundID for each point according to the number of elements in its cluster*/
//			for(unsigned i=0;i<clustPoints->size();i++)
//			{
//				unsigned curClustID=clustPoints->at(i).clustAssignment;
//				unsigned nCount=(*clustInfo)[curClustID].ncount;
//				if(nCount>=nFragments)
//				{
//					curWindow.at(i)->compoundID=curClustID+1;//increase the original clustID by 1 to avoid the the ClustID with zero value
//					curWindow.at(i)->windowID=nWindowCount;
//				}
//				else
//					curWindow.at(i)->compoundID=-3;//not enough fragments after clustering
//			}
//			
//		}
//		else//set compounID=-2 if not enough fragments after denoising
//		{
//			for_each(curWindow.begin(), curWindow.end(), assignCompIDforCurWindow);	
//		}
//		//move forward to the next undetected feature
//		lower_it=find_if(featurelist.begin(), featurelist.end(),findZeroCompID);
//		
//		
//	}
//	
//}