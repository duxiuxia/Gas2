/*
 *  peakPicking.cpp
 *  Janus
 *
 *  Created by mike on 2/3/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */
#include "peakPicking.h"
bool compPkDetectPoints(pkDetectPoint p1,pkDetectPoint p2)
{
	return p1.val<p2.val;
}

/*read peak list from csv*/
void loadFeatureList(std::vector<FEATURE_POINT> & featurelist,std::string outputDir,std::string filename)
{
	std::cout<<"loading peak list ..."<<std::endl;;
	// time_t time1=time(NULL);
	filename.insert(0,"peakpicking/");
	filename.insert(0,outputDir);
	filename.append("PeakList.csv");

    /*read file*/
   	
	std::ifstream inputfile;
	inputfile.open(filename.c_str(), std::ios::in);
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
		ss>>fpoint.pkInd;
		ss.clear();
		
		std::getline(csvStream, csvElement, ',');
		ss.str(csvElement);
		ss>>fpoint.Intensity;
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
		ss>>fpoint.area;
		ss.clear();
	
		featurelist.push_back(fpoint);
		}
	inputfile.close();
	
}

/*save the peak list into csv*/
void saveFeatureList(FeatureType & featurelist,std::string outputDir,std::string filename)
{
	filename.insert(0,"peakpicking/");
	filename.insert(0,outputDir);
	filename.append("PeakList.csv");
	std::ofstream outputfile;
	outputfile.open(filename.c_str(), std::ios::out);
	outputfile<<"curMz,pkInd,Intensity,lboundInd,rboundInd,area"<<std::endl;//isRm_SmArea,isRm_ComBoundary"<<std::endl;

	FeatureType::iterator it;
	for(it=featurelist.begin();it!=featurelist.end();++it)
	{
		unsigned curMz=it->first;
		char sMz[4];
		sprintf(sMz,"%i",curMz);

		std::vector<FEATURE_POINT> * fpoints=&(it->second);
		std::vector<FEATURE_POINT>::iterator it1;
		for(it1=fpoints->begin();it1!=fpoints->end();it1++)
		{
		
			if(!it1->isRm_ComBoundary&&!it1->isRm_SmArea&&!it1->isRm_SmHeight)
				outputfile<<curMz<<","<<it1->pkInd<<","<<it1->Intensity<<","<<it1->lboundInd<<","<<it1->rboundInd<<","<<it1->area<<std::endl;
			//outputfile<<curMz<<","<<it1->pkInd<<","<<it1->lboundInd<<","<<it1->rboundInd<<","<<it1->area<<","<<it1->isRm_SmArea<<","<<it1->isRm_ComBoundary<<std::endl;

		}
		
	}
	
	outputfile.close();
	//std::cout << "*** SUCCESS writing cdf file!" << std::endl;
		
}

/*local maxima detection, if sign<0, then detect the local minima*/
void peakDetect(const int sign,const std::vector<int>  x,std::vector<int>& resultInd,unsigned nPkSpan=3)
{
	/*get the intensity vector size*/
	unsigned n=x.size();
	  
	if ((nPkSpan < 1) | (nPkSpan > n))throw 11;

	std::vector<int> ind;//peak index
	std::vector<pkDetectPoint> points(n);
	std::vector<pkDetectPoint>::iterator it;
	
	/*init point value*/
	for(unsigned i=0;i<n;i++)
	{
		points.at(i).val=sign*x.at(i);//fill the itensity value into point
		points.at(i).counter=0;//init the counter which records how many times each point was detected as maxima within the moving window
	}	
	int nWindows=nPkSpan/2+1;//moving window of half size of span
	for(it=points.begin();it!=points.end();it++)
	{
		if(it+nWindows!=points.end())
		   {
			/*select maxima within window range*/
			std::vector<pkDetectPoint>::iterator maxPoint=max_element(it, it+nWindows,compPkDetectPoints);
			maxPoint->counter++;//increment of maxima counter
		   }
		else
		   break;
	}
	for(unsigned i=0;i<n;i++)
	{
		/*only when the point is detected as maxima half of span (nWindows) times, it is considered as peak apex within this span*/
		if(points.at(i).counter==nWindows)
			resultInd.push_back(i);
	}
	
}
/*search for smallest peak within local window as the local noise, and use SNR to denoise the small peaks within local window*/
void DenoiseByHeight(std::vector<FEATURE_POINT>& vFeature,unsigned nLocalWindowSize,unsigned EICPkRatio)
{
	if(nLocalWindowSize<600)
		throw 2;
	int minInd=vFeature.front().pkInd;
	int maxInd=vFeature.back().pkInd;
	unsigned nWindows=(maxInd-minInd)/nLocalWindowSize;
	for(unsigned i=0;i<nWindows;i++)
	{
		FEATURE_POINT val;//create  feature point object to assign the boundary information
		
		int lboundInd=i*nLocalWindowSize;//get the lbound index
		val.pkInd=lboundInd;// assign it to the object
		std::vector<FEATURE_POINT>::iterator lbound=lower_bound(vFeature.begin(), vFeature.end(),val,compFeaturePointRT);//get the lbound iterator 
		
		int rboundInd=(i+1)*nLocalWindowSize;
		val.pkInd=rboundInd;
		std::vector<FEATURE_POINT>::iterator rbound=upper_bound(vFeature.begin(), vFeature.end(),val,compFeaturePointRT)-1;//get the rbound iterator 
		if(rbound>lbound)
		{
			std::vector<FEATURE_POINT>::iterator minFP=min_element(lbound, rbound,compFeaturePointHeight);//get the min area
			int minPkHeight=minFP->Intensity;
			for(std::vector<FEATURE_POINT>::iterator it=lbound;it!=rbound+1;it++)//transpass the local window to set flag
				if(it->Intensity/minPkHeight<(int)EICPkRatio)
					it->isRm_SmHeight=1;
		}
		else//if there is no more than two pk in the local window skip to next without denoising.
			continue;
	}
}
void DenoiseByArea(std::vector<FEATURE_POINT>& vFeature,int nLocalWindowSize,int EICPkRatio)
{
	if(nLocalWindowSize<600)
		throw 2;
	int minInd=vFeature.front().pkInd;
	int maxInd=vFeature.back().pkInd;
	int nWindows=(maxInd-minInd)/nLocalWindowSize;
	for(int i=0;i<nWindows;i++)
	{
		FEATURE_POINT val;//create  feature point object to assign the boundary information
		
		int lboundInd=i*nLocalWindowSize;//get the lbound index
		val.pkInd=lboundInd;// assign it to the object
		std::vector<FEATURE_POINT>::iterator lbound=lower_bound(vFeature.begin(), vFeature.end(),val,compFeaturePointRT);//get the lbound iterator 
		
		int rboundInd=(i+1)*nLocalWindowSize;
		val.pkInd=rboundInd;
		std::vector<FEATURE_POINT>::iterator rbound=upper_bound(vFeature.begin(), vFeature.end(),val,compFeaturePointRT)-1;//get the rbound iterator 
		if(rbound>lbound)
		{
			std::vector<FEATURE_POINT>::iterator minFP=min_element(lbound, rbound,compFeaturePointArea);//get the min area
			int minPkArea=minFP->area;
			for(std::vector<FEATURE_POINT>::iterator it=lbound;it!=rbound+1;it++)//transpass the local window to set flag
				if(it->area/minPkArea<EICPkRatio)
					it->isRm_SmArea=1;
		}
		else//if there is no more than two pk in the local window skip to next without denoising.
			continue;
	}
	}

/*remove the peak apex if it share the boundary with another bigger peak apex*/
void mergeLowerAjacentPk(int i,std::vector<FEATURE_POINT> & vFeature)
{
	
	if(i>0)
	{
		if(vFeature.at(i-1).lboundInd==vFeature.at(i).lboundInd)
		{
			int lowerScan=vFeature.at(i-1).Intensity<vFeature.at(i).Intensity?(i-1):i;
			vFeature.at(lowerScan).isRm_ComBoundary=1;
		}
		
	}
	
}

/*detect boundary and save the peak meta info into feature point vector
 size of vFeature vector =size of peak index vector*/
void getFeatureList(std::vector<FEATURE_POINT>* vFeature,const std::vector<int> & vValleyIndex, const std::vector<int>& vPeakIndex,const std::vector<int>& vecInt)
{
	
	int nPeak=vPeakIndex.size();
	int nValley=vValleyIndex.size();
	int j=0;
	int nScan=vecInt.size();
		
	for(int i=0; i<nPeak;i++)
	{
		/*each featuer save the info of one whole peak profile (mz,apex time,height,edges and area)*/
		FEATURE_POINT * curFeature=&vFeature->at(i);
		if(j==nValley-1)//check if reach the end of valley vector
		{
			curFeature->lboundInd=vValleyIndex.at(j-1);
			curFeature->rboundInd=nScan-1;
		}
		else
		{
			/*search for left closest valley as left bound*/	
			for(;j<nValley;j++)
			{
				if(vValleyIndex.at(j)<vPeakIndex.at(i)) //if current valley index less than peak index,move forward to find the next closed valley
				{
					if(j==nValley-1)//when reach the valley array upbound
					{
						curFeature->lboundInd=vValleyIndex.at(j-1);
						curFeature->rboundInd=nScan-1;
						break; //finish the current peak ,jump to next peak
					}
					continue;//move forward along the valley array until the valley index larger than the current peak ind
				}
				else //if current valley index larger than peak index,then rbound detected
				{
					if(j>0)
					{
						curFeature->lboundInd=vValleyIndex.at(j-1);
						curFeature->rboundInd=vValleyIndex.at(j);
					}
					else//if the first valley is the right bound of the current peak, 
					{
						curFeature->lboundInd=0;;//set the start point of EIC array as the lbound
						curFeature->rboundInd=vValleyIndex.at(j);
					}
					break;//finish the boundary detection for current peak and jump to the next peak
				}
			}
			
		}
		curFeature->pkInd=vPeakIndex.at(i);
		curFeature->Intensity=vecInt.at(vPeakIndex.at(i));
		curFeature->area=accumulate(vecInt.begin()+curFeature->lboundInd, vecInt.begin()+curFeature->rboundInd,0);
		mergeLowerAjacentPk(i,*vFeature);
	}
	
}


