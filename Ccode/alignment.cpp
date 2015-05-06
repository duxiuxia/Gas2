/*
 *  alignment.cpp
 *  Janus
 *
 *  Created by mike on 2/24/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */

#include "alignment.h"
#include "kmedoid.h"
#if !defined(max)
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

double sumSpecSimilarity(double init,COMPOUND p)
{
	init+=p.specSimilarity;
	return init;
}

unsigned sumCmpPkInd(unsigned init,COMPOUND p)
{
	init+=p.pkInd;
	return init;
}

unsigned sumFeaturePkInd(unsigned init,FEATURE_POINT p)
{
	init+=p.pkInd;
	return init;
}
unsigned sumFeatureWeightedInt(float init,FEATURE_POINT p)
{
	init+=p.weightedInt;
	return init;
}

bool compWinIDandCompID(FEATURE_POINT p,FEATURE_POINT val)
{
	if(p.windowID<val.windowID)
		return 1;
	else
	{
		if(p.windowID==val.windowID)
		{
			return p.compoundID<val.compoundID;
		}
		else
			return 0;
	}
		
}



bool compRT(COMPOUND val,COMPOUND p)
{
	return val.pkInd<p.pkInd;
}

bool compInt(FEATURE_POINT val,FEATURE_POINT p)
{
	return val.Intensity<p.Intensity;
}

bool compMz(FEATURE_POINT val,FEATURE_POINT p)
{
	return val.mz<p.mz;
}


bool compSpecSimilarity(COMPOUND val,COMPOUND p)
{
	return val.specSimilarity<p.specSimilarity;
}
bool compSimilarity(COMPOUND val,COMPOUND p)
{
	return val.score<p.score;
}

//get the earliest compound across all samples as the tentative reference cmp
std::list<COMPOUND>::iterator getRfComp(std::vector<std::list<COMPOUND> > & compoundVec)
{
	
	std::list<COMPOUND>::iterator minRefCompIt=compoundVec.at(0).begin();
	std::list<COMPOUND>::iterator curCompIt;
	int minPkInd=INT_MAX;
	for(std::vector<std::list<COMPOUND> >::iterator cmplistIt=compoundVec.begin();cmplistIt!=compoundVec.end();cmplistIt++)
	{
		if(!cmplistIt->empty())//only check the non-empty compound list
		{
			curCompIt=cmplistIt->begin();//get current front compound
			int curPkInd=curCompIt->pkInd;//get the pk RT index of current front compound
			minRefCompIt=curPkInd<minPkInd?curCompIt:minRefCompIt;//compare with the previous pk RT, choose the earlier compound
			
			//update the previous compound information
			//PreRefCompIt=refCompIt;
			minPkInd=minRefCompIt->pkInd;
		}
	}
	//COMPOUND refComp(*minRefCompIt);//copy the earliest compound content to new object as the reference compound
	//unsigned refFileID=refComp.FileID;//get the fileid from the reference compound
	//compoundVec.at(refFileID).erase(minRefCompIt);//remove the reference compound from the respective compound list
	return minRefCompIt;
}

/*convert the deconvoluted peaklist into compound (component) structure*/
void getCompoundList(std::vector<FEATURE_POINT> &featurelist,std::list<COMPOUND> * compounds,unsigned FileID)
{
	stable_sort(featurelist.begin(),featurelist.end(),compWinIDandCompID);//sort by windowID,compoundID
	
	//pach feature points into compoundlist
	unsigned preWindowID=8888;
	unsigned preCompoundID=8888;	
	unsigned nCompoundCounter=0;//unique compID accross windows
	//COMPOUND * curCompound;
	for(std::vector<FEATURE_POINT>::iterator it=featurelist.begin();it!=featurelist.end();it++)
	{
		unsigned curWindowID=it->windowID;
		unsigned curCompoundID=it->compoundID;	
		if(curWindowID==preWindowID&&curCompoundID==preCompoundID)//if the same winID and compID, then within the same compound group
			compounds->back().fragments.push_back(*it);//simply append the current feature point the current compound group
		else// if either winID or compID differnt from the previous ones, then reach the new compound group
		{
			COMPOUND newComp;
			compounds->push_back(newComp);//create a new compound,push into vector
			nCompoundCounter++;
			compounds->back().ID=nCompoundCounter;//assign unique compID accross windows
			compounds->back().FileID=FileID;
			compounds->back().specSimilarity=0;//init the Similarity
			compounds->back().fragments.push_back(*it);//append the current point 
			preWindowID=curWindowID;
			preCompoundID=curCompoundID;
		}

	}
}

//get the intensity vector by uinioned mz vector(fill in zero for non-existant mz)
void getIntVec(std::vector<FEATURE_POINT> &val, std::vector<FEATURE_POINT>& unionPoints,std::vector<double> &data,bool isWeighted)
{
	for(unsigned i=0;i<unionPoints.size();i++)
	{
		FEATURE_POINT &curFpoint=unionPoints.at(i);
		std::vector<FEATURE_POINT>::iterator it=lower_bound(val.begin(), val.end(),curFpoint,compMz);//find the element >=curent mz

		if(it!=val.end())
		{
			if(it->mz==curFpoint.mz)//if lbound== curMz
				data.at(i)=isWeighted?it->weightedInt:it->Intensity;
			else//if lbound>cur Mz
			{
				if(it!=val.begin()) //if not the beginning
				{
					it--;//move backward one step
					if(it->mz==curFpoint.mz)//if this mz ==curMz
						data.at(i)=isWeighted?it->weightedInt:it->Intensity;
				}
			}
		}
	}
}
//scale the intensity values of the current spectrum
void normIntPoints(std::vector<FEATURE_POINT> & curPoints)
{
	std::vector<FEATURE_POINT>::iterator maxIntPoint=max_element(curPoints.begin(), curPoints.end(),compInt);

	for (unsigned i = 0; i <curPoints.size(); i++)
	{ 
		curPoints.at(i).weightedInt=((double)(curPoints.at(i).Intensity)/(maxIntPoint->Intensity))*999;
	}
}

//dotproduct between two spectrum
double score(std::vector<FEATURE_POINT> & curPoints,std::vector<FEATURE_POINT> & refPoints,unsigned nSimilarityType,bool isPure)	
{
	std::vector<FEATURE_POINT> result(curPoints.size()+refPoints.size());
	std::vector<FEATURE_POINT>::iterator it;
	if(isPure)
		//get intersecion vector of mz
		it=set_intersection(curPoints.begin(), curPoints.end(), refPoints.begin(), refPoints.end(),result.begin(),compMz);
	else//get union vector of mz
		it=set_union(curPoints.begin(), curPoints.end(), refPoints.begin(), refPoints.end(),result.begin(),compMz);

	std::vector<FEATURE_POINT> trimedResult(result.begin(),it);
	if(trimedResult.size()==0)
		return 0;//if not common fragments when do pure score, then return maximum angle
	else
	{
		//get the intensity vector by uinioned mz vector(fill in zero for non-existant mz)
		std::vector<double>data1(trimedResult.size(),0);
		std::vector<double>data2(trimedResult.size(),0);
		bool isWeighted=1;
		getIntVec(curPoints,trimedResult,data1,isWeighted);
		getIntVec(refPoints,trimedResult,data2,isWeighted);

		
		switch (nSimilarityType) 
		{
			case 1:
				return dotProduct(data1,data2);
				break;
			case 2:
				return NistDotProduct(data1,data2);
				break;
			default:
				return 0;
				break;
		}
	}
	
}

//calculcate the distance between two spectrum by combine pure and impure scores
double specDistance(std::vector<FEATURE_POINT> & curPoints,std::vector<FEATURE_POINT> & refPoints,unsigned nSimilarityType)	
{
	//sort feature points by mz;and refPoints has already been sorted ouside of this function
	stable_sort(curPoints.begin(), curPoints.end(),compMz);
	
	double impureScore=score(curPoints,refPoints,nSimilarityType,0);
	double pureScore=score(curPoints,refPoints,nSimilarityType,1);	
	
	/*pure score*/
	
	return impureScore*0.3+pureScore*0.7;
	
}

//check if the component vector in all samples are empty, which indicates the alignment process is finished
bool isEmptyCompVec(std::vector<std::list<COMPOUND> > &compoundVec)
{
       // add some debug variables to determine the maximum sample number limit
       int limit_th =0;
       for(std::vector<std::list<COMPOUND> >::iterator it=compoundVec.begin();it!=compoundVec.end();it++)
	{
		if(it->empty())
                {
                  std::cout<<"The maximum sample number ADAP is" << limit_th<< "\n";
		  return true;
                }
                limit_th++;
	}
	return false;
}

//count total number of components left in the vector, which are unaligned
unsigned countCompVec(std::vector<std::list<COMPOUND> > &compoundVec)
{
	unsigned nCount=0;
	for(std::vector<std::list<COMPOUND> >::iterator it=compoundVec.begin();it!=compoundVec.end();it++)
	{
		nCount+=it->size();

	}
	return nCount;
}
/*component match by the spectrum similarity and RT similarity*/
void compMatch(std::vector<std::list<COMPOUND> > & compoundVec,std::list<COMPOUND>::iterator refCompIt,ALIGNED_COMPOUNDS &curCluster,int Max_Compound_ET_shift,double minSpecSimilarity,double scanInterval)
{
	int maxScanShift=Max_Compound_ET_shift/(scanInterval*60);//the threshold to specify alignment window

	curCluster.avgSpecSimilarity=0;//init the currrent cluster
	curCluster.alignedClustID=INT_MAX;//init the clusterid
	
	refCompIt->specSimilarity=0;//reset the angle for reference compound
	curCluster.refComp=COMPOUND(*refCompIt);//save a reference copy to clust
	curCluster.refCompIt=refCompIt;//save the reference pointer for the future erase
	
	COMPOUND & refComp=curCluster.refComp;
	
	unsigned refFileID=refComp.FileID;//get current reference file ID
	std::vector<FEATURE_POINT> & refPoints=refComp.fragments;//get feature points of reference compound
	stable_sort(refPoints.begin(), refPoints.end(),compMz);//sort the reference compound by fragment mz
	
	//get the upper bound value
	COMPOUND up;
	up.pkInd=refComp.pkInd+maxScanShift;

	//search the non-empty and non-reference-file compound list for the similar component that matches the reference compound 
	for(unsigned i=0;i<compoundVec.size();i++)
	{
		if(i!=refFileID&&!compoundVec.at(i).empty())
		{
			//get the match window 
			std::list<COMPOUND> &curCompList=compoundVec.at(i);
			
			/*the upper_bound only return the last element which is > val, 
			 so there is chance that the first element is the upperbound, 
			 which means no elements in the list satisfy the searching range */
			std::list<COMPOUND>::iterator upper_it=upper_bound(curCompList.begin(), curCompList.end(),up,compRT);
			//std::list<COMPOUND>::iterator lower_it=lower_bound(curCompList.begin(), curCompList.end(),low,compRT);
			std::list<COMPOUND>::iterator lower_it=curCompList.begin();
			//calculate the spec similarity by dot product
			if(lower_it!=curCompList.end()&&upper_it!=curCompList.begin())//if the range is not empty
			{
				for(std::list<COMPOUND>::iterator it=lower_it;it!=upper_it;it++)//loop within the window
				{
					std::vector<FEATURE_POINT> & curPoints=it->fragments;//get all the fragments
					it->specSimilarity=specDistance(curPoints,refPoints,2);	//calculate the spectral similarity
					it->RTdiff=((double)abs(it->pkInd-refComp.pkInd))/maxScanShift;//calcalute the RT difference
					it->score=0.9*it->specSimilarity+0.1*(1-it->RTdiff);//combine the spectral and RT distance
				}
				//get the best match within the range [start,end) 
				std::list<COMPOUND>::iterator bestMatchComp=max_element(curCompList.begin(), upper_it,compSimilarity);
				
				if(bestMatchComp->specSimilarity>=minSpecSimilarity)//check the minmum spec similarity threshold
				{
					curCluster.matchedComps.push_back(*bestMatchComp);//save a copy of the best match into current cluster
					curCluster.matchedCompsIt.push_back(bestMatchComp);//save the pointer of the best match for the future erase
				}
			}
			
		}
	}
	if(curCluster.matchedComps.size()>0)
	{
		//get current average angle
		double init=0;
		//cacluate the average spectral similarity for current reference
		curCluster.avgSpecSimilarity=accumulate(curCluster.matchedComps.begin(), curCluster.matchedComps.end(),init,sumSpecSimilarity)/curCluster.matchedComps.size();
	}
}

//remove the aligned components from the component vector
void eraseComp(ALIGNED_COMPOUNDS & bestCluster,std::vector<std::list<COMPOUND> > & compoundVec)		
{
	//remove the matched compound 
	for(std::vector<std::list<COMPOUND>::iterator >::iterator eraseCmpIt=bestCluster.matchedCompsIt.begin();eraseCmpIt!=bestCluster.matchedCompsIt.end();eraseCmpIt++)
	{
		unsigned fileID=(*eraseCmpIt)->FileID;
		compoundVec.at(fileID).erase(*eraseCmpIt);
	}
	//remove  reference compound 
	unsigned fileID=bestCluster.refComp.FileID;
	compoundVec.at(fileID).erase(bestCluster.refCompIt);

}
//scaling intensity and weighting by mass
void weightIntensity(std::vector<FEATURE_POINT> & curPoints)
{
	double a=0.5;
	double init=0;
	double sumInt=accumulate(curPoints.begin(), curPoints.end(),init,sumFeatureWeightedInt);
	double weightFactor1=1/(a+sumInt-1);
	
	
	double IntPower=1;
	unsigned MassPower=1;
	
	for (std::vector<FEATURE_POINT>::iterator it = curPoints.begin(); it!=curPoints.end(); it++)
	{ 
		double weightFactor2=1/(1+weightFactor1*it->weightedInt);
		it->weightedInt=pow(it->weightedInt,IntPower)*pow(it->mz,MassPower)*weightFactor2;
	}
	
	
}


//only one global moving pointer, earliest header compound across the samples is picked as reference
void alignment(std::vector<ALIGNED_COMPOUNDS> &alignedCompList,std::vector<std::list<COMPOUND> > & compoundVec,int Max_Compound_ET_shift,double minSpecSimilarity1,double minSpecSimilarity2,bool isTwoPhase,double scanInterval)
{
	for(unsigned i=0;i<compoundVec.size();i++)
	{
		std::list<COMPOUND> &curCompoundList=compoundVec.at(i);
		for(std::list<COMPOUND>::iterator it=curCompoundList.begin();it!=curCompoundList.end();it++)
		{
			/*calcluate mean pk RT for each compound*/
			unsigned init=0;
			it->pkInd=accumulate(it->fragments.begin(), it->fragments.end(),init,sumFeaturePkInd)/it->fragments.size();
			/*normalize (or rescale) the peak intensity*/
			normIntPoints(it->fragments);
			/*add mass weight to intensity*/
			weightIntensity(it->fragments);

		}
		//sort the compound list by RT
		curCompoundList.sort(compRT);
	}
	
	unsigned nCluster=0;//aligned group ID
	double nTotal=countCompVec(compoundVec);//count total number of components left in the vector, which are unaligned
	double nAligned=0;
	while(!isEmptyCompVec(compoundVec))		
	{
		ALIGNED_COMPOUNDS tCluster, bestCluster;
		
		//get the earliest compound across all samples as the tentative reference cmp
		std::list<COMPOUND>::iterator refCompIt=getRfComp(compoundVec);
		
		//init the best cluster by matching the tentative reference component with the second threshold
		compMatch(compoundVec,refCompIt,bestCluster,Max_Compound_ET_shift,minSpecSimilarity2,scanInterval);
		
		if(isTwoPhase)
		{
			/*first phase: tentative match in order to find all the candidte reference 
			 save the candidate references to tCluster*/
			compMatch(compoundVec,refCompIt,tCluster,Max_Compound_ET_shift,minSpecSimilarity1,scanInterval);
			
			//second phase:loop through all the candidate reference 
			for(std::vector<std::list<COMPOUND>::iterator>::iterator matchedCompIt=tCluster.matchedCompsIt.begin();matchedCompIt!=tCluster.matchedCompsIt.end();matchedCompIt++)
			{
				ALIGNED_COMPOUNDS curCluster;
				//get matched compounds for current candidate reference
				compMatch(compoundVec,*matchedCompIt,curCluster,Max_Compound_ET_shift,minSpecSimilarity2,scanInterval);
				//update the bestCluster 
				if(curCluster.avgSpecSimilarity>bestCluster.avgSpecSimilarity)
					bestCluster=curCluster;
				
			}
		}
		
		nCluster++;//increase the global counter
		bestCluster.alignedClustID=nCluster;//update it into current cluster
		alignedCompList.push_back(bestCluster);//save it as the seperate clust
				
		eraseComp(bestCluster,compoundVec);	//remove the aligned components from the vector
		
		nAligned+=bestCluster.matchedComps.size()+1;//increase the number of aligned components
		
		std::cout<<std::fixed<<std::setprecision(0)<<100*(nAligned/nTotal)<<"%\r"<<std::flush;	

	}	
	std::cout<<"\n";
}

//write the reference fragments to csv
void outputComp(unsigned alignedClustID,double alignedET,COMPOUND &cmp,std::ofstream & outputfile,double scanInterval,std::vector<double> delaytime,bool isRef)
{
	//each compound
	unsigned FileID=cmp.FileID+1;//compatible with R counter from fileid
	unsigned compoundID=cmp.ID;
	double compoundRT=(double)(cmp.pkInd-1)*scanInterval+delaytime[cmp.FileID];
	//double specSimilarity=cmp.specSimilarity;
	for(std::vector<FEATURE_POINT>::iterator it=cmp.fragments.begin();it!=cmp.fragments.end();++it)
	{
		//each fragments
		double ET=(double)(it->pkInd-1)*scanInterval+delaytime[cmp.FileID];
		outputfile<<ET<<","<<it->lboundInd<<","<<it->rboundInd<<","<<compoundRT<<","<<it->Intensity<<","<<it->area<<","<<compoundID<<","<<FileID<<","<<it->mz<<","<<alignedClustID<<","<<alignedET<<","<<isRef<< ","<< it->isModel<<","<< it->isUnique <<","<< it->factor <<std::endl;
	}
}

void saveAlignedFeatureList(std::vector<ALIGNED_COMPOUNDS> &alignedCompList,std::string outputDir,std::string projectName)
{
	std::cout << "*** start writing aligned FeatureList file.." << std::endl;
	/*read the RT for the first scan of each sample, for the purpose of converting scan number into RT*/
	std::vector<double> delaytime;
	std::vector<double> scanIntervalVec;
	LoadFirstScanRT(scanIntervalVec,delaytime,outputDir,projectName);
	
	/*at this point, sample rate is directly calculated from scaninteveral
	 * and the assume each sample has the same sample rate
	 * later it needs to be extended to the samples with different sample rate
	 */
	double scanInterval=scanIntervalVec.at(0);

	projectName.insert(0,outputDir);
	projectName.append("_aligned.csv");

	std::ofstream outputfile;
	outputfile.open(projectName.c_str(), std::ios::out);
	/*write the header*/
	outputfile<<"ET,lboundInd,rboundInd,compoundET,int,pkArea,CompoundID,FileID,mz,clustID,alignedET,isRef, isModel, isUnique, factor"<<std::endl;//isRm_SmArea,isRm_ComBoundary"<<std::endl;
	for(std::vector<ALIGNED_COMPOUNDS>::iterator clustIt=alignedCompList.begin();clustIt!=alignedCompList.end();clustIt++)
	{		
		//each aligned compounds cluster
		unsigned init=0;
		unsigned alignedPkInd;
		if(clustIt->matchedComps.size()>0)//calculate the average aligned peak index
		{
			unsigned avgPkInd=accumulate(clustIt->matchedComps.begin(),clustIt->matchedComps.end(),init,sumCmpPkInd)/clustIt->matchedComps.size();
			alignedPkInd=(avgPkInd+clustIt->refComp.pkInd)/2;
		}
		else
			alignedPkInd=clustIt->refComp.pkInd;
		double RTinit=0;
		//calculate the aligned RT by average peak index and average delaytime
		double alignedET=(double)(alignedPkInd-1)*scanInterval+accumulate(delaytime.begin(), delaytime.end(), RTinit)/delaytime.size();
		unsigned alignedClustID= clustIt->alignedClustID;
		//write the reference fragments
		outputComp(alignedClustID,alignedET,clustIt->refComp,outputfile, scanInterval, delaytime,1);
		//write all the other aligned compnonents information
		for(std::vector<COMPOUND>::iterator cmpIt=clustIt->matchedComps.begin();cmpIt!=clustIt->matchedComps.end();cmpIt++)
			outputComp(alignedClustID,alignedET,*cmpIt,outputfile, scanInterval, delaytime,0);

	}
	outputfile.close();
	std::cout << "*** SUCCESS writing aligned FeatureList file!" << std::endl;
	
}

