#include "readEIC.h"
#include "readCDF.h"
#include "peakPicking.h"
#include "deconvolution.h"
#include "alignment.h"
#include <stdio.h>
#include <stdlib.h>  /* The standard C libraries */
#include <sys/stat.h>
#include <mpi.h>
#include <sqlite3.h>
#include <dirent.h>


void checkOutputFolder(std::string dir)
{
	std::string outputDir=dir.append("/output");
	mkdir(outputDir.c_str(), S_IRWXU);
	mkdir((outputDir+"/decon").c_str(), S_IRWXU);
	mkdir((outputDir+"/EIC").c_str(), S_IRWXU);
	mkdir((outputDir+"/peakpicking").c_str(), S_IRWXU);
	mkdir((outputDir+"/TIC").c_str(), S_IRWXU);
    //  mkdir((outputDir+"/EIC_Denoise").c_str(), S_IRWXU);
    //	mkdir((outputDir+"/TIC_Denoise").c_str(), S_IRWXU);
}


unsigned splitEIC(std::string DataFile,EICtype& EIC,std::string outputDir,double &delaytime,double &scanInterval)//read CDF ,split into EIC and save EIC
{
    int nTotalScans,nTotalPoints;
    std::string filename;
	parseFileName (DataFile,filename);
	
	//double *arrayRT;
    int *arrayCount;
    int *arrayInt;
    int *arrayMz;
	
	readCDF(DataFile.c_str(),nTotalScans,nTotalPoints,arrayCount,arrayInt,arrayMz,delaytime,scanInterval);
	getEIC(&EIC,nTotalScans,arrayCount,arrayInt,arrayMz);
	saveEIC1(EIC,nTotalScans,outputDir,filename.c_str());//could be turn off since it is time consuming to output the entire EIC

	delete []arrayCount;
	delete []arrayInt;
	delete []arrayMz;
	//delete []arrayRT;
	
	return nTotalScans;

}
void peakpicking(std::string DataFile,EICtype & EIC,const unsigned nTotalScans,PARAMS myParams)//load EIC,peak detect and edge detect, save the result
{
	/*get the file name from the absolute path*/
	std::string filename;
	parseFileName (DataFile,filename);
		
	time_t time1;
	time(&time1);
	
	EICtype::iterator EICiter;
	FeatureType featurelists;
	
	unsigned pkSpan=9;//maxima detection span,each peak apex has to be the maxima within this consecutive span
	unsigned valleySpan=3;//minima detection span
	
	std::cout<<"start peak picking...";
	
	/*loop through each EIC (mz)*/
	for(EICiter=EIC.begin();EICiter!=EIC.end();++EICiter)
	{
		unsigned curMz=EICiter->first;
//		std::cout << "curMZ=" << curMz << std::endl;
//		if(curMz ==800)
//			std::cout<<"Notice curMz 800"<<std::endl;
		std::vector<int> vPeakIndex,vValleyIndex;//vector of peak and value index
		std::vector<int> vecInt(EICiter->second,EICiter->second+nTotalScans);//vector of EIC intensity
		peakDetect(1,vecInt,vPeakIndex,pkSpan);//peak detect
		peakDetect(-1,vecInt,vValleyIndex,valleySpan);//valley detect
		
		// Added wenchao
		if (vValleyIndex.size()<2) continue;
		// Add over
		FEATURE_POINT Feature;
		initFeature(&Feature);//initialize the feature value
		featurelists[curMz]=std::vector<FEATURE_POINT>(vPeakIndex.size(),Feature);//add new element to feature list map
		std::vector<FEATURE_POINT> *vFeature=&featurelists[curMz];
		
		//detect edge and calculate area
		getFeatureList(vFeature,vValleyIndex,vPeakIndex,vecInt);
		//DenoiseByArea(*vFeature,nLocalWindowSize,EICPkRatio);
		DenoiseByHeight(*vFeature,myParams.LocalWindowSize,myParams.EICPkRatio);
		
	}
	std::cout<<"peak picking is done!"<<std::endl;
	time_t time2;
	time(&time2);
	std::cout<<difftime(time2, time1)<<std::endl;
	/*save the peak list into csv*/
	saveFeatureList(featurelists,myParams.outputDir,filename.c_str());
		
}
void decon(std::string DataFile,EICtype & EIC,PARAMS myParams)//load peakpicking results and perform deconvolution, save the results
{
	std::string filename;
	parseFileName (DataFile,filename);
	
	std::vector<FEATURE_POINT> featurelist;
	/*read peak list from csv*/
	loadFeatureList(featurelist,myParams.outputDir,filename);
	
	
	time_t time1;
	time(&time1);
	std::cout<<"start deconvolution..."<<std::endl;
	deconvolution(featurelist,EIC,myParams.MaxScanShift,myParams.FragPkRatio,myParams.Fragments,myParams.MaxIntraDist);
	/*save the deconvolued peak list into csv*/
	saveFeatureListDecon(featurelist,myParams.outputDir,filename);

		
	time_t time2;
	time(&time2);
	
	std::cout<<"deconvolution is done!"<<std::endl;
	std::cout<<difftime(time2, time1)<<std::endl;
}
unsigned count_filenum(std::string Folder_path, std::vector<std::string> &Decon_File_Vec, std::string extend_name)
{
   unsigned file_count = 0;
   DIR * dirp;
   struct dirent * entry;

   if((dirp = opendir(Folder_path.c_str()))==NULL)
   {
     std::cout<< "can't open file folder:" << Folder_path << std::endl; 
     closedir(dirp);
     return 0;  
   } 

   while ((entry = readdir(dirp))!= NULL) 
   {
     std::string file_name = entry->d_name;
     std::string Trim_file_name;
     size_t found= file_name.find(extend_name);
     if((entry->d_type == DT_REG)&&(found!= std::string::npos)) /* If the entry is a regular file */
     { 
        file_count++;
        Trim_file_name =file_name.substr(0,found);
        Decon_File_Vec.push_back(Trim_file_name);
     }
   }
   closedir(dirp);
   return file_count;
}

bool Find_InVector(std::vector<std::string> &Decon_File_Vec, std::string File_Element)
{
   std::vector<std::string>::iterator it;
   for(it=Decon_File_Vec.begin(); it< Decon_File_Vec.end(); it++)
   {
      if(0==File_Element.compare(*it))
       return true;
   }
   return false;
}
void align(PARAMS &myParams,std::vector<double>& scanIntervalVec)//load peakpicking results and perform deconvolution, save the results
{
	unsigned nFile=myParams.DataFiles.size();
        // Add some codes for considering xx_decon.CSV missing case
        std::vector<std::string>Decon_File_Vec;
        unsigned nFile_Decon = count_filenum(myParams.outputDir+ "decon/", Decon_File_Vec,"Decon.csv"); 
        std::vector<std::list<COMPOUND> > compoundsVec(nFile_Decon);
	/*each element of the vector*/
	//std::vector<std::list<COMPOUND> > compoundsVec(nFile);
        int j=0;
	 
	for(unsigned i=0;i<nFile;i++)  // before for(unsigned i=0;i<nFile;i++)
	{

		std::string DataFile=myParams.DataFiles.at(i);
		std::string filename;
		parseFileName(DataFile,filename);
                //Add the following code to use the real decon csv file
                if(!Find_InVector(Decon_File_Vec,filename))
                  continue;
                //add over
		std::vector<FEATURE_POINT> featurelist;
		//load deconvolution results from all files into single featurelist
		if(LoadFeatureListDecon(featurelist,myParams.outputDir,filename))
                {
		   /*convert the deconvoluted peaklist into compound (component) structure*/
		   getCompoundList(featurelist,&(compoundsVec.at(j)),j);
                   j++;
                }
	}
       
        if(j!=nFile_Decon) 
           std::cout << "Load feature number must eaqul to Decon File Number " << std::endl;

	time_t time1;
	time(&time1);
	std::cout<<"start alignment..."<<std::endl;;
	std::vector<ALIGNED_COMPOUNDS> alignedCompList;
	double scanInterval=scanIntervalVec.at(0);
	alignment(alignedCompList,compoundsVec,myParams.MaxCompoundETShift,myParams.MinSpecSimilarity1,myParams.MinSpecSimilarity2,myParams.TwoPhase,scanInterval);
	/*save the aligned compound into csv*/
	saveAlignedFeatureList(alignedCompList,myParams.outputDir,myParams.JobName);
	
	time_t time2;
	time(&time2);
	
	std::cout<<"alignment is done!"<<std::endl;
	std::cout<<difftime(time2, time1)<<std::endl;
}

void preProcess(const std::string & rawFile,PARAMS myParams,double &delaytime,double &scanInterval)
{
	std::cout<<rawFile<<std::endl;
	std::string filename;
	parseFileName (rawFile,filename);
	EICtype  EIC;
	unsigned nTotalScans;
	if(myParams.isSplitEIC)
	{
		nTotalScans=splitEIC(rawFile,EIC,myParams.outputDir,delaytime,scanInterval);//get eci
	}

	if(myParams.isPeakpicking | myParams.isDeconvolution)
	{
		loadEIC1(EIC,nTotalScans,myParams.outputDir,filename);
	}

	if(myParams.isPeakpicking)
	{
		peakpicking(rawFile,EIC,nTotalScans,myParams);//peak picking
	}
	if(myParams.isDeconvolution)
	{
		decon(rawFile,EIC,myParams);//decon
	}
	dispose_EIC(EIC);
}
/*append the dir to each raw file name to make thme as absolute paths*/
void addfullpath(std::vector<std::string> &DataFiles,std::string dir)
{
	for(std::vector<std::string>::iterator it=DataFiles.begin();it!=DataFiles.end();it++)
		it->insert(0,dir);
}
//fill parameters from object into myParams struct
void getPara(sqlite3_stmt *ppStmt,PARAMS &myParams)
{
	const char* pname=reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,0));
	if (strcmp(pname,"JobName")==0) 		
			myParams.JobName=reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,1));
	
	else if(strcmp(pname,"Owner")==0)
			myParams.Owner=reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,1));
	
	else if(strcmp(pname,"JobID")==0)
		myParams.JobID=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"WorkDir")==0)
		myParams.WorkDir=reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,1));
	
	else if(strcmp(pname,"DataDir")==0)
		myParams.DataDir=reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,1));
	
	else if(strcmp(pname,"File")==0) 
		myParams.DataFiles.push_back(reinterpret_cast<const char*>(sqlite3_column_text(ppStmt,1)));
	
	else if(strcmp(pname,"LocalWindowSize")==0)
		myParams.LocalWindowSize=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"EICPkRatio")==0)
		myParams.EICPkRatio=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"MaxScanShift")==0)
		myParams.MaxScanShift=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"FragPkRatio")==0)
		myParams.FragPkRatio=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"Fragments")==0)
		myParams.Fragments=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"MaxIntraDist")==0)
		myParams.MaxIntraDist=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"TwoPhase")==0)
		myParams.TwoPhase=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"isSmoothing")==0)
			myParams.isSmoothing=sqlite3_column_int(ppStmt,1);

	else if(strcmp(pname,"isBaseline")==0)
				myParams.isBaseline=sqlite3_column_int(ppStmt,1);

	else if(strcmp(pname,"isSplitEIC")==0)
					myParams.isSplitEIC=sqlite3_column_int(ppStmt,1);

	else if(strcmp(pname,"isPeakpicking")==0)
				myParams.isPeakpicking=sqlite3_column_int(ppStmt,1);

	else if(strcmp(pname,"isDeconvolution")==0)
					myParams.isDeconvolution=sqlite3_column_int(ppStmt,1);

	else if(strcmp(pname,"isAlignment")==0)
						myParams.isAlignment=sqlite3_column_int(ppStmt,1);
	//else if(strcmp(pname,"SampleRate")==0)
		//myParams.SampleRate=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"MaxCompoundETShift")==0)
		myParams.MaxCompoundETShift=sqlite3_column_int(ppStmt,1);
	
	else if(strcmp(pname,"MinSpecSimilarity1")==0)
	{
		myParams.MinSpecSimilarity1=sqlite3_column_double(ppStmt,1);
		myParams.MinSpecSimilarity2=myParams.MinSpecSimilarity1;
	}	

	else if(strcmp(pname,"nNode")==0)
			myParams.nNode=sqlite3_column_int(ppStmt,1);
			
}
//read parameters from sqlite database
void readPara_db(const char *dbpath,PARAMS & myParams)
{
	int rc;
	sqlite3* db;
	sqlite3_stmt *ppStmt;
	//char* db_err;

	rc=sqlite3_open(dbpath, &db);
	if( rc ){
		fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
		sqlite3_close(db);
		exit(1);
	}
	rc=sqlite3_prepare(db, "select * from options", -1,&ppStmt,0);
	if( rc!=SQLITE_OK ){
		fprintf(stderr, "Could not prepare statement\n");
		exit(1);
	}
	while(sqlite3_step(ppStmt) == SQLITE_ROW)
	{
		//fill parameters from object into myParams struct
		getPara(ppStmt,myParams);
	}
	
	sqlite3_finalize(ppStmt);
	sqlite3_close(db);
}

//paralell preprocessing stages
void ParPipleline(int argc, char **argv,PARAMS &myParams)
{
	time_t time1;
	time(&time1);
	


	double scanInterval;
	double * scanIntervalArray=NULL;//the vector to collect delaytime from each computing node
	double delaytime;//the local variable of each computing node for storing Retention time of the first scan, usually signal starts from 5 minutes
	double * delaytimeArray=NULL;//the vector to collect delaytime from each computing node


	unsigned nSample=myParams.DataFiles.size();
	int rank, nNode;//rank is the node id; nNode is the total number of nodes
	MPI_Status status;

	MPI_Init(&argc, &argv);//init and start the MPI, from now on,all the vairables above are broadcasted to each computing node 
	MPI_Comm_size(MPI_COMM_WORLD, &nNode);//get the total number of nodes
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);//get the current node ID
	unsigned nWoker=nNode-1;//total number of computing node, root node doesn't involve computing
	unsigned root=0;//specify root node ID
	
	if(rank==root)//code executed at root node
	{
		if(myParams.isSplitEIC)
		{
			delaytimeArray=new double[nSample];//claim storage space for collecting delaytime from each
			scanIntervalArray=new double[nSample];//claim storage space for collecting delaytime from each
		}
		unsigned nProcessed=0;//counter of compeleted samples
		unsigned fileindex;// sample file ID
		time_t time2;
		time(&time2);
		std::cout<<"Total preprocess time:"<<difftime(time2, time1)<<std::endl;
		
		/*start the iteration as long as there is any unprocessed sample;
		 once the while loop is done, all the samples are processed
		 */
		while(nProcessed<nSample)
		{
			/*loop through all the worker nodes until all the samples are processed*/
			for(unsigned source=1;source<nNode&&nProcessed<nSample;source++)
			{
				/*waiting for message of delaytime and fileindex from the respective worker node(nodeid=source);
				 MPI_Recv function keep the root node waiting until the current worker node (nodeid=source) finish
				 the preprocessing and send back the message
				 */

				if(myParams.isSplitEIC)
				{
					MPI_Recv(&delaytime, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&scanInterval, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
				}
				MPI_Recv(&fileindex, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
				if(myParams.isSplitEIC)
				{
					/*save the message into array*/
					delaytimeArray[fileindex]=delaytime;
					scanIntervalArray[fileindex]=scanInterval;
				}
				nProcessed++;//increment of sample counter 
			}
		}
		
		if(myParams.isSplitEIC)
		{
			/*convert the array into vector*/
			std::vector<double> delaytimeVec(delaytimeArray,delaytimeArray+nSample);
			std::vector<double> scanIntervalVec(scanIntervalArray,scanIntervalArray+nSample);
			for(unsigned i=0;i<delaytimeVec.size();i++)
			{
				std::cout<<"file"<<i<<": "<<delaytimeVec.at(i)<<std::endl;
				std::cout<<"file"<<i<<": "<<scanIntervalVec.at(i)<<std::endl;
			}
			/*save the delaytime of each sample into csv file*/
			saveFirstScanRT(scanIntervalVec,delaytimeVec,myParams.DataFiles,myParams.outputDir,myParams.JobName);
			if(myParams.isAlignment)
			{
				/*perform alignment*/
				align(myParams,scanIntervalVec);
			}
		}
		else
		{
			std::vector<double> delaytimeVec;
			std::vector<double> scanIntervalVec;
			LoadFirstScanRT(scanIntervalVec,delaytimeVec,myParams.outputDir,myParams.JobName);
			if(myParams.isAlignment)
			{
				/*perform alignment*/
				align(myParams,scanIntervalVec);
			}
		}
		time_t time3;
		time(&time3);
		std::cout<<"Alignment time:"<<difftime(time3, time2)<<std::endl;
	}
	else//code executed at worker nodes;distribute different files to worker nodes according to the file ID and node ID 
	{	/*loop through fileid*/
		for(unsigned FileIndex=0;FileIndex<myParams.DataFiles.size();FileIndex++)
		{
			/*only process the file ID which satisfy the module equation, so that files are evenly and exlusively assigned to the respective node*/
			if((FileIndex%nWoker+1)==(unsigned)rank)
			{	/*perform EIC splitting, peak picking and deconvolution*/

				preProcess(myParams.DataFiles.at(FileIndex),myParams,delaytime,scanInterval);
				/*after the preprocessing, send the message of delaytime and current processed sample ID to root node;
				 once the root recieved the message,the MPI_Recv function is finished holding
				 */
				if(myParams.isSplitEIC)
				{
					MPI_Send(&delaytime,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
					MPI_Send(&scanInterval,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
				}
				MPI_Send(&FileIndex,1,MPI_INT,0,1,MPI_COMM_WORLD);

			}
		}
	}
		
	MPI_Finalize();	
	
}

//sequential processing
void pipleline(PARAMS &myParams)
{
	time_t time1;
	time(&time1);

	std::vector<double> delaytime(myParams.DataFiles.size());
	std::vector<double> scanInterval(myParams.DataFiles.size());

	for(unsigned i=0;i<myParams.DataFiles.size();i++)
	{
		preProcess(myParams.DataFiles.at(i),myParams,delaytime.at(i),scanInterval.at(i));//3 preprocessing steps:EIC splitting,peak picking,deconvolution
		
	}
	if(myParams.isSplitEIC)
	{
		saveFirstScanRT(scanInterval,delaytime,myParams.DataFiles,myParams.outputDir,myParams.JobName);
	}
	if(myParams.isAlignment)
	{
		align(myParams,scanInterval);
	}

	time_t time2;
	time(&time2);
	std::cout<<"Total preprocess time:"<<difftime(time2, time1)<<std::endl;
	

}
int main(int argc, char **argv)
{
	/*read parameters and source file information from sqlite database*/
	PARAMS myParams;
	readPara_db(argv[1],myParams);

	/*specify output subdir*/
	myParams.outputDir=myParams.WorkDir+"output/";

	/*append the dir to each raw file name to make thme as absolute paths*/
	addfullpath(myParams.DataFiles,myParams.DataDir);	
	/*check if the respective output subdir exists, if not,then create them*/
	checkOutputFolder(myParams.WorkDir);
	/*launch the parallel pipeline with all procedures*/
    if(myParams.nNode>=2)
    {
    	std::cout<<"Parallel process! "<<std::endl;
    	ParPipleline(argc,argv,myParams);
    }
    else
    {
	    std::cout<<"Serialization process! "<<std::endl;
    	pipleline(myParams);
    }

	return(0);
	
}

