#include <fstream>
#include <time.h>
#include "readEIC.h"
#include <netcdfcpp.h>
static const int NC_ERR = 2;

//void getEIC (EICtype* EIC,int nTotalScans,const double*  arrayRT, const int* arrayCount,const int* arrayInt,const int* arrayMz)
void getEIC (EICtype* EIC,int nTotalScans,const int* arrayCount,const int* arrayInt,const int* arrayMz)
{

    /*get EICs*/
    printf("splitting into EICs...\n");
    unsigned nPos=0;//abslute postion of mz/int array
    for(int curScanID=0;curScanID<nTotalScans;curScanID++)
    {
        int nSpecPoints=arrayCount[curScanID];//number of points of the current scan
        for(int j=0;j<nSpecPoints;j++)
        {
            unsigned curMz=arrayMz[nPos];//get mz
            int curInt=arrayInt[nPos];//get int
            if(EIC->find(curMz)==EIC->end())
            {
                int * curIntArray=new int[nTotalScans];
				memset(curIntArray,0,nTotalScans*sizeof(int));//init all elements with 0
                (*EIC)[curMz]=curIntArray;//if current mz does not exist,create one vector for EIC
            }
            (*EIC)[curMz][curScanID]=curInt;
            nPos++;

        }

    }

     printf("done! %i EICs...\n",(int)EIC->size());


}

void dispose_EIC(EICtype& EIC)
{
     EICtype::iterator it;

    for(it=EIC.begin();it!=EIC.end();++it)
    {
        delete []it->second;
    }
}
void printEIC(EICtype & EIC,int nTotalScans)
{
    EICtype::iterator it;
	for(it=EIC.begin();it!=EIC.end();++it)
    {
		
        char sMz[4];
        sprintf(sMz,"%i",it->first);
		
        printf("%s\n",sMz);
        for(int i=0;i<nTotalScans;i++)
        {
            printf("%i ",it->second[i]);
        }
    }
}
/*save all EICs into multiple int vectors:each mz has one int vector
 it is intuitive but slow down the EIC file processing*/
int saveEIC(EICtype & EIC,unsigned nTotalScans,std::string outputDir,std::string filename)
{
	
	filename.insert(0,"EIC/");
	filename.insert(0,outputDir);
	filename.append("EIC.cdf");
	
	std::cout << "start saving EIC file:" << filename << std::endl;

	NcFile dataFile(filename.c_str(), NcFile::Replace);
	if (!dataFile.is_valid())
	{
		std::cout << "Couldn't open file!\n";
		return NC_ERR;
	}
	
	int nDim=nTotalScans;
	NcDim* Dim = dataFile.add_dim("totalScans",nDim);
	
	EICtype::iterator it;
	for(it=EIC.begin();it!=EIC.end();++it)
    {
		unsigned curMz=it->first;
		char sMz[4];
		sprintf(sMz,"%i",curMz);
		NcVar *varInt = dataFile.add_var(sMz, ncInt, Dim);
		int * arrayInt=it->second;
		varInt->put(arrayInt,nDim);
	}
	std::cout << "*** SUCCESS writing cdf file!" << std::endl;
	return 0;
	//std::ofstream outputfile;
//	outputfile.open("EIC.dat", std::ios::out);
//    EICtype::iterator it;
//	for(it=EIC.begin();it!=EIC.end();++it)
//    {
//		unsigned curMz=it->first;
//		outputfile.write((char*)(&curMz),sizeof(typeof(curMz)));
//		outputfile<<std::endl;
//		int * arrayInt=it->second;				 
//		outputfile.write((char*)arrayInt,nTotalScans*sizeof(typeof(*arrayInt)));
//		outputfile<<std::endl;
//	}
//	outputfile.close();
}
/*add all EIC vectors into one long vector:intVec in order to speed up EIC file processing*/
int saveEIC1(EICtype & EIC,unsigned nTotalScans,std::string outputDir,std::string filename)
{
	
	filename.insert(0,"EIC/");
	filename.insert(0,outputDir);
	filename.append("EIC.cdf");
	
	std::cout << "start saving EIC file:" << filename << std::endl;
	
	NcFile dataFile(filename.c_str(), NcFile::Replace);
	if (!dataFile.is_valid())
	{
		std::cout << "Couldn't open file!\n";
		return NC_ERR;
	}
	
	EICtype::iterator it;
	int nMz=EIC.size();
	int * mzArray=new int[nMz];
	int nlength=nTotalScans*nMz;
	int * IntArray=new int[nlength];
	unsigned i=0;
	for(it=EIC.begin();it!=EIC.end();++it)
    {
		mzArray[i]=it->first;
		memcpy(IntArray+i*nTotalScans, it->second,nTotalScans*sizeof(int));
		i++;
	}
	
	NcDim* Dim1 = dataFile.add_dim("mz",nMz);
	NcVar *varMz = dataFile.add_var("mzVec", ncInt, Dim1);
	varMz->put(mzArray,nMz);

	NcDim* Dim2 = dataFile.add_dim("int",nlength);
	NcVar *varInt = dataFile.add_var("intVec", ncInt, Dim2);
	varInt->put(IntArray,nlength);
	
	//NcDim* Dim = dataFile.add_dim("totalScans",nDim);
	
	
	std::cout << "*** SUCCESS writing cdf file!" << std::endl;
	delete []mzArray;
	delete []IntArray;


	return 0;
}
/*load EIC data from old structure*/

int loadEIC(EICtype & EIC,unsigned & nTotalScans ,unsigned startMz,unsigned endMz,std::string outputDir,std::string filename)
{
	filename.insert(0,"EIC/");
	filename.insert(0,outputDir);
	filename.append("EIC.cdf");
	std::cout<<"loading EIC:"<<filename<<std::endl;;
	// time_t time1=time(NULL);
	
    /*read file*/
    NcFile dataFile(filename.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid())
    {
		printf("Couldn't open file!\n");
		return NC_ERR;
	}
	
	nTotalScans=dataFile.get_dim("totalScans")->size();//get total number of scans

	for(unsigned curMz=startMz;curMz<=endMz;curMz++)
    {
		char sMz[4];
		sprintf(sMz,"%i",curMz);
		NcVar * varInt = dataFile.get_var(sMz);
		int * curIntArray=new int[nTotalScans];
		varInt->get(curIntArray,nTotalScans);
		EIC[curMz]=curIntArray;
	}

   	
	//std::ifstream inputfile;
//	inputfile.open("EIC.dat", std::ios::in);
//    EICtype::iterator it;
//	for(int i=0;i<nEIC;i++)
//    {
//		int curMz;
//		inputfile.getline((char*)(&curMz), nTotalScans);
//		int * curIntArray=new int[nTotalScans];
//		memset(curIntArray,0,nTotalScans);//init all elements with 0
//		inputfile.getline((char*)curIntArray, nTotalScans);
//		EIC[curMz]=curIntArray;//if current mz does not exist,create one vector for EIC
//	}
//	inputfile.close();
		return 0;
}
/*load EIC from new structure*/
int loadEIC1(EICtype & EIC,unsigned & nTotalScans,std::string outputDir,std::string filename)
{
	filename.insert(0,"EIC/");
	filename.insert(0,outputDir);
	filename.append("EIC.cdf");
	std::cout<<"loading EIC:"<<filename<<std::endl;;
	// time_t time1=time(NULL);

    /*read file*/
    NcFile dataFile(filename.c_str(), NcFile::ReadOnly);
    if (!dataFile.is_valid())
    {
		printf("Couldn't open file!\n");
		return NC_ERR;
	}

    //get mz vector
    unsigned nMzVec=dataFile.get_dim("mz")->size();//get total number of mzs
    int *arrayMz=new int[nMzVec];
    NcVar * varMz=dataFile.get_var("mzVec");
    varMz->get(arrayMz,nMzVec);

    //get int vector for all mzs
    unsigned nIntVec=dataFile.get_dim("int")->size(); //get total number of points
    int *arrayInt=new int[nIntVec];
    NcVar * varInt=dataFile.get_var("intVec");
    varInt->get(arrayInt,nIntVec);

	nTotalScans=nIntVec/nMzVec;//get total scan number

	for(unsigned i=0;i<nMzVec;i++)
    {
		char sMz[4];
		unsigned curMz=arrayMz[i];
		sprintf(sMz,"%i",curMz);
		unsigned offset=i*nTotalScans;
		//unsigned endInd=startInd+nTotalScans-1;

		int * curIntArray=new int[nTotalScans];
		memcpy(curIntArray,arrayInt+offset,sizeof(int)*nTotalScans);
		EIC[curMz]=curIntArray;
	}

	delete [] arrayMz;
	delete [] arrayInt;
	return 0;
}

