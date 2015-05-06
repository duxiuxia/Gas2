#include "readCDF.h"

void readCDF(const char * sourcefile,int& nTotalScans,int& nTotalPoints,int*& arrayCount,int*& arrayInt,int*& arrayMz,double &delaytime,double &scanInterval)
{
	double *  arrayRT;
	std::string FileName;

    printf("inside readCDF.h");
    parseFileName(sourcefile,FileName);
	printf("%s reading data...\n",FileName.c_str());
   // time_t time1=time(NULL);

    /*read file*/
    NcFile dataFile(sourcefile, NcFile::ReadOnly);
    if (!dataFile.is_valid())
    {
     printf("Couldn't open file!\n");
		exit(-1);
      }

    NcVar * varRT = dataFile.get_var("scan_acquisition_time");//read cdf file
    NcVar * varInt = dataFile.get_var("intensity_values");
    NcVar * varMz = dataFile.get_var("mass_values");
    NcVar *varCount = dataFile.get_var("point_count");//point count for each scan

    nTotalScans=dataFile.get_dim("scan_number")->size();//get total number of scans
    nTotalPoints=dataFile.get_dim("point_number")->size(); //get total number of points
    /*get arrays*/
    arrayRT =new double[nTotalScans];
    arrayCount=new int[nTotalScans];
    arrayInt=new int[nTotalPoints];
    arrayMz=new int[nTotalPoints];

    varRT->get(arrayRT,nTotalScans);
    varCount->get(arrayCount,nTotalScans);

    varMz->get(arrayMz,nTotalPoints);
    varInt->get(arrayInt,nTotalPoints);

	delaytime=arrayRT[0];
	scanInterval=arrayRT[1]-arrayRT[0];
	delete []arrayRT;

}
