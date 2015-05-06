/*
 *  kmediod.cpp
 *  Janus
 *
 *  Created by mike on 2/11/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */
#include "kmedoid.h"

#if !defined(max)
#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
/* ********************************************************************* */
/* Find the center */
void getclustermedoids(unsigned nclusters, unsigned nelements, double** distance,
					   unsigned clusterid[], unsigned centroids[], double errors[])
{ unsigned i, j, k;
	for (j = 0; j < nclusters; j++) errors[j] = DBL_MAX;
	for (i = 0; i < nelements; i++)
	{ double d = 0.0;
		j = clusterid[i];
		for (k = 0; k < nelements; k++)
		{ if (i==k || clusterid[k]!=j) continue;
			d += (i < k ? distance[k][i] : distance[i][k]);
			if (d > errors[j]) break;
		}
		if (d < errors[j])
		{ errors[j] = d;
			centroids[j] = i;
		}
	}
}
//calculate the distance between two EICs
double distFeaturePoint(FEATURE_POINT& p1, FEATURE_POINT& p2, EICtype & EIC)
{
	//union two RT range
	int lboundInd=min(p1.lboundInd,p2.lboundInd);
	int rboundInd=max(p1.rboundInd,p2.rboundInd);
	
	int *e1=EIC[p1.mz];
	int * e2=EIC[p2.mz];	
	std::vector<double> data1(e1+lboundInd,e1+rboundInd+1);
	std::vector<double> data2(e2+lboundInd,e2+rboundInd+1);
	return dotproduct2Angle(dotProduct(data1,data2));
}


double dotproduct2Angle(double dotproduct)
{
	return (dotproduct-1)>0?0:acos(dotproduct)*180/M_PI;
		
}
double dotProduct (const std::vector<double> &data1,const std::vector<double> &data2)
{ 
	
	double result = 0.;
	double prd=0,denom1=0,denom2=0;
	//normalization
	double max1=*max_element(data1.begin(), data1.end());
	double max2=*max_element(data2.begin(), data2.end());
	double maxInt=max(max1,max2);
	//dot product
	for (unsigned i = 0; i <data1.size(); i++)
	{ 
		double x=(double)(data1.at(i))/maxInt;
		double y=(double)(data2.at(i))/maxInt;		


		prd +=x*y;
		denom1 += x*x;
		denom2 += y*y;			
	}
	result=prd/(sqrt(denom1)*sqrt(denom2));

	return result;
}
double NistDotProduct(const std::vector<double> &data1,const std::vector<double> &data2)
{ 
	
	double result = 0.;
	double prd=0,denom1=0,denom2=0;
	//normalization
	double max1=*max_element(data1.begin(), data1.end());
	double max2=*max_element(data2.begin(), data2.end());
	double maxInt=max(max1,max2);
	
	//NIST dot product
	for (unsigned i = 0; i <data1.size(); i++)
	{ 
		double x=(double)(data1.at(i))/maxInt;
		double y=(double)(data2.at(i))/maxInt;		

		
		prd +=sqrt(x*y);
		denom1 += x;
		denom2 += y;			
	}
	result=pow(prd,2)/(denom1*denom2);
	
	

	return result;
}



double** distancematrix (std::vector<FEATURE_POINT*> curWindow, EICtype & EIC)
{ 
	
	/* First determine the size of the distance matrix */
	const unsigned n = curWindow.size();
	unsigned i,j;
	double** matrix;
	
	if (n < 2) return NULL;
	
	/* Set up the ragged array */
	matrix = new double*[n];
	matrix[0] = NULL;
	/* The zeroth row has zero columns. We allocate it anyway for convenience.*/
	for (i = 1; i < n; i++)
		matrix[i] = new double[i];

	/* Calculate the distances and save them in the ragged array */
	for (i = 1; i < n; i++)
		for (j = 0; j < i; j++)
		{
			
			matrix[i][j]=distFeaturePoint(*curWindow.at(i),*curWindow.at(j),EIC);
			//std::cout<<i<<" to "<<j<<std::endl;
		}
	

		return matrix;
}
void silhouette(ClusterSolution& clusterResult,double **distmatrix)
{

	unsigned nPoints=clusterResult.point.size();
	
	//calculate the total intraCluster-distance for each cluster
	for (unsigned pointInd = 0; pointInd<nPoints; pointInd++)
	{ 

		unsigned curCentID = clusterResult.point.at(pointInd).clustAssignment;//get centroid ID
		if(pointInd ==curCentID)//if the center itself,skip the dist calculation
		   continue;
		double d =(pointInd < curCentID ? distmatrix[curCentID][pointInd] : distmatrix[pointInd][curCentID]);//get the distance between current element and its center

		if(clusterResult.clustInfo.find(curCentID)==clusterResult.clustInfo.end())//if current silhoutee does not exist,then add it
		{
			silhouetteType curSilhoutte;
			curSilhoutte.ai=d;//add the first distn to the silhoutte
			curSilhoutte.bi=0;//init the bi
			curSilhoutte.ncount=2;//init the element counter with two:one center and the first non-center elements
			curSilhoutte.si=0;//init si
			clusterResult.clustInfo[curCentID]=curSilhoutte;// add the current silhoutte to the respective cluster
		}
		else//if current silhoutee already exists,accumulate the distance and counter
		{
			clusterResult.clustInfo[curCentID].ai+=d;//add current dist to the respective cluster
			clusterResult.clustInfo[curCentID].ncount++;//element counter of respective cluster increment 
		}
	}
	
	std::map<unsigned,silhouetteType> * clustInfo=&(clusterResult.clustInfo);//get the cluster information pointer

	if(clustInfo->size()==1)
		clustInfo->begin()->second.bi=clustInfo->begin()->second.ai;//if only one cluster, then set the inter-dist=intro-dist,so that si will be zero after later calculation
	else
	{
		//calcuate the shortest inter-cluster distance for each point
		for (unsigned pointInd = 0; pointInd<nPoints; pointInd++)
		{
			clustPointType * curPoint=&(clusterResult.point.at(pointInd));
			unsigned curCentID = curPoint->clustAssignment;//get centID for current element
			
			for(std::map<unsigned,silhouetteType>::iterator it=clustInfo->begin();it!=clustInfo->end();it++)//loop to find the nearest neibgour cluster
			{
				unsigned neighbourCentID=it->first;//get center ID
				if(neighbourCentID!=curCentID)//check if the current center is the neighbour or itself
				{
					double d=(pointInd < neighbourCentID ? distmatrix[neighbourCentID][pointInd] : distmatrix[pointInd][neighbourCentID]);
					curPoint->neigbourDist=min(curPoint->neigbourDist,d);
				}
				else//if it is the same cluster,then skip
					continue;
			}
		}
		
		//calcuate the total inter-cluster distance for each cluster
		for (unsigned pointInd = 0; pointInd<nPoints; pointInd++)
		{
			clustPointType * curPoint=&(clusterResult.point.at(pointInd));
			unsigned curCentID = curPoint->clustAssignment;//get centID for current element
			clusterResult.clustInfo[curCentID].bi+=curPoint->neigbourDist;//add current neighbour dist to the respective cluster
		}
	}
	//calculate the average ai,bi,si for each cluster
	for(std::map<unsigned,silhouetteType>::iterator it=clustInfo->begin();it!=clustInfo->end();it++)
	{
		silhouetteType * curSilhoutte=&(it->second);
		curSilhoutte->ai=curSilhoutte->ai/curSilhoutte->ncount;
		curSilhoutte->bi=curSilhoutte->bi/curSilhoutte->ncount;
		curSilhoutte->si=(curSilhoutte->bi-curSilhoutte->ai)/max(curSilhoutte->bi,curSilhoutte->ai);
	}
	
}
void kmedoids (unsigned nclusters, unsigned nelements, double** distmatrix,
			   unsigned npass, unsigned clusterid[], double* error, int* ifound)

{ unsigned i, j, icluster;
	unsigned* tclusterid;
	unsigned* saved;
	unsigned* centroids;
	double* errors;
	unsigned ipass = 0;
	
	if (nelements < nclusters)
	{ *ifound = 0;
		return;
	} /* More clusters asked for than elements available */
	
	*ifound = -1;
	
	/* We save the clustering solution periodically and check if it reappears */
	saved = new unsigned[nelements];
	if (saved==NULL) return;
	
	centroids = new unsigned[nclusters];
	if(!centroids)
	{ free(saved);
		return;
	}
	
	errors =new double[nclusters];
	if(!errors)
	{ free(saved);
		free(centroids);
		return;
	}
	
	/* Find out if the user specified an initial clustering */
	if (npass<=1) tclusterid = clusterid;
	else
	{ tclusterid = new unsigned[nelements];
		if(!tclusterid)
		{ free(saved);
			free(centroids);
			free(errors);
			return;
		}
	}
	
	*error = DBL_MAX;
	do /* Start the loop */
	{ double total = DBL_MAX;
		int counter = 0;
		int period = 10;
		
		if (npass!=0) randomassign (nclusters, nelements, tclusterid);
		while(1)
		{ double previous = total;
			total = 0.0;
			
			if (counter % period == 0) /* Save the current cluster assignments */
			{ for (i = 0; i < nelements; i++) saved[i] = tclusterid[i];
				if (period < INT_MAX / 2) period *= 2;
			}
			counter++;
			
			/* Find the center */
			getclustermedoids(nclusters, nelements, distmatrix, tclusterid,
							  centroids, errors);
			
			for (i = 0; i < nelements; i++)
			/* Find the closest cluster */
			{ double distance = DBL_MAX;
				for (icluster = 0; icluster < nclusters; icluster++)
				{ double tdistance;
					j = centroids[icluster];
					if (i==j)
					{ distance = 0.0;
						tclusterid[i] = icluster;
						break;
					}
					tdistance = (i > j) ? distmatrix[i][j] : distmatrix[j][i];
					if (tdistance < distance)
					{ distance = tdistance;
						tclusterid[i] = icluster;
					}
				}
				total += distance;
			}
			if (total>=previous) break;
			/* total>=previous is FALSE on some machines even if total and previous
			 * are bitwise identical. */
			for (i = 0; i < nelements; i++)
				if (saved[i]!=tclusterid[i]) break;
			if (i==nelements)
				break; /* Identical solution found; break out of this loop */
		}
		
		for (i = 0; i < nelements; i++)
		{ 
			if (clusterid[i]!=centroids[tclusterid[i]])
			{ 
				if (total < *error)//when total >*error, the clusterid is not updated by previous centroids.Thus centroids basically are not necessary cluster id
				{ *ifound = 1;								 //but final Cluster ids is the final centroids
					*error = total;
					/* Replace by the centroid in each cluster. */
					for (j = 0; j < nelements; j++)
						clusterid[j] = centroids[tclusterid[j]];
				}
				break;
			}
		}
		if (i==nelements) (*ifound)++; /* break statement not encountered */
	} while (++ipass < npass);
	
	/* Deallocate temporarily used space */
	if (npass > 1) delete [] tclusterid;
	
	delete [] saved;
	delete [] centroids;
	delete [] errors;
	
	return;
}
/* ************************************************************************ */

void randomassign (unsigned nclusters, unsigned nelements, unsigned clusterid[])
/*
 Purpose
 =======
 
 The randomassign routine performs an initial random clustering, needed for
 k-means or k-median clustering. Elements (genes or microarrays) are randomly
 assigned to clusters. The number of elements in each cluster is chosen
 randomly, making sure that each cluster will receive at least one element.
 
 
 Arguments
 =========
 
 nclusters  (input) int
 The number of clusters.
 
 nelements  (input) int
 The number of elements to be clustered (i.e., the number of genes or microarrays
 to be clustered).
 
 clusterid  (output) int[nelements]
 The cluster number to which an element was assigned.
 
 ============================================================================
 */
{ unsigned i, j;
	unsigned k = 0;
	double p;
	unsigned n = nelements-nclusters;
	/* Draw the number of elements in each cluster from a multinomial
	 * distribution, reserving ncluster elements to set independently
	 * in order to guarantee that none of the clusters are empty.
	 */
	for (i = 0; i < nclusters-1; i++)
	{ p = 1.0/(nclusters-i);
		j = binomial(n, p);
		n -= j;
		j += k+1; /* Assign at least one element to cluster i */
		for ( ; k < j; k++) clusterid[k] = i;
	}
	/* Assign the remaining elements to the last cluster */
	for ( ; k < nelements; k++) clusterid[k] = i;
	
	/* Create a random permutation of the cluster assignments */
	for (i = 0; i < nelements; i++)
	{ j = (unsigned) (i + (nelements-i)*uniform());
		k = clusterid[j];
		clusterid[j] = clusterid[i];
		clusterid[i] = k;
	}
	
	return;
}

/* ************************************************************************ */

int binomial(int n, double p)
/*
 Purpose
 =======
 
 This routine generates a random number between 0 and n inclusive, following
 the binomial distribution with probability p and n trials. The routine is
 based on the BTPE algorithm, described in:
 
 Voratas Kachitvichyanukul and Bruce W. Schmeiser:
 Binomial Random Variate Generation
 Communications of the ACM, Volume 31, Number 2, February 1988, pages 216-222.
 
 
 Arguments
 =========
 
 p          (input) double
 The probability of a single event. This probability should be less than or
 equal to 0.5.
 
 n          (input) int
 The number of trials.
 
 
 Return value
 ============
 
 An integer drawn from a binomial distribution with parameters (p, n).
 
 ============================================================================
 */
{ const double q = 1 - p;
	if (n*p < 30.0) /* Algorithm BINV */
	{ const double s = p/q;
		const double a = (n+1)*s;
		double r = exp(n*log(q)); /* pow() causes a crash on AIX */
		int x = 0;
		double u = uniform();
		while(1)
		{ if (u < r) return x;
			u-=r;
			x++;
			r *= (a/x)-s;
		}
	}
	else /* Algorithm BTPE */
	{ /* Step 0 */
		const double fm = n*p + p;
		const int m = (int) fm;
		const double p1 = floor(2.195*sqrt(n*p*q) -4.6*q) + 0.5;
		const double xm = m + 0.5;
		const double xl = xm - p1;
		const double xr = xm + p1;
		const double c = 0.134 + 20.5/(15.3+m);
		const double a = (fm-xl)/(fm-xl*p);
		const double b = (xr-fm)/(xr*q);
		const double lambdal = a*(1.0+0.5*a);
		const double lambdar = b*(1.0+0.5*b);
		const double p2 = p1*(1+2*c);
		const double p3 = p2 + c/lambdal;
		const double p4 = p3 + c/lambdar;
		while (1)
		{ /* Step 1 */
			int y;
			int k;
			double u = uniform();
			double v = uniform();
			u *= p4;
			if (u <= p1) return (int)(xm-p1*v+u);
			/* Step 2 */
			if (u > p2)
			{ /* Step 3 */
				if (u > p3)
				{ /* Step 4 */
					y = (int)(xr-log(v)/lambdar);
					if (y > n) continue;
					/* Go to step 5 */
					v = v*(u-p3)*lambdar;
				}
				else
				{ y = (int)(xl+log(v)/lambdal);
					if (y < 0) continue;
					/* Go to step 5 */
					v = v*(u-p2)*lambdal;
				}
			}
			else
			{ const double x = xl + (u-p1)/c;
				v = v*c + 1.0 - fabs(m-x+0.5)/p1;
				if (v > 1) continue;
				/* Go to step 5 */
				y = (int)x;
			}
			/* Step 5 */
			/* Step 5.0 */
			k = abs(y-m);
			if (k > 20 && k < 0.5*n*p*q-1.0)
			{ /* Step 5.2 */
				double rho = (k/(n*p*q))*((k*(k/3.0 + 0.625) + 0.1666666666666)/(n*p*q)+0.5);
				double t = -k*k/(2*n*p*q);
				double A = log(v);
				if (A < t-rho) return y;
				else if (A > t+rho) continue;
				else
				{ /* Step 5.3 */
					double x1 = y+1;
					double f1 = m+1;
					double z = n+1-m;
					double w = n-y+1;
					double x2 = x1*x1;
					double f2 = f1*f1;
					double z2 = z*z;
					double w2 = w*w;
					if (A > xm * log(f1/x1) + (n-m+0.5)*log(z/w)
						+ (y-m)*log(w*p/(x1*q))
						+ (13860.-(462.-(132.-(99.-140./f2)/f2)/f2)/f2)/f1/166320.
						+ (13860.-(462.-(132.-(99.-140./z2)/z2)/z2)/z2)/z/166320.
						+ (13860.-(462.-(132.-(99.-140./x2)/x2)/x2)/x2)/x1/166320.
						+ (13860.-(462.-(132.-(99.-140./w2)/w2)/w2)/w2)/w/166320.)
						continue;
					return y;
				}
			}
			else
			{ /* Step 5.1 */
				int i;
				const double s = p/q;
				const double aa = s*(n+1);
				double f = 1.0;
				for (i = m; i < y; f *= (aa/(++i)-s));
				for (i = y; i < m; f /= (aa/(++i)-s));
				if (v > f) continue;
				return y;
			}
		}
	}
	/* Never get here */
	return -1;
}

/* ************************************************************************ */
/* *********************************************************************  */

double uniform(void)
/*
 Purpose
 =======
 
 This routine returns a uniform random number between 0.0 and 1.0. Both 0.0
 and 1.0 are excluded. This random number generator is described in:
 
 Pierre l'Ecuyer
 Efficient and Portable Combined Random Number Generators
 Communications of the ACM, Volume 31, Number 6, June 1988, pages 742-749,774.
 
 The first time this routine is called, it initializes the random number
 generator using the current time. First, the current epoch time in seconds is
 used as a seed for the random number generator in the C library. The first two
 random numbers generated by this generator are used to initialize the random
 number generator implemented in this routine.
 
 
 Arguments
 =========
 
 None.
 
 
 Return value
 ============
 
 A double-precison number between 0.0 and 1.0.
 ============================================================================
 */
{ int z;
	static const int m1 = 2147483563;
	static const int m2 = 2147483399;
	const double scale = 1.0/m1;
	
	static int s1 = 0;
	static int s2 = 0;
	
	if (s1==0 || s2==0) /* initialize */
	{ unsigned int initseed = (unsigned int) time(0);
		srand(initseed);
		s1 = rand();
		s2 = rand();
	}
	
	do
	{ int k;
		k = s1/53668;
		s1 = 40014*(s1-k*53668)-k*12211;
		if (s1 < 0) s1+=m1;
		k = s2/52774;
		s2 = 40692*(s2-k*52774)-k*3791;
		if(s2 < 0) s2+=m2;
		z = s1-s2;
		if(z < 1) z+=(m1-1);
	} while (z==m1); /* To avoid returning 1.0 */
	
	return z*scale;
}

/* ************************************************************************ */

