/*
 *  kmediod.h
 *  Janus
 *
 *  Created by mike on 2/11/10.
 *  Copyright 2010 ucc. All rights reserved.
 *
 */



#ifndef KMEDOID_H_INCLUDED
#define KMEDOID_H_INCLUDED
#include "common.h"
#include <stdlib.h>
#include <float.h>
#include <limits.h>

struct silhouetteType {
	double si;
	double ai;//avg intra-cluster distance of each cluster
	double bi;//avg inter-cluster distance of each cluster
	unsigned ncount;//number of elements of each cluster
	
};
struct clustPointType {
	unsigned clustAssignment;
	double neigbourDist;
};
struct ClusterSolution
{
	
	std::map<unsigned,silhouetteType> clustInfo;//pair of <clustID,sillhouette>
	std::vector<clustPointType> point;

};

void getclustermedoids(int nclusters, int nelements, double** distance,int clusterid[], int centroids[], double errors[]);
void normalization(std::vector<int> &data1,std::vector<int> &data2);
double dotproduct2Angle(double dotproduct);
double dotProduct (const std::vector<double> &data1,const std::vector<double> &data2);
double NistDotProduct(const std::vector<double> &data1,const std::vector<double> &data2);

double distFeaturePoint(FEATURE_POINT& p1, FEATURE_POINT& p2, EICtype & EIC);
double** distancematrix (std::vector<FEATURE_POINT*> curWindow, EICtype & EIC);
void kmedoids (unsigned nclusters, unsigned nelements, double** distmatrix,unsigned npass, unsigned clusterid[], double* error, int* ifound);
void randomassign (unsigned nclusters, unsigned nelements, unsigned clusterid[]);
double uniform(void);
int binomial(int n, double p);
void silhouette(ClusterSolution& clusterResult,double **distmatrix);

#endif
