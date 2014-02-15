#include <stdlib.h>
#include <stdio.h>
#include "graph.h"
#include "niiImage.h"
#include "Mat.h"
typedef unsigned char uchar;

/*------------------------------------------------------------------------------------------*/
template<typename T>
class myedge
{
public:
	node_id id;
	T weight;
};
/*------------------------------------------------------------------------------------------*/
// edge[6] are 6 edges for a node: (x,y,z) + (-1,0,0)/(1,0,0)/(0,-1,0)/(0,1,0)/(0,0,-1)/(0,0,1)
// we can define a shift array so that we can compute the node id of the opposite side of the edge
template<typename T>
class mynode
{
public:
	node_id id;
	myedge<T> edge[6];
};
/*------------------------------------------------------------------------------------------*/
template<typename imgtype>
class Segment3D{

public:
	Segment3D(int dims[], niiImage<imgtype> *img, Mat3D<bool> *gt=NULL);
	~Segment3D();
	void setStarShapeConstraint2D(Mat2D<int> &seeds2D);
	void setStarShapeConstarint3D(Mat2D<int> &seeds3D);
	//function for setting bin number
	void setBinNum(int num);
	void computeTWeight();
	void computeNWeight();
	void setMask(int box_fg[]); 

private:
	void computeAppearanceModel2D();
	void computeAppearanceModel3D();
	float computeImageVariance();
	void add_edge(node_id i,node_id j,captype cap, captype rev_cap); //used for smoothterm
	void add_tweight(node_id i, tcaptype cap_source, tcaptype cap_sink); // used for bkcost and fgcost
	node_id getPixelID(int pos[],int dims[]);
	int getBinID(float x);
	void addNLinkForSingleDirection(float variance, float lambda, float dist, int start[], int end[], int shift[]);
	void seeds2DFromCentroid();

private:
	Graph<float,float,float> *graph;
	niiImage<imgtype> *image;
	Mat3D<uchar> *mask;
	Mat3D<float> *smoothimg;
	Mat3D<float> *dataimg;
	//Mat3D<float> *bkcost;
	//Mat3D<float> *fgcost;	
	//Mat3D<mynode<float> > *smoothterm;
	Mat3D<bool> *labeling;
	Mat3D<bool> *gt;
	Mat3D<int> *binID;
	Mat2D<int> *seeds3D;
	Mat2D<int> *seeds2D;
	int spacing[3];
	int dims[3];
	int binNum;
	imgtype tmin,tmax;

};