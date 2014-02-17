#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.VC80.CRT' version='8.0.50727.762' processorArchitecture='X86' publicKeyToken='1fc8b3b9a1e18e3b' language='*'\"")
#include <stdio.h>
#include <vector>
#include <math.h>
#include <queue> 
#include "segment_3d.hpp"
#include "starshape.hpp"
#include "fileOperation.h"
#include "Matrix3D.hpp"
#include <time.h>
#include <iostream>
#include <numeric>      // std::accumulate
#include <string>
using namespace std;

vector<string> path(3);
string out_path="D:/workspace/project/result/result_0216/paraTuning/";
//vector<float> fmeasure,recall,precision;  //to store the stastics


char image_indice[88][30]; 
int bin_num; 
float bmin, bmax;
float spacing[3];
int dims[3];
int box_fg[6];
int tightBox[6];
float ***img=NULL;
uchar ***gt=NULL;
int index;

//used for parameters tuning
float ***bkgCost=NULL;
float ***fgCost=NULL;
float lambdaxy[15]={0.1,0.3,0.5,0.7,0.9,1,3,5,7,9,11,13,15,17,19};
float lambdaz[15]={0.1,0.3,0.5,0.7,0.9,1,3,5,7,9,11,13,15,17,19};
int BinNum[20]={10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200};

float ctBias[20]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
float radiusRatio[10]={0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1};

enum MASKVALUE
{
	BACKGROUND=0,FOREGROUND,HDBACKGROUND,HDFOREGROUND,UNKNOW
};

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const string currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://www.cplusplus.com/reference/clibrary/ctime/strftime/
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

int getBinID(float x,float binmin,float binmax,int binnum)
{
	bin_id id=(int)((binnum-1)*(x-binmin)/(binmax-binmin));
	//if(id==binnum) id=binnum-1;
	return id;
}



inline node_id getPixelID(int pos[3], int dims[3]) 
{
	int s=pos[0]; int r=pos[1]; int c=pos[2];
	int slices=dims[0]; int rows=dims[1]; int cols=dims[2];
	return(s*rows*cols+r*cols+c);

}

void computePixelID(int***pixelID,int box[6],int newDims[3])
{
	int pos[3]; node_id id;
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				pos[0]=s-box[0]; pos[1]=r-box[1]; pos[2]=c-box[2];
				id=getPixelID(pos,newDims);
				pixelID[s][r][c]=id;
			}

}

float computeImageVariance(float ***img, int box[6])
{
	float DEFAULT_VARIANCE = (float) 0.0001;
	float v = (float) 0.0;
	int total = 0;
    for(int s=box[0]+1;s<=box[3];s++)
		for(int r=box[1]+1;r<=box[4];r++)
			for(int c=box[2]+1;c<=box[5];c++)
			{
				v=v+abs(img[s][r][c]-img[s][r-1][c])+abs(img[s][r][c]-img[s][r][c-1])+abs(img[s][r][c]-img[s-1][r][c]);
				total=total+3;
			}
	if (v/total < DEFAULT_VARIANCE ) 
		return(DEFAULT_VARIANCE);
	else return( v/total);
}

float computeImageVarianceXY(float ***img, int box[6])
{
	float DEFAULT_VARIANCE = (float) 0.0001;
	float v = (float) 0.0;
	int total = 0;
    for(int s=box[0];s<=box[3];s++)
		for(int r=box[1]+1;r<=box[4];r++)
			for(int c=box[2]+1;c<=box[5];c++)
			{
				v=v+abs(img[s][r][c]-img[s][r-1][c])+abs(img[s][r][c]-img[s][r][c-1]);
				total=total+2;
			}
	if (v/total < DEFAULT_VARIANCE ) 
		return(DEFAULT_VARIANCE);
	else return( v/total);
}

float computeImageVarianceZ(float ***img, int box[6])
{
	float DEFAULT_VARIANCE = (float) 0.0001;
	float v = (float) 0.0;
	int total = 0;
    for(int s=box[0]+1;s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				v=v+abs(img[s][r][c]-img[s-1][r][c]);
				total=total+1;
			}
	if (v/total < DEFAULT_VARIANCE ) 
		return(DEFAULT_VARIANCE);
	else return( v/total);
}

inline float computeNLinkCost(float diff,float variance, float lambda)
{
	return lambda*exp(-(diff*diff)/(variance*variance));

}




template<typename T>
void setValueInsideCube(T ***matrix,int dims[3], int box[6],T value)
{
	if(matrix==NULL) { printf("error in function setVauleInsideCube, matrix is empty\n"); exit(0); }
	if(box[0]<0 || box[1]<0 || box[2]<0 ||box[3]>dims[0]-1 ||box[4]>dims[1]-1 || box[5]>dims[2]-1)
	{ printf("error in function segVauleInsideCube, cube exceeds the range of data\n"); exit(0); }
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
				matrix[s][r][c]=value;

}

//getCubeIntersection used for shifting the box, but the box should not go out of the user provided bounding box
void getCubeIntersection(int box1[6],int box2[6],int *intersection)
{
	intersection[0]=box1[0]>box2[0] ? box1[0]:box2[0];
	intersection[1]=box1[1]>box2[1] ? box1[1]:box2[1];
	intersection[2]=box1[2]>box2[2] ? box1[2]:box2[2];
	intersection[3]=box1[3]<box2[3] ? box1[3]:box2[3];
	intersection[4]=box1[4]<box2[4] ? box1[4]:box2[4];
	intersection[5]=box1[5]<box2[5] ? box1[5]:box2[5];
}

void setMaskShiftBoxWithFgHdCst(uchar ***mask,int dims[3],int box_fg[6],vector<vector<int> > &seeds, int num_layer)
{
	int box_bk[6]={0,0,0,dims[0]-1,dims[1]-1,dims[2]-1};
	if(mask[0][0][0]!=HDBACKGROUND) setValueInsideCube<uchar>(mask,dims,box_bk,HDBACKGROUND);
	setValueInsideCube<uchar>(mask,dims,box_fg,BACKGROUND);
	int num_seed=seeds.size();
	int ct[2]={dims[1]/2,dims[2]/2};
//	int box_shift[6],box_hdshift[6];
	for(int i=0;i<num_seed;i++)
	{
		int dr=seeds[i][1]-ct[0]; 
		int dc=seeds[i][2]-ct[1];
		int box_shift[6]={seeds[i][0],box_fg[1]+dr,box_fg[2]+dc,seeds[i][0],box_fg[4]+dr,box_fg[5]+dc};
		int box_inter[6];
		getCubeIntersection(box_fg,box_shift,box_inter);
	    setValueInsideCube<uchar>(mask,dims,box_inter,FOREGROUND);  //set inside to be forground
	}

	//set hard constraint for foreground
	int size_fg[3]={box_fg[3]-box_fg[0]+1,box_fg[4]-box_fg[1]+1,box_fg[5]-box_fg[2]+1};
	int halfsize_hdfg[3]={size_fg[0]/6,size_fg[1]/6,size_fg[2]/6};
	int center[3]={dims[0]/2,dims[1]/2,dims[2]/2};
	int box_hdfg[6]={center[0]-halfsize_hdfg[0],center[1]-halfsize_hdfg[1],center[2]-halfsize_hdfg[2],center[0]+halfsize_hdfg[0],center[1]+halfsize_hdfg[1],center[2]+halfsize_hdfg[2]};
//	printf("hdbox %d %d %d %d %d %d\n",box_hdfg[0],box_hdfg[1],box_hdfg[2],box_hdfg[3],box_hdfg[4],box_hdfg[5]);
	if(size_fg[0]==1) { box_hdfg[0]=box_fg[0]; box_hdfg[3]=box_fg[3];}
	for(int i=box_hdfg[0];i<=box_hdfg[3];i++)
	{
		int s=i-1;
		int dr=seeds[s][1]-center[1]; 
		int dc=seeds[s][2]-center[2];
		int box_shift[6]={seeds[s][0],box_hdfg[1]+dr,box_hdfg[2]+dc,seeds[s][0],box_hdfg[4]+dr,box_hdfg[5]+dc};
		int box_inter[6];
		getCubeIntersection(box_fg,box_shift,box_inter);
	   setValueInsideCube<uchar>(mask,dims,box_inter,HDFOREGROUND);	   
	   
	}

}

// setMaskForAppearance is used for set label of the mask to construct appearance model
void setMaskShiftBox(uchar ***mask, int dims[3], int box_fg[6],vector<vector<int> > &seeds)
{
	int box_bk[6]={0,0,0,dims[0]-1,dims[1]-1,dims[2]-1};
	if(mask[0][0][0]!=HDBACKGROUND) setValueInsideCube<uchar>(mask,dims,box_bk,HDBACKGROUND);
	setValueInsideCube<uchar>(mask,dims,box_fg,BACKGROUND);
	int num_seed=seeds.size();
	int ct[2]={dims[1]/2,dims[2]/2};
//	int box_shift[6],box_hdshift[6];
	for(int i=0;i<num_seed;i++)
	{
		int dr=seeds[i][1]-ct[0]; 
		int dc=seeds[i][2]-ct[1];
		int box_shift[6]={seeds[i][0],box_fg[1]+dr,box_fg[2]+dc,seeds[i][0],box_fg[4]+dr,box_fg[5]+dc};
		int box_inter[6];
		getCubeIntersection(box_fg,box_shift,box_inter);
	    setValueInsideCube<uchar>(mask,dims,box_inter,FOREGROUND);  //set inside to be forground
	}

}

void setMask(uchar ***mask, int dims[3], int box_fg[3],vector<vector<int> > &seedsSide,int num_layer)
{
	int box_bk[6]={0,0,0,dims[0]-1,dims[1]-1,dims[2]-1};
	if(mask[0][0][0]!=HDBACKGROUND) setValueInsideCube<uchar>(mask,dims,box_bk,HDBACKGROUND);
	setValueInsideCube<uchar>(mask,dims,box_fg,FOREGROUND);
}

void computeBinID(float ***img,int ***BinID,int dims[3],float binmin, float binmax,float binnum)
{
	for(int s=0;s<dims[0];s++)
		for(int r=0;r<dims[1];r++)
			for(int c=0;c<dims[2];c++)
			{
				BinID[s][r][c]=getBinID(img[s][r][c],binmin,binmax,binnum);

			}
}

void computeHist(float ***img,uchar ***mask,int *** BinID,int dims[3],int box[6],vector<float> &fgProb ,vector<float> &bkProb)
{
	 
	 int bin_size=bkProb.size();
	 vector<float> fgHist(bin_size,0);
	 vector<float> bkHist(bin_size,0);
	 int countFg=0, countBkg=0;
	 //get bin id for every pixel in the image
	 for(int s=box[0];s<=box[3];s++)
	 {
		 for(int r=box[1];r<=box[4];r++)
			 for(int c=box[2];c<=box[5];c++)
			 {
				 bin_id binID=BinID[s][r][c];
				 if(mask[s][r][c]==FOREGROUND || mask[s][r][c]==HDFOREGROUND){
					 fgHist[binID]+=1;
					 countFg++;
				 }
				 else {
					 bkHist[binID]+=1;
					 countBkg++;
				 }
			 }
	 }
	 //to avoid division by 0
	 if(countFg==0) countFg=1;  
	 if(countBkg==0) countBkg=1;
	 for(int i=0;i<bin_size;i++)
	 {
		 fgProb[i]=fgHist[i]/countFg;
		 bkProb[i]=bkHist[i]/countBkg;
	 }

}

void computeDataterm(int ***BinID, int*** pixelID,int box[6],float ***fgDataCost,float ***bkDataCost,vector<float> &fgProb, vector<float> &bkProb )
{
	 float SMALLCONST=1.0e-3;

	 for(int s=box[0];s<=box[3];s++)
		 for(int r=box[1];r<=box[4];r++)
			 for(int c=box[2];c<=box[5];c++)
			 {
				node_id id=pixelID[s][r][c];
				bin_id binID=BinID[s][r][c];
				float fgCost=-log(fgProb[binID]+SMALLCONST);
				float bkCost=-log(bkProb[binID]+SMALLCONST);
                fgDataCost[s][r][c]=fgCost;
				bkDataCost[s][r][c]=bkCost;
			}
	 
}

void constructDatatermFromMatrix(GraphType *g, int box[6],int ***pixelID,float ***fgDataCost, float ***bkDataCost)
{
	 for(int s=box[0];s<=box[3];s++)
		 for(int r=box[1];r<=box[4];r++)
			 for(int c=box[2];c<=box[5];c++)
			 {
                g->add_tweights(pixelID[s][r][c],bkDataCost[s][r][c],fgDataCost[s][r][c]);
			}

}

void constructDataterm(GraphType *g,uchar ***mask,int ***BinID, int*** pixelID,int box[6],float ***dataterm,vector<float> &fgProb, vector<float> &bkProb )
{
	 float SMALLCONST=1.0e-3;
	 if(dataterm) {
	 for(int s=box[0];s<=box[3];s++)
		 for(int r=box[1];r<=box[4];r++)
			 for(int c=box[2];c<=box[5];c++)
			 {
				node_id id=pixelID[s][r][c];
				bin_id binID=BinID[s][r][c];
				float fgCost=-log(fgProb[binID]+SMALLCONST);
				float bkCost=-log(bkProb[binID]+SMALLCONST);
				dataterm[s][r][c]=bkCost-fgCost;
                g->add_tweights(id,bkCost,fgCost);
			}
	 }
	 else{
	 for(int s=box[0];s<=box[3];s++)
		 for(int r=box[1];r<=box[4];r++)
			 for(int c=box[2];c<=box[5];c++)
			 {
				node_id id=pixelID[s][r][c];
				bin_id binID=BinID[s][r][c];
				float fgCost=-log(fgProb[binID]+SMALLCONST);
				float bkCost=-log(bkProb[binID]+SMALLCONST);
                g->add_tweights(id,bkCost,fgCost);
			}
	 }
}



void addNLink(GraphType *g, float ***img,int ***pixelID,float ***smoothterm,float variance,float lambda, float dist,int box[6], int shift[3])
{
	//used for parameters tuning, remove after tuning

	//-----------------------------------------------
   node_id pixelID1, pixelID2;
   float diff,weight,weight1,weight2;
   if(smoothterm!=NULL)  // using the if we can use this function without saving smoothterm, this is used in parameter tuning
   {
   for(int s=box[0];s<=box[3];s++)
	   for(int r=box[1];r<=box[4];r++)
		   for(int c=box[2];c<=box[5];c++)
		   {
			   pixelID1=pixelID[s][r][c];
			   pixelID2=pixelID[s+shift[0]][r+shift[1]][c+shift[2]];
			   diff=img[s][r][c]-img[s+shift[0]][r+shift[1]][c+shift[2]];
			   weight=computeNLinkCost(diff,variance,lambda)/dist;
			   weight1=weight;
			   weight2=weight;
			   //if(diff>0) { weight1=weight; weight2=1;}
			   //else { weight1=1; weight2=weight; }
			   g->add_edge(pixelID1,pixelID2,weight1,weight2);
			   smoothterm[s][r][c]+=weight1;
			   smoothterm[s+shift[0]][r+shift[1]][c+shift[2]]+=weight2;
		   }
   }

   else{
   for(int s=box[0];s<=box[3];s++)
	   for(int r=box[1];r<=box[4];r++)
		   for(int c=box[2];c<=box[5];c++)
		   {
			   pixelID1=pixelID[s][r][c];
			   pixelID2=pixelID[s+shift[0]][r+shift[1]][c+shift[2]];
			   diff=img[s][r][c]-img[s+shift[0]][r+shift[1]][c+shift[2]];
			   weight=computeNLinkCost(diff,variance,lambda)/dist;
			   weight1=weight;
			   weight2=weight;
			   g->add_edge(pixelID1,pixelID2,weight1,weight2);
			   
		   }
   }

}

void computeNLink(float ***img,float ***smoothterm,float variance,float lambda, float dist,int box[6], int shift[3])
{
  float diff,weight;
  for(int s=box[0];s<=box[3];s++)
	   for(int r=box[1];r<=box[4];r++)
		   for(int c=box[2];c<=box[5];c++)
		   {
			   diff=img[s][r][c]-img[s+shift[0]][r+shift[1]][c+shift[2]];
			   weight=computeNLinkCost(diff,variance,lambda)/dist;
			   smoothterm[s][r][c]=weight;
		   }
}

void addNLinkFromMatrix(GraphType *g, int ***pixelID,float ***smoothterm, float lambda,int box[6], int shift[3])
{
	//used for parameters tuning, remove after tuning
	//float ***sm=NULL;
	//if(1==shift[0]) sm=Sx;
	//if(1==shift[1]) sm=Sy;
	//if(1==shift[2]) sm=Sz;
	//-----------------------------------------------
   node_id pixelID1, pixelID2;
   float weight;

   for(int s=box[0];s<=box[3];s++)
	   for(int r=box[1];r<=box[4];r++)
		   for(int c=box[2];c<=box[5];c++)
		   {
			   pixelID1=pixelID[s][r][c];
			   pixelID2=pixelID[s+shift[0]][r+shift[1]][c+shift[2]];
			   weight=lambda*smoothterm[s][r][c];
			   g->add_edge(pixelID1,pixelID2,weight,weight);
		   }

}

void constructSmoothtermFromMatrix(GraphType *g,int ***pixelID,float ***Sx, float ***Sy,float ***Sz,float lambda_xy,float lambda_z,int box[6],bool direction[3])
{
	int sbox[6],shift[3];   //sbox: range of pixel1 for adding smoothterm ( g->add_edge(pixel1,pixel2,weight1,weight2) )

	// n link for direction - on same slice
	if(direction[2] && lambda_xy!=0){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]; sbox[5]=box[5]-1;
	shift[0]=0; shift[1]=0; shift[2]=1;
	addNLinkFromMatrix(g,pixelID,Sx,lambda_xy,sbox,shift);
	}

	// n link for direction | on same slice
	if(direction[1] && lambda_xy!=0){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]-1; sbox[5]=box[5];
	shift[0]=0; shift[1]=1; shift[2]=0;
	addNLinkFromMatrix(g,pixelID,Sy,lambda_xy,sbox,shift);
	}

	// n link for direction | on different slices
	if(direction[0] && lambda_z!=0){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]-1; sbox[4]=box[4]; sbox[5]=box[5];
	shift[0]=1; shift[1]=0; shift[2]=0;
	addNLinkFromMatrix(g,pixelID,Sz,lambda_z,sbox,shift);
	}
}

void computeSmoothterm(float ***img,float ***Sx, float ***Sy, float ***Sz, float spacing[3],float variance, float lambda_xy,float lambda_z,int box[6],bool direction[3])
{
	int sbox[6],shift[3];   //sbox: range of pixel1 for adding smoothterm ( g->add_edge(pixel1,pixel2,weight1,weight2) )
	float dist;

	// n link for direction - on same slice
	if(direction[2]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]; sbox[5]=box[5]-1;
	shift[0]=0; shift[1]=0; shift[2]=1;
	dist=spacing[2];
	computeNLink(img,Sx,variance,lambda_xy,dist,sbox,shift);
	}

	// n link for direction | on same slice
	if(direction[1]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]-1; sbox[5]=box[5];
	shift[0]=0; shift[1]=1; shift[2]=0;
	dist=spacing[1];
	computeNLink(img,Sy,variance,lambda_xy,dist,sbox,shift);
	}

	// n link for direction | on different slices
	if(direction[0]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]-1; sbox[4]=box[4]; sbox[5]=box[5];
	shift[0]=1; shift[1]=0; shift[2]=0;
	dist=spacing[0];
	computeNLink(img,Sz,variance,lambda_z,dist,sbox,shift);
	}
}

void constructSmoothterm(GraphType *g,float ***img,int ***pixelID,float ***smoothterm,float spacing[3],float variance,float lambda_xy,float lambda_z,int box[6],bool direction[3])
{
	int sbox[6],shift[3];   //sbox: range of pixel1 for adding smoothterm ( g->add_edge(pixel1,pixel2,weight1,weight2) )
	float dist;

	// n link for direction - on same slice
	if(direction[2]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]; sbox[5]=box[5]-1;
	shift[0]=0; shift[1]=0; shift[2]=1;
	dist=spacing[2];
	addNLink(g,img,pixelID,smoothterm,variance,lambda_xy,dist,sbox,shift);
	}

	// n link for direction | on same slice
	if(direction[1]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]; sbox[4]=box[4]-1; sbox[5]=box[5];
	shift[0]=0; shift[1]=1; shift[2]=0;
	dist=spacing[1];
	addNLink(g,img,pixelID,smoothterm,variance,lambda_xy,dist,sbox,shift);
	}

	// n link for direction | on different slices
	if(direction[0]){
	sbox[0]=box[0]; sbox[1]=box[1]; sbox[2]=box[2];
	sbox[3]=box[3]-1; sbox[4]=box[4]; sbox[5]=box[5];
	shift[0]=1; shift[1]=0; shift[2]=0;
	dist=spacing[0];
	addNLink(g,img,pixelID,smoothterm,variance,lambda_xy,dist,sbox,shift);
	}
}

//get the labeling and return the volume of the labeling
int getLabeling(GraphType *g, int ***pixelID,uchar ***labeling,int box[6])
{
	int volume=0;
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				node_id id=pixelID[s][r][c];
				if(g->what_segment(id)== GraphType::SOURCE)
					{
						labeling[s][r][c]=1;
						volume++;
				}
				else 
					labeling[s][r][c]=0;
			}


	return volume;

}


void hardConstraint(GraphType *g, int ***pixelID, int box[6],uchar ***mask)
{
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				node_id id=pixelID[s][r][c];
				if(mask[s][r][c]==HDFOREGROUND) g->add_tweights(id,HDHUGE,0);

				else if(mask[s][r][c]==HDBACKGROUND) g->add_tweights(id,0,HDHUGE);

			}
}

// add 3d starshape constraint to graph, here there's just one seed
void setStarShapeConstraint3D(GraphType *g,int ***pixelID,int box[6],vector<vector<int> > &seeds3D)
{   int dims_star3D[3]={box[3]-box[0]+1,box[4]-box[1]+1,box[5]-box[2]+1};
	int num_pixels=dims_star3D[0]*dims_star3D[1]*dims_star3D[2];

	vector<int> idArr(num_pixels);
	int arraypos=0;
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
		
				idArr[arraypos]=pixelID[s][r][c];  // shift is node id shift for different layer		
				arraypos++;
			}

	int seed_star3D[3];
	int num_seeds=seeds3D.size();
	for(int i=0;i<num_seeds;i++)
	{
		for(int j=0;j<3;j++)
		{
		seed_star3D[j]=seeds3D[i][j]-box_fg[j];
		}
        setStarTlinksAllPixels3D(g,dims_star3D,seed_star3D,idArr);
	}
}

// add 2d star shape constraint to some slices
void setStarShapeConstraint2D(GraphType *g,int ***pixelID,int box[6],vector<vector<int> > &seeds)
{
	int num_seeds=seeds.size();
	int dims_star2D[2]={box[4]-box[1]+1,box[5]-box[2]+1}; //dimension: (r,c)
	int num_pixels=dims_star2D[0]*dims_star2D[1];
	vector<int> idArr(num_pixels);
	for(int i=0;i<num_seeds;i++)
	{
		int s=seeds[i][0];
		if(s<box[0] || s>box[3]) continue;
	//	if(s==0) continue;   //changed 1014
		int arraypos=0;
		int seed_star2D[2]={seeds[i][1]-box[2],seeds[i][2]-box[2]};
	    for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				idArr[arraypos]=pixelID[s][r][c];
				arraypos++;
			}
	//	printf("seed: %d %d %d\n",pos_seed[0],pos_seed[1],pos_seed[2]);
		setStarTlinksAllPixels2D(g,idArr,dims_star2D,seed_star2D);
	}
}


void setMask(uchar ***mask, uchar ***labeling,int box[6])
{
    for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				if(mask[s][r][c]==HDFOREGROUND || mask[s][r][c]==HDBACKGROUND ) continue;
				mask[s][r][c]=labeling[s][r][c];
			}
}


void computeCentroid(uchar ***labeling, int dims[3], int box_fg[3], vector<vector<int> > &seeds)
{
	int num_seeds=box_fg[3]-box_fg[0]+1;
	for(int i=0;i<num_seeds;i++)
	{
		int s=box_fg[0]+i;
	    float sum[2]={0,0};
	    int num_fgpixels=0;
		seeds[i][0]=s;
		for(int r=box_fg[1];r<=box_fg[4];r++)
			for(int c=box_fg[2];c<=box_fg[5];c++)
			{
				if(labeling[s][r][c]==1)
				{
					num_fgpixels++;
					sum[0]+=(float)r;
					sum[1]+=(float)c;
				}
			}
	    // if there's no tumor classified on this slice, set the seed to the center of the coordinate
			if(num_fgpixels==0) {seeds[i][1]=dims[1]/2; seeds[i][2]=dims[2]/2; continue;}  // skip
		else {
		seeds[i][1]=(int)(sum[0]/num_fgpixels);
		seeds[i][2]=(int)(sum[1]/num_fgpixels);
		}
	}
}



void initialSeeds2D(uchar ***gt,vector<vector<int> > &seeds2D,int seed3D[3],int dims[3],int box_fg[6])
{
	int num_seeds=seeds2D.size();
	if(num_seeds==1) {seeds2D[0][0]=dims[0]/2; seeds2D[0][1]=dims[1]/2; seeds2D[0][2]=dims[2]/2; return;}
	computeCentroid(gt,dims,box_fg,seeds2D);
	int radiusR=4;
    int radiusC=4;
	int dr,dc,r,c,s;
	int t=0;
	//do{
 //   if(t>50) break;
	//dr=rand()%radiusR;
	//dc=rand()%radiusC;
	//r=seeds2D[0][0]+dr; c=seeds2D[0][1]+dc; s=seeds2D[0][2];
	//t++;
	//}
	//while(gt[r][c][s]==0 || gt[r-1][c][s]==0 || gt[r+1][c][s]==0 || gt[r][c-1][s]==0 || gt[r][c+1][s]==0);
	//if(t<=50) 
	//seeds2D[0][0]=r; seeds2D[0][1]=c; seeds2D[0][2]=s; 
	//t=0;
	//do{
	//if(t>50) break;
	//dr=2+rand()%radiusR;
	//dc=2+rand()%radiusC;
	//r=seeds2D[num_seeds-1][0]+dr; c=seeds2D[num_seeds-1][1]+dc; s=seeds2D[num_seeds-1][2];
	//t++;
	//}
	//while(gt[r][c][s]==0 || gt[r-1][c][s]==0 || gt[r+1][c][s]==0 || gt[r][c-1][s]==0 || gt[r][c+1][s]==0 );
	//if(t<=50) 
	//seeds2D[num_seeds-1][0]=r; seeds2D[num_seeds-1][1]=c; seeds2D[num_seeds-1][2]=s; 

	// seeds for 3d
	int cti=dims[0]/2-seeds2D[0][0]; 
	float kr,kc; int ds;
	seeds2D[cti][0]=dims[0]/2; seeds2D[cti][1]=dims[1]/2; seeds2D[cti][2]=dims[2]/2;
	//seeds for left slices
	kc=(float)(seeds2D[0][2]-seed3D[2])/(seed3D[0]-seeds2D[0][0]);
	kr=(float)(seeds2D[0][1]-seed3D[1])/(seed3D[0]-seeds2D[0][0]);
	ds=1;
	for(int i=seed3D[0]-seeds2D[0][0]-1;i>0;i--)
	{
		//ds=seed3D[2]-seeds2D[0][2]-i;
		r=seed3D[1]+floor(ds*kr);
		c=seed3D[2]+floor(ds*kc);
		s=seed3D[0]-ds;
		seeds2D[i][0]=s; seeds2D[i][1]=r; seeds2D[i][2]=c;
		ds++;
	}
    //seeds for right slices
	kc=(float)(seeds2D[num_seeds-1][2]-seed3D[2])/(seeds2D[num_seeds-1][0]-seed3D[0]);
	kr=(float)(seeds2D[num_seeds-1][1]-seed3D[1])/(seeds2D[num_seeds-1][0]-seed3D[0]);
	ds=1;
	for(int i=seed3D[0]-seeds2D[0][0]+1;i<num_seeds-1;i++)
	{
		//ds=seeds2D[num_seeds-1][2]-seed3D[2]-i;
		r=seed3D[1]+floor(ds*kr);
		c=seed3D[2]+floor(ds*kc);
		s=seed3D[0]+ds;
		seeds2D[i][0]=s; seeds2D[i][1]=r; seeds2D[i][2]=c;
		ds++;
	}
 //   int num_seeds=seeds2D.size();
	//for(int i=1;i<num_seeds/2;i++)
	//	seeds2D[i]=seeds2D[0];
	//for(int i=num_seeds-2;i>num_seeds/2;i--)
	//	seeds2D[i]=seeds2D[num_seeds-1];
}

float computeRecall(uchar ***gt, uchar ***labeling, int box_fg[6])
{
	float tp=0;
	float relevant=0.0;
	for(int s=box_fg[0];s<=box_fg[3];s++)
	{
		for(int r=box_fg[1];r<=box_fg[4];r++)
			for(int c=box_fg[2];c<=box_fg[5];c++)
			{
				if(gt[s][r][c]!=0) {
				relevant+=1;
				if(labeling[s][r][c]==1) tp+=1;
				}
			}
	}
	if(relevant==0) return 0;
	return tp/relevant;
}

float computePrecision(uchar ***gt, uchar ***labeling, int box_fg[6])
{
	float tp=0, retrieved=0;
	for(int s=box_fg[0];s<=box_fg[3];s++)
		for(int r=box_fg[1];r<=box_fg[4];r++)
			for(int c=box_fg[2];c<=box_fg[5];c++)
			{
				if(labeling[s][r][c]==1) {
				retrieved+=1;
				if(gt[s][r][c]!=0) tp+=1;
				}
			}
	if(retrieved==0) return 0;
	return tp/retrieved;
}


template<typename T>
void writeVector2D(char *file,vector<vector<T> > &src,char *format)
{
   int n=src.size();
   int m=src[0].size();
   FILE *fid=fopen(file,"w");
   fprintf(fid,"%d %d\n",n,m); //print the size of the 2D vector
   for(int i=0;i<n;i++)
   {
	   for(int j=0;j<m;j++)
		   fprintf(fid,format,src[i][j]);
	   fprintf(fid,"\n");
   }
   fclose(fid);

}

// used for parameters tuning
template<typename T>
void readVector2D(char *file,vector<vector<T> > &dst,char *format)
{
   int n,m;
   FILE *fid=fopen(file,"r");
   if(fid==NULL) {printf("file %s does not exist\n",file); exit(0); }
   fscanf(fid,"%d %d\n",&n,&m);
   dst=vector<vector<T> >(n,vector<T>(m,0));
   T v;
   for(int i=0;i<n;i++)
   {
	   for(int j=0;j<m;j++)
		  {
			  fscanf(fid,format,&v);
			  dst[i][j]=v;
	   }
	   
   }
   fclose(fid);

}
//-----------------------------------------------------------

template<typename T>
void writeVector1D(char *file,vector<T> &src,char *format)
{
	int n=src.size();
	FILE *fid=fopen(file,"w");
	for(int i=0;i<n;i++)
		fprintf(fid,format,src[i]);
	fclose(fid);

}

void writeStastistic(char *file,uchar ***gt, uchar ***labeling,int box_fg[3])
{
	FILE *fid=fopen(file,"a+");
	fprintf(fid,"%d th image: %s\n",index,image_indice[index]);
	float recall=computeRecall(gt,labeling,box_fg);  //recall
	fprintf(fid,"recall is: %.3f%\n",recall);
	float precision=computePrecision(gt,labeling,box_fg); //precision
	fprintf(fid,"precision is: %.3f%\n",precision);
	fprintf(fid,"fmeasure is: %.3f\n",2*recall*precision/(recall+precision));
	fprintf(fid,"\n");
	fclose(fid);
}

float computeFmeasure(uchar ***gt, uchar ***labeling,float &re, float &pre,int box[6])
{
	float tp=0;
	float relevant=0; float retrieved=0;
	for(int s=box[0];s<=box[3];s++)
	{
		for(int r=box[1];r<=box[4];r++)
		{
			for(int c=box[2];c<=box[5];c++)
			{
				if(gt[s][r][c]!=0) {
					relevant+=1;
					if(labeling[s][r][c]!=0) {
					        tp+=1;
					}
				}
				if(labeling[s][r][c]!=0) retrieved+=1;
			}
		}
	}
	if(tp==0)
	{ 
		printf("no pixel classified to foreground, forget to put hard constraint?\n");
		re=0; pre=0;
		return 0;
	}
	re=tp/relevant;
	pre=tp/retrieved;
	float fm=2*re*pre/(re+pre);
	return fm;
}

bool isInRange(int pos[3],int dims[3])
{
	if(pos[0]>=0 && pos[1]>=0 && pos[2]>=0 && pos[0]<dims[0] && pos[1]<dims[1] && pos[2]<dims[2])
		return true;
	return false;
}


void setMaskCircularSeeds(uchar ***mask,vector<vector<int> > &seeds,int radius,int val)
{
	int diameter=1+2*radius;
	vector<vector<int> > circleMask(diameter,vector<int>(diameter,0));
	int ct[2]={radius,radius};
	for(int i=0;i<diameter;i++)
		for(int j=0;j<diameter;j++)
		{
			if((i-ct[0])*(i-ct[0])+(j-ct[1])*(j-ct[1])<=radius*radius) circleMask[i][j]=1;
		}

	int num_seeds=seeds.size();
	for(int index_seeds=0;index_seeds<num_seeds;index_seeds++)
	{
		int cts=seeds[index_seeds][0]; int ctr=seeds[index_seeds][1]; int ctc=seeds[index_seeds][2];
		for(int i=-radius;i<=radius;i++)
			for(int j=-radius;j<=radius;j++)
			{
				if(circleMask[radius+i][radius+j]==1)
				{
					int r=ctr+i; int c=ctc+j; int s=cts;
					mask[s][r][c]=val;
				}
			}
	}
}



void paraTuningForLocalHist(float ***img,int dims[3],int segBox[6],float spacing[3],float variance,
	vector<vector<int> > &seeds3D, vector<vector<int> > &seeds2D, 
	int ***BinID,int ***pixelID,uchar ***mask,float ***smoothterm,float ***dataterm,uchar ***labeling)
{
	float ***Sx=create3DMat<float>(dims,0);
	float ***Sy=create3DMat<float>(dims,0);
	float ***Sz=create3DMat<float>(dims,0);
	float ***bkDataCost=create3DMat<float>(dims,0);
	float ***fgDataCost=create3DMat<float>(dims,0);

	//pixelID can be different because of different segmentation regions
	//pixelID is computed for every different region segmentation
	//    int appearanceBox[6]={segBox[0],0,0,segBox[3],dims[1]-1,dims[2]-1};   //box for appearance model
	int graphDims[3]={segBox[3]-segBox[0]+1,segBox[4]-segBox[1]+1,segBox[5]-segBox[2]+1};
	int num_nodes=graphDims[0]*graphDims[1]*graphDims[2];
	int num_edges=10*num_nodes;  //set larger number 
	computePixelID(pixelID,segBox,graphDims);  //compute pixelID for segBox
	GraphType *g=new GraphType(num_nodes,num_edges);
	bool direction[3]={true,true,true};
	if(graphDims[0]==1) direction[0]=false;

	//smoothterm cost for lambda_xy=1, lambda_z=1
	computeSmoothterm(img,Sx,Sy,Sz,spacing,variance,1,1,segBox,direction);

	float fmeasure[20][15][15],recall[20][15][15],precision[20][15][15];
	//begin of the tuning of bin num
	int binNum; float re,pre,fm;
	float curr_lambda_xy,curr_lambda_z, last_lambda_xy=0,last_lambda_z=0;

	for(int iter_bin=0;iter_bin<20;iter_bin++){
		binNum=BinNum[iter_bin];
		computeBinID(img,BinID,dims,bmin,bmax,binNum);
		vector<float> fgProb(binNum,0); //histogram for foreground
		vector<float> bkgProb(binNum,0); //histogram for background

		for(int i=segBox[0];i<=segBox[3];i++){
			int appearanceBox[6]={i,0,0,i,dims[1]-1,dims[2]-1};   //box for appearance model
			int localBox[6]={i,segBox[1],segBox[2],i,segBox[4],segBox[5]};
			computeHist(img,mask,BinID,dims,appearanceBox,fgProb,bkgProb);
			computeDataterm(BinID,pixelID,localBox,fgDataCost,bkDataCost,fgProb,bkgProb);
		}

		for(int iter_lz=0;iter_lz<15;iter_lz++){
			g->add_node(num_nodes);
			//unary term
			hardConstraint(g,pixelID,segBox,mask);
			constructDatatermFromMatrix(g,segBox,pixelID,fgDataCost,bkDataCost);

			//pair wise term
			//star shape constraint
			if(graphDims[0]!=1){
				setStarShapeConstraint3D(g,pixelID,segBox,seeds3D);
			}
			setStarShapeConstraint2D(g,pixelID,segBox,seeds2D);

			//smoothness term
			last_lambda_xy=0;
			curr_lambda_z=lambdaz[iter_lz];
			constructSmoothtermFromMatrix(g,pixelID,Sx,Sy,Sz,0,curr_lambda_z,segBox,direction);

			for(int iter_lxy=0;iter_lxy<15;iter_lxy++){
				curr_lambda_xy=lambdaxy[iter_lxy];
				float add_lambda_xy=curr_lambda_xy-last_lambda_xy;
				constructSmoothtermFromMatrix(g,pixelID,Sx,Sy,Sz,add_lambda_xy,0,segBox,direction);
				float flow=g->maxflow();
				//     printf("flow is %f\n",flow);
				getLabeling(g,pixelID,labeling,segBox);

				fm=computeFmeasure(gt,labeling,re,pre,segBox);
				fmeasure[iter_bin][iter_lxy][iter_lz]=fm; 
				recall[iter_bin][iter_lxy][iter_lz]=re; 
				precision[iter_bin][iter_lxy][iter_lz]=pre;

				last_lambda_xy=curr_lambda_xy;
			}
			g->reset();
		}

	}

	delete g;
	free3DMat<float>(Sx);
	free3DMat<float>(Sy);
	free3DMat<float>(Sz);
	free3DMat<float>(fgDataCost);
	free3DMat<float>(bkDataCost);

	string fmfile=out_path+"fmeasure/"+string(image_indice[index])+".txt";
	string refile=out_path+"recall/"+string(image_indice[index])+".txt";
	string prefile=out_path+"precision/"+string(image_indice[index])+".txt";
	int statDim[3]={20,15,15};

	write3DMat<float>(fmfile.c_str(),(float *)fmeasure,dims);
	write3DMat<float>(refile.c_str(),(float *)recall,dims);
	write3DMat<float>(prefile.c_str(),(float *)precision,dims);

}



void segmentGlobalHist(float ***img,int dims[3],int segBox[6],float spacing[3],float variance,int binNum,
	float lambda_xy, float lambda_z,int iter_times, vector<vector<int> > &seeds3D, vector<vector<int> > &seeds2D, 
	int ***BinID,int ***pixelID,uchar ***mask,float ***smoothterm,float ***dataterm,uchar ***labeling)
{
	    //pixelID can be different because of different segmentation regions
	    //pixelID is computed for every different region segmentation
	    int appearanceBox[6]={segBox[0],0,0,segBox[3],dims[1]-1,dims[2]-1};   //box for appearance model
	    int graphDims[3]={segBox[3]-segBox[0]+1,segBox[4]-segBox[1]+1,segBox[5]-segBox[2]+1};
		int num_nodes=graphDims[0]*graphDims[1]*graphDims[2];
	    int num_edges=10*num_nodes;  //set larger number 
		vector<float> fgProb(binNum,0); //histogram for foreground
		vector<float> bkgProb(binNum,0); //histogram for background
		computePixelID(pixelID,segBox,graphDims);  //compute pixelID for segBox
	    GraphType *g=new GraphType(num_nodes,num_edges);
	    bool direction[3]={true,true,true};
		if(graphDims[0]==1) direction[0]=false;

		float re,pre,fm;
		int iter_i=0;
		while(iter_i<iter_times){
	    g->add_node(num_nodes);
		hardConstraint(g,pixelID,segBox,mask);
		if(graphDims[0]!=1){
		  setStarShapeConstraint3D(g,pixelID,segBox,seeds3D);
		}
		setStarShapeConstraint2D(g,pixelID,segBox,seeds2D);
		computeHist(img,mask,BinID,dims,appearanceBox,fgProb,bkgProb);
		constructDataterm(g,mask,BinID,pixelID,segBox,dataterm,fgProb,bkgProb);
		constructSmoothterm(g,img,pixelID,smoothterm,spacing,variance,lambda_xy,lambda_z,segBox,direction);
		float flow=g->maxflow();
        printf("flow is %f\n",flow);
		getLabeling(g,pixelID,labeling,segBox);
		fm=computeFmeasure(gt,labeling,re,pre,segBox);
		printf("fm:%f re:%f pre:%f \n",fm,re,pre);
		g->reset();
		setMask(mask,labeling,segBox);
		writeData(path[iter_i],string(image_indice[index]),dims,segBox,labeling,dataterm,smoothterm,seeds2D,true);
		iter_i++;
		}

		delete g;

}

void segmentLocalHist(float ***img,int dims[3],int segBox[6],float spacing[3],float variance,int binNum,
	float lambda_xy, float lambda_z,int iter_times, vector<vector<int> > &seeds3D, vector<vector<int> > &seeds2D, 
	int ***BinID,int ***pixelID,uchar ***mask,float ***smoothterm,float ***dataterm,uchar ***labeling)
{
	    //pixelID can be different because of different segmentation regions
	    //pixelID is computed for every different region segmentation
	    
	    int graphDims[3]={segBox[3]-segBox[0]+1,segBox[4]-segBox[1]+1,segBox[5]-segBox[2]+1};
		int num_nodes=graphDims[0]*graphDims[1]*graphDims[2];
	    int num_edges=10*num_nodes;  //set larger number 
		vector<float> fgProb(binNum,0); //histogram for foreground
		vector<float> bkgProb(binNum,0); //histogram for background
		computePixelID(pixelID,segBox,graphDims);  //compute pixelID for segBox
	    GraphType *g=new GraphType(num_nodes,num_edges);
		bool direction[3]={true,true,true};
		if(graphDims[0]==1) direction[0]=false;

		int iter_i=0;
		while(iter_i<iter_times){
	    g->add_node(num_nodes);
		hardConstraint(g,pixelID,segBox,mask);
		if(graphDims[0]!=1){
		  setStarShapeConstraint3D(g,pixelID,segBox,seeds3D);
		}
		setStarShapeConstraint2D(g,pixelID,segBox,seeds2D);

		for(int i=segBox[0];i<=segBox[3];i++){
			int appearanceBox[6]={i,0,0,i,dims[1]-1,dims[2]-1};   //box for appearance model
			int localBox[6]={i,segBox[1],segBox[2],i,segBox[4],segBox[5]};
		    computeHist(img,mask,BinID,dims,appearanceBox,fgProb,bkgProb);
		    constructDataterm(g,mask,BinID,pixelID,localBox,dataterm,fgProb,bkgProb);
		}

		constructSmoothterm(g,img,pixelID,smoothterm,spacing,variance,lambda_xy,lambda_z,segBox,direction);
		float flow=g->maxflow();
        printf("flow is %f\n",flow);
		getLabeling(g,pixelID,labeling,segBox);
		g->reset();
		setMask(mask,labeling,segBox);
		writeData(path[iter_i],string(image_indice[index]),dims,segBox,labeling,dataterm,smoothterm,seeds2D,true);
		iter_i++;
		}

		delete g;

}

void markAllNodes(GraphType *g, int num_nodes)
{
	for(int i=0;i<num_nodes;i++) g->mark_node(i);
}

//data term and starshape constraint should have been added into the graph *g
//variance is the variance 
void adaptiveLambdaXY(GraphType *g,float ***img,int box[6], int targetVolume,float spacing[3],float stepLambdaXY,
      float maxLambdaXY,float stopLevel,int ***pixelID,uchar ***labeling, int &bestVolume, float &bestLambdaXY)
{	
	//the first iteration is the first graph cut
	int num_nodes=(box[3]-box[0]+1)*(box[4]-box[1]+1)*(box[5]-box[2]+1);
	float lambdaXY=stepLambdaXY;
	float variance=computeImageVarianceXY(img,box);
	bool direction[3]={false,true,true};
	constructSmoothterm(g,img,pixelID,NULL,spacing,variance,lambdaXY,0,box,direction);
    float flow=g->maxflow();
	int volume=getLabeling(g,pixelID,labeling,box);
	int minDiff=abs(volume-targetVolume); bestLambdaXY=lambdaXY; bestVolume=volume;
	int threshold=stopLevel*targetVolume;
//	printf("targetVolume %d, lambdaXY: %f, volume: %d\n",targetVolume,lambdaXY,volume);
	while(lambdaXY<=maxLambdaXY){
		lambdaXY+=stepLambdaXY;
		constructSmoothterm(g,img,pixelID,NULL,spacing,variance,stepLambdaXY,0,box,direction);
		markAllNodes(g,num_nodes);
		flow=g->maxflow(true);
		volume=getLabeling(g,pixelID,labeling,box);
//		printf("targetVolume %d, lambdaXY: %f, volume: %d\n",targetVolume,lambdaXY,volume);
		if(abs(volume-targetVolume)<minDiff) {
			minDiff=abs(volume-targetVolume);
			bestLambdaXY=lambdaXY;
			bestVolume=volume;
		}
	//	else break;
	}

}

float segmentAdaptiveLambdaXY(float ***img,int dims[3],int box[6],int tightBox[6],float spacing[3],float variance,int binNum,
	float lambda_xy, float lambda_z,int iter_times, vector<vector<int> > &seeds3D, vector<vector<int> > &seeds2D, 
	int ***BinID,int ***pixelID,uchar ***mask,float ***smoothterm,float ***dataterm,uchar ***labeling)
{
	    //pixelID can be different because of different segmentation regions
	    //pixelID is computed for every different region segmentation
	    
	    int appearanceBox[6]={box[0],0,0,box[3],dims[1]-1,dims[2]-1};   //box for appearance model
		
	    int graphDims[3]={box[3]-box[0]+1,box[4]-box[1]+1,box[5]-box[2]+1};
		int num_nodes=graphDims[0]*graphDims[1]*graphDims[2];
	    int num_edges=10*num_nodes;  //set larger number 
		vector<float> fgProb(binNum,0); //histogram for foreground
		vector<float> bkgProb(binNum,0); //histogram for background
		computePixelID(pixelID,box,graphDims);  //compute pixelID for segBox
	    GraphType *g=new GraphType(num_nodes,num_edges);
	    bool direction[3]={true,true,true};
		if(graphDims[0]==1) direction[0]=false;

	    g->add_node(num_nodes);
		hardConstraint(g,pixelID,box,mask);
		//setStarShapeConstraint3D(g,pixelID,segBox,seeds3D);
		setStarShapeConstraint2D(g,pixelID,box,seeds2D);
		computeHist(img,mask,BinID,dims,appearanceBox,fgProb,bkgProb);
		constructDataterm(g,mask,BinID,pixelID,box,dataterm,fgProb,bkgProb);
		int targetVolume=(tightBox[4]-tightBox[1]+1)*(tightBox[5]-tightBox[2]+1);
		int bestVolume=0; float bestLambdaXY=0;
		adaptiveLambdaXY(g,img,box,targetVolume,spacing,0.1,10,0,pixelID,labeling,bestVolume,bestLambdaXY);
		delete g;
		printf("targetVolume %d, bestLambdaXY: %f, bestVolume: %d\n",targetVolume,bestLambdaXY,bestVolume);
		return bestLambdaXY;
}


void addDataterm(GraphType *g, int dims[3], int ***pixelID)
{
	for(int r=0;r<dims[0];r++)
		for(int c=0;c<dims[1];c++)
			for(int s=0;s<dims[2];s++)
			{
				g->add_tweights(pixelID[r][c][s],bkgCost[r][c][s],fgCost[r][c][s]);
			}
}



//---------------------------------------------------------------------------

void segmentation(float ***img,uchar ***gt, int dims[3],int box_fg[6],float spacing[3],float lambda_xy,float lambda_z, int binNum)
{
	//first create some matrix for intermediate result
	//construct mask for hard constraint and appearance model
	uchar ***mask=create3DMat<uchar>(dims,HDBACKGROUND);
	//matrix for storing histogram id for each pixel
	int ***BinID=create3DMat<int>(dims,0);
	//matrix for storing pixel's node_id
	int ***pixelID=create3DMat<int>(dims,-1); //-1 represents for non-graph node
	//matrix for smoothness term
	float ***smoothterm=create3DMat<float>(dims,0);
	//matrix for data term
	float ***dataterm=create3DMat<float>(dims,0);
	//matrix for the labeling
	uchar ***labeling=create3DMat<uchar>(dims,0);

	bin_num=binNum;
	computeBinID(img,BinID,dims,bmin,bmax,binNum);  //get BinID
	//write3DMat<int>("binID",BinID,dims);

	int box[6]={0,0,0,dims[0]-1,dims[1]-1,dims[2]-1};
	float variance=computeImageVariance(img,box);

	// star shape constraint
	vector<vector<int> > seeds3D(1,vector<int>(3,0));
	seeds3D[0][0]=dims[0]/2; seeds3D[0][1]=dims[1]/2; seeds3D[0][2]=dims[2]/2;
	int num_seeds=box_fg[3]-box_fg[0]+1;
	vector<vector<int> > seeds2D(num_seeds,vector<int>(3,0));
	vector<vector<int> > seedsSide(2,vector<int>(3,0));
	int center[3]={dims[0]/2,dims[1]/2,dims[2]/2};
	initialSeeds2D(gt,seeds2D,center,dims,box_fg);
	seedsSide[0][0]=seeds2D[0][0]; seedsSide[0][1]=seeds2D[0][1]; seedsSide[0][2]=seeds2D[0][2]; 
	seedsSide[1][0]=seeds2D[num_seeds-1][0]; seedsSide[1][1]=seeds2D[num_seeds-1][1]; seedsSide[1][2]=seeds2D[num_seeds-1][2]; 

	int seedsize=2>sqrt(dims[0]*dims[1]*0.0041/3.1417) ? 2:sqrt(dims[0]*dims[1]*0.0041/3.1417);

	int segBox[6]={box_fg[0],box_fg[1]-1,box_fg[2]-1,box_fg[3],box_fg[4]+1,box_fg[5]+1};

	setMaskShiftBoxWithFgHdCst(mask,dims,box_fg,seeds2D,1); //no shift, no hd
	setMaskCircularSeeds(mask,seedsSide,seedsize,HDFOREGROUND);
//	segmentGlobalHist(img,dims,segBox,spacing,variance,binNum,lambda_xy,lambda_z,2,seeds3D,seeds2D,BinID,pixelID,mask,smoothterm,dataterm,labeling);
	
	
	//int adaptiveSegBox[6]={tightBox[0],segBox[1],segBox[2],tightBox[3],segBox[4],segBox[5]};
	//lambda_xy=segmentAdaptiveLambdaXY(img,dims,adaptiveSegBox,tightBox,spacing,variance,binNum,lambda_xy,lambda_z,2,seeds3D,seeds2D,BinID,pixelID,mask,smoothterm,dataterm,labeling);
	//segmentGlobalHist(img,dims,segBox,spacing,variance,binNum,lambda_xy,lambda_z,2,seeds3D,seeds2D,BinID,pixelID,mask,smoothterm,dataterm,labeling);
	//

//	paraTuningForBinNumGlobalHist(img,dims,segBox,spacing,variance,lambda_xy,lambda_z,seeds3D,seeds2D,BinID,pixelID,mask,smoothterm,dataterm,labeling);
	paraTuningForLocalHist(img,dims,segBox,spacing,variance,seeds3D,seeds2D,BinID,pixelID,mask,smoothterm,dataterm,labeling);
	
	//float re,pre,fm;
	//fm=computeFmeasure(gt,labeling,re,pre,box);
	//printf("recall: %f precision: %f fm: %f \n",re,pre,fm);
	//recall.push_back(re); precision.push_back(pre); fmeasure.push_back(fm);

	free3DMat<uchar>(mask);
	free3DMat<int>(BinID);
	free3DMat<int>(pixelID);
	free3DMat<float>(smoothterm);
	free3DMat<float>(dataterm);
}




void paraTuningComputeDataCost(vector<float> &fgProb, vector<float> &bkProb, int ***BinID,float ***fgCost, float ***bkgCost, int box[6])
{
	for( int r=box[0];r<=box[3];r++)
		for(int c=box[1];c<=box[4];c++)
			for(int s=box[2];s<=box[5];s++)
			{
				fgCost[r][c][s]=-log(fgProb[BinID[r][c][s]]);
				bkgCost[r][c][s]=-log(bkProb[BinID[r][c][s]]);

			}

}

void paraTuningAddDataTerm(GraphType *g, float ***fgCost, float ***bkgCost,int dims[3], int box[6])
{
	node_id id; int pos[3];
	for( int r=box[0];r<=box[3];r++)
		for(int c=box[1];c<=box[4];c++)
			for(int s=box[2];s<=box[5];s++)
			{
				pos[0]=r; pos[1]=c; pos[2]=s;
				id=getPixelID(pos,dims);
				g->add_tweights(id,bkgCost[r][c][s],fgCost[r][c][s]);
			}

}




int main(int argc, char **argv)
{

	//int st,ed;
   int st=atoi(argv[1]); int ed=atoi(argv[2]);
   readImgIndice("D:/workspace/project/crop3D/recrop/TOMOB2INT_CROP/image_name.txt");
   

   int img_num=0;
   float lambdaBias=1; float lambda_xy=5, lambda_z=7;

   //string out_file(argv[3]);
   //string outFM=out_path+out_file+"_fm.txt";
   //string outRE=out_path+out_file+"_re.txt";
   //string outPE=out_path+out_file+"_pe.txt";

   //strcpy(out_path,argv[3]);
   char fmpath[300]="D:/workspace/project/result/result_0115/single_global_no_hd_no_bias/fmeasure/";
   char repath[300]="D:/workspace/project/result/result_1201/test_global_lambda_5_7_volume_04_19/recall/";
   char prepath[300]="D:/workspace/project/result/result_1201/test_global_lambda_5_7_volume_04_19/precision/";
 //  char out_path[300]="D:/workspace/project/result/result_1106/single_layer_local_hist_linear/";
   path[0]=out_path+"iteration1/";
   path[1]=out_path+"iteration2/";

//   st=0; ed=10;
   for(int i=st;i<ed;i++)
   {
	index=i;
	img_num++;
	printf("%d th image: %s\n",i,image_indice[i]);

	readInData(string(image_indice[i]));
	printf("dimension is: %d %d %d\n",dims[0],dims[1],dims[2]);
	if(dims[0]*dims[1]*dims[2]>100000)  //trained on 65 volumes  
	{ free3DMat<float>(img); free3DMat<uchar>(gt); printf("image is too large, skip\n"); continue;}

	int box[6]={0,0,0,dims[0]-1,dims[1]-1,dims[2]-1};
	bmin=min3D<float>(img,box); 
	bmax=max3D<float>(img,box);
	//int binNum=(int)(bmax-bmin+1);
	int binNum=70;

	clock_t t1,t2;
	t1=clock();
	segmentation(img,gt,dims,box_fg,spacing,lambda_xy,lambda_z,binNum);
	t2=clock();
    float diff ((float)t2-(float)t1);
	printf("running time is: %f sec\n\n",diff/CLOCKS_PER_SEC);
   
	// free allocated memory
	free3DMat<float>(img);
	free3DMat<uchar>(gt);
    
   }

   /*
   FILE* fidFM=fopen(outFM.c_str(),"wt");
   FILE* fidRE=fopen(outRE.c_str(),"wt");
   FILE* fidPE=fopen(outPE.c_str(),"wt");
   for(int j=0;j<88;j++){
	   for(int k=0;k<20;k++){
		   fprintf(fidFM,"%f ",fmeasure[j][k]);
		   fprintf(fidRE,"%f ",recall[j][k]);
		   fprintf(fidPE,"%f ",precision[j][k]);
	   }
	   fprintf(fidFM,"\n");
	   fprintf(fidRE,"\n");
	   fprintf(fidPE,"\n");
   }
   fclose(fidFM);
   fclose(fidRE);
   fclose(fidPE);
   */

   //float avgRecall=std::accumulate(recall.begin(),recall.end(),0.0)/recall.size();
   //float avgPrecision=std::accumulate(precision.begin(),precision.end(),0.0)/precision.size();
   //float avgFmeasure=std::accumulate(fmeasure.begin(),fmeasure.end(),0.0)/fmeasure.size();
   //printf("average recall: %f, precision: %f, fmeasure: %f\n",avgRecall, avgPrecision, avgFmeasure);

	return 1;
}

