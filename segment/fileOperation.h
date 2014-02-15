#ifndef _FILEOPERATION_H_
#define _FILEOPERATION_H_


#include"segment_3d.hpp"
#include "Matrix3D.hpp"
#include<stdio.h>
#include<vector>
#include<string>
using namespace std;


extern float ***img;
extern uchar ***gt;
extern int box_fg[6];
extern int tightBox[6];
extern int dims[3];
extern float spacing[3];
extern char image_indice[88][30];

void readInData(string filename)
{
	// input path that remain the same for different methods and parameters
    string imgPath="D:/workspace/project/crop3D/cropBinary/VOLUME/";
    string gtPath="D:/workspace/project/crop3D/cropBinary/GT/";
    string infoPath="D:/workspace/project/crop3D/cropBinary/INFO/";
	//get full name for input files
	string imgName,gtName,infoName;
	imgName=imgPath+filename+".txt";
	gtName=gtPath+filename+".txt";
	infoName=infoPath+filename+".txt";


	//read foreground box
	FILE *fid=fopen(infoName.c_str(),"rt");
	if(fid==NULL) { printf("failed to read info\n"); exit(0);}
	fscanf(fid,"%d %d %d\n",&dims[0],&dims[1],&dims[2]);
	fscanf(fid,"%f %f %f\n",&spacing[0],&spacing[1],&spacing[2]);
	fscanf(fid,"%d %d %d %d %d %d",&box_fg[0],&box_fg[1],&box_fg[2],&box_fg[3],&box_fg[4],&box_fg[5]);
	fscanf(fid,"%d %d %d %d %d %d",&tightBox[0],&tightBox[1],&tightBox[2],&tightBox[3],&tightBox[4],&tightBox[5]);
	fclose(fid);
	for(int ib=0;ib<6;ib++) { box_fg[ib]=box_fg[ib]-1; tightBox[ib]=tightBox[ib]-1;}
	img=read3DMat<float>(imgName.c_str(),dims);
	gt=read3DMat<uchar>(gtName.c_str(),dims);
}

void writeData(string path,string filename,int dims[3],int segBox[6],uchar ***labeling,float ***dataterm, float ***smoothterm, vector<vector<int> > &seeds2D, bool uint8)
{
	//write result: segmentation, smoothterm, dataterm
	string smoothName, dataName,segName,seedsName;
	smoothName=path+"smoothterm/"+filename+".txt";
	dataName=path+"dataterm/"+filename+".txt";
	segName=path+"segmentation/"+filename+".txt";
	seedsName=path+"seeds/"+filename+".txt";

	if(smoothterm){
		if(uint8) {
			normalize3D<float>(smoothterm,segBox,0,255);
			uchar*** smoothUint8=create3DMat<uchar>(dims,0);
			set3DMat<uchar,float>(smoothUint8,smoothterm,segBox);
			write3DMat<uchar>(smoothName.c_str(),smoothUint8,dims);
			free3DMat<uchar>(smoothUint8);
		}
		else write3DMat<float>(smoothName.c_str(),smoothterm,dims);
	}
	if(dataterm){
		if(uint8) {
			normalize3D<float>(dataterm,segBox,0,255);
			uchar*** dataUint8=create3DMat<uchar>(dims,0);
			set3DMat<uchar,float>(dataUint8,dataterm,segBox);
			write3DMat<uchar>(dataName.c_str(),dataUint8,dims);
			free3DMat<uchar>(dataUint8);
		}
		else write3DMat<float>(dataName.c_str(),dataterm,dims);
	}
	if(labeling)
	write3DMat<uchar>(segName.c_str(),labeling,dims);

	FILE *fid=fopen(seedsName.c_str(),"wt");
	fprintf(fid,"%d 3 0\n",seeds2D.size());
	for(int ei=0;ei<seeds2D.size();ei++)
	{
       fprintf(fid,"%d %d\n",seeds2D[ei][0],seeds2D[ei][1]);
	}
	fclose(fid);
}

void readImgIndice(char *file)
{
	FILE *fid=fopen(file,"rt");
	for(int i=0;i<88;i++)
		 fscanf(fid,"%s\n",image_indice[i]);
	fclose(fid);

}

#endif