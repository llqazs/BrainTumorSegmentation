#include"segment_3d.hpp"
#include "Matrix3D.hpp"
#include<stdio.h>
#include<vector>
using namespace std;


extern float ***img;
extern bool ***gt;
extern int box_fg[6];
extern float spacing[3];
extern char image_indice[88][30];

void getFileFullName(char *filefullname,char *path, char *filename, char *extension)
{
	strcpy(filefullname,path);
	strcat(filefullname,filename);
	strcat(filefullname,extension);
}

void readIMG(char *file,char *format,int *dims)
{
	FILE *fid=fopen(file,"rt");
	if(fid==NULL) { printf("file %s does not exist\n",file); exit(0); }
	int rows,cols,slices;
	fscanf(fid,"%d %d %d\n",&rows,&cols,&slices);
	dims[0]=rows; dims[1]=cols; dims[2]=slices; // ????
	img=CreateMatrix3D<float>(dims,0);
	float v;
	for(int r=0;r<rows;r++)
	{
		for(int s=0;s<slices;s++)
		{
			for(int c=0;c<cols;c++)
			{
			//	printf("%d %d %d\n",r,c,s);
				fscanf(fid,format,&v);
				img[r][c][s]=v;
			}
		}
	}
	fclose(fid);
}

void readGT(char *file,char *format,int *dims)
{
	FILE *fid=fopen(file,"rt");
	if(fid==NULL) { printf("file %s does not exist\n",file); exit(0); }
	int rows,cols,slices;
	fscanf(fid,"%d %d %d\n",&rows,&cols,&slices);
	dims[0]=rows; dims[1]=cols; dims[2]=slices;
	gt=CreateMatrix3D<bool>(dims,0);
	int v;
	for(int r=0;r<rows;r++)
	{
		for(int s=0;s<slices;s++)
		{
			for(int c=0;c<cols;c++)
			{
			//	printf("%d %d %d\n",r,c,s);
				fscanf(fid,format,&v);
				gt[r][c][s]=(bool)v;
			}
		}
	}
	fclose(fid);
}

void readInData(char *filename)
{
	// input path that remain the same for different methods and parameters
	char img_path[300]="D:/workspace/project/crop3D/recrop/TOMOB2INT_CROP/";
	char gt_path[300]="D:/workspace/project/crop3D/recrop/SEG_CROP/";
	char fgbx_path[300]="D:/workspace/project/crop3D/recrop/Foreground_box/";
	char space_path[300]="D:/workspace/project/data/SPACING/";
	//get full name for input files
	char img_name[300], gt_name[300], fgbx_name[300],space_name[300];
	getFileFullName(img_name,img_path,filename,".txt");
	getFileFullName(gt_name,gt_path,filename,".txt");
	getFileFullName(fgbx_name,fgbx_path,filename,".txt");
	getFileFullName(space_name,space_path,filename,".txt");

	int dims[3];
	readIMG(img_name,"%f ",dims);  //read image
	readGT(gt_name,"%d ",dims);   //read ground truth
	//read foreground box
	FILE *fid=fopen(fgbx_name,"rt");
	if(fid==NULL) { printf("failed to read %s\n",fgbx_name); exit(0);}
	fscanf(fid,"%d %d %d %d %d %d",&box_fg[0],&box_fg[1],&box_fg[2],&box_fg[3],&box_fg[4],&box_fg[5]);
	fclose(fid);
	for(int ib=0;ib<6;ib++) box_fg[ib]=box_fg[ib]-1;

	// read spacing
	fid=fopen(space_name,"rt");
	fscanf(fid,"%f %f %f",&spacing[0],&spacing[1],&spacing[2]);
    fclose(fid);
}

void writeData(char *path,char *filename,int dims[],bool ***labeling,float ***dataterm, float ***smoothterm, vector<vector<int> > &seeds2D)
{
	//write result: segmentation, smoothterm, dataterm
	char smooth_name[300], data_name[300], seg_name[300], seeds_name[300];
	char smooth_path[300], data_path[300], seg_path[300], seeds_path[300];
	strcpy(smooth_path,path);
	strcat(smooth_path,"smoothterm/");
	getFileFullName(smooth_name,smooth_path,filename,".txt");

	strcpy(data_path,path);
	strcat(data_path,"dataterm/");
	getFileFullName(data_name,data_path,filename,".txt");

	strcpy(seg_path,path);
	strcat(seg_path,"segmentation/");
	getFileFullName(seg_name,seg_path,filename,".txt");

	strcpy(seeds_path,path);
	strcat(seeds_path,"seeds/");
	getFileFullName(seeds_name,seeds_path,filename,".txt");

    writeMatrix3D<float>(smooth_name,smoothterm,dims,"%.5f ");
	writeMatrix3D<float>(data_name,dataterm,dims,"%.5f ");
	writeMatrix3D<bool>(seg_name,labeling,dims,"%d ");

	FILE *fid=fopen(seeds_name,"wt");
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

//void readFilePath()
//{
//   char img_path[300]="D:/workspace/project/crop3D/recrop/TOMOB2INT_CROP/";
//   char gt_path[300]="D:/workspace/project/crop3D/recrop/SEG_CROP/";
//   char fgbx_path[300]="D:/workspace/project/crop3D/recrop/Foreground_box/";
//   char space_path[300]="D:/workspace/project/data/SPACING/";
//   //parzenwin  no_small_const encourage_edge
//   char smooth_path[300]="D:/workspace/project/result/result_0804/fixed_bin/smoothterm/";
//   char data_path[300]="D:/workspace/project/result/result_0804/fixed_bin/dataterm/";
//   char seg_path[300]="D:/workspace/project/result/result_0804/fixed_bin/segmentation/";
//   char seeds_path[300]="D:/workspace/project/result/result_0804/fixed_bin/seeds/";
//
//   FILE *fid_indice=fopen("D:/workspace/project/crop3D/recrop/TOMOB2INT_CROP/image_name.txt","rt");
//   for(int i=0;i<88;i++)
//   {
//	   fscanf(fid_indice,"%s\n",image_indice[i]);
//	   strcpy(img_name[i],img_path);
//	   strcat(img_name[i],image_indice[i]);
//	   strcat(img_name[i],".txt");
//
//	   strcpy(gt_name[i],gt_path);
//	   strcat(gt_name[i],image_indice[i]);
//	   strcat(gt_name[i],".txt");
//
//	   strcpy(fgbx_name[i],fgbx_path);
//	   strcat(fgbx_name[i],image_indice[i]);
//	   strcat(fgbx_name[i],".txt");
//
//	   strcpy(space_name[i],space_path);
//	   strcat(space_name[i],image_indice[i]);
//	   strcat(space_name[i],".txt");
//
//	   strcpy(smooth_name[i],smooth_path);
//	   strcat(smooth_name[i],image_indice[i]);
//	   strcat(smooth_name[i],".txt");
//
//	   strcpy(data_name[i],data_path);
//	   strcat(data_name[i],image_indice[i]);
//	   strcat(data_name[i],".txt");
//
//	   strcpy(seg_name[i],seg_path);
//	   strcat(seg_name[i],image_indice[i]);
//	   strcat(seg_name[i],".txt");
//
//	   strcpy(seeds_name[i],seeds_path);
//	   strcat(seeds_name[i],image_indice[i]);
//	   strcat(seeds_name[i],".txt");
//   }
//   fclose(fid_indice);
//
//}