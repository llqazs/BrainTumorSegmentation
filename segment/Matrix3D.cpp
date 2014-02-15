#include"Matrix3D.h"
#include<stdlib.h>
#include<stdio.h>
template<typename T>
void writeMatrix3D(char *file,T ***matrix,int dims[],char *format)
{
	if(matrix==NULL) { printf("data to be write dose not exist\n"); exit(0);}
	FILE* out=fopen(file,"w");
	if(out==NULL) { printf("%s does not exist\n"); exit(0); }
	fprintf(out,"%d %d %d\n",dims[0],dims[1],dims[2]);
	for(int r=0;r<dims[0];r++)
	{
		for(int s=0;s<dims[2];s++)
			for(int c=0;c<dims[1];c++)
			{
				fprintf(out,format,matrix[r][c][s]);
			}

	    fprintf(out,"\n");
	}
	fclose(out);
}

template <typename T>
T ***CreateMatrix3D(int dims[],T value)
{
   int rows=dims[0]; int cols=dims[1]; int slices=dims[2];
   T ***matrix=NULL;
   matrix=(T ***)malloc(sizeof(T**)*rows);
   for(int r=0;r<rows;r++)
   {
	   matrix[r]=(T **)malloc(sizeof(T*)*cols);
	   for(int c=0;c<cols;c++)
	   {
		   matrix[r][c]=(T *)malloc(sizeof(T)*slices);
		   for(int s=0;s<slices;s++)
		   {
			   matrix[r][c][s]=value;
		   }
	   }
   }

   return matrix;

}


template <typename T>
void FreeMatrix3D(T ***matrix, int dims[])
{
	if(matrix!=NULL)
	{
    int rows=dims[0]; int cols=dims[1]; int slices=dims[2];
	if(matrix!=NULL)
	{
		for(int r=0;r<rows;r++)
		{
			for(int c=0;c<cols;c++)
			{
				free(matrix[r][c]);
			}
			free(matrix[r]);
		}
		free(matrix);
	}
	}
}

void normalize3D(float ***matrix,int size[],float vmin, float vmax)
{
	if(matrix==NULL) {
		printf("empty matrix!\n");
		exit(0);
	}

	int rows=size[0]; int cols=size[1]; int slices=size[2];
	float min=1; float max=0;
	for(int r=0;r<rows;r++)
		for(int c=0;c<cols;c++)
			for(int s=0;s<slices;s++)
			{
				if(matrix[r][c][s]>max) max=matrix[r][c][s];
				if(matrix[r][c][s]<min) min=matrix[r][c][s];
			}
    float k=(vmax-vmin)/(max-min);
	for(int r=0;r<rows;r++)
		for(int c=0;c<cols;c++)
			for(int s=0;s<slices;s++)
			{
                matrix[r][c][s]=k*(matrix[r][c][s]-min)+vmin;
			}

}