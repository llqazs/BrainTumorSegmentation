// file that handles 3D matrix creation, reading and writing
/*mat[s][r][c], s: index of slice
                r: index of row
				c: index of column
				*/


#ifndef _MATRIX3D_HPP_
#define _MATRIX3D_HPP_

#include<stdlib.h>
#include<stdio.h>


template<typename T>
T ***build3DMat(T * arr, int dims[3])
{
	int nRows=dims[0]*dims[1];
	T **rows=(T **) malloc(sizeof(T*)*nRows);
	T *pRow=arr;
	for(int i=0;i<nRows;i++)
	{  
		rows[i]=pRow;
		pRow=pRow+dims[2]; 
	}

	T ***slices=(T***) malloc(sizeof(T**)*dims[0]);
	T **pSli=rows;
	for(int i=0;i<dims[0];i++)
	{
		slices[i]=pSli;
		pSli=pSli+dims[1];
	}

	return slices;

}

template<typename T>
void free3DMat(T ***mat)
{
	if(mat[0][0][0]) free(mat[0][0]);
	if(mat[0][0]) free(mat[0]);
	if(mat) free(mat);
}

template<typename T>
T*** create3DMat(int dims[3],T val=0)
{
	int count=dims[0]*dims[1]*dims[2];
	T* arr=(T*) malloc(sizeof(T)*count);
	for(int i=0;i<count;i++) arr[i]=val;
	T*** mat=build3DMat(arr,dims);
	return mat;

}

template<typename T>
T ***read3DMat(const char *file,int dims[3])
{

 FILE * pFile;
  long lSize;
  char * buffer;
  size_t result;

  pFile = fopen ( file , "rb" );
  if (pFile==NULL) {fputs ("File error",stderr); exit (0);}

  // obtain file size:
  fseek (pFile , 0 , SEEK_END);
  lSize = ftell (pFile);
  rewind (pFile);

  // allocate memory to contain the whole file:
  buffer = (char*) malloc (sizeof(char)*lSize);
  if (buffer == NULL) {fputs ("Memory error",stderr); exit (1);}

  // copy the file into the buffer:
  result = fread (buffer,1,lSize,pFile);
  if (result != lSize) {fputs ("Reading error",stderr); exit (2);}

  /* the whole file is now loaded in the memory buffer. */

  // terminate
  fclose (pFile);
  T *arr;

  arr=(T*)buffer;

  T ***mat=build3DMat(arr,dims);

	return mat;
}


template<typename T>
void write3DMat(const char *file,T* mat, int dims[3])
{
	if(file==NULL) {printf("No file for written\n"); exit(0);}
	if(mat==NULL) {printf("empty matrix\n"); exit(0);}

	int size=dims[0]*dims[1]*dims[2];

	FILE* pFile=fopen ( file , "wb" );
    if (pFile==NULL) {printf("File error\n"); exit (0);}

	fwrite(mat,sizeof(T),(unsigned int)(dims[0]*dims[1]*dims[2]),pFile);
	fclose(pFile);

}


template<typename T>
T min3D(T ***matrix,int box[6])
{
	T min=matrix[box[0]][box[1]][box[2]];
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				if(matrix[s][r][c]<min) min=matrix[s][r][c];
			}

	return min;

}


template<typename T>
T max3D(T ***matrix,int box[6])
{
	T max=matrix[box[0]][box[1]][box[2]];
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<=box[5];c++)
			{
				if(matrix[s][r][c]>max) max=matrix[s][r][c];
			}

	return max;

}


template<typename T>
void normalize3D(T ***matrix, int box[6], T outMin, T outMax)
{
	T min=min3D(matrix,box);
	T max=max3D(matrix,box);
	T range=max-min;
	if(range==0) { printf("max==min in normalize3D, nothing is done\n"); return; }
	T outRange=outMax-outMin;
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<box[5];c++)
			{
				matrix[s][r][c]=outRange*(matrix[s][r][c]-min)/range+outMin;
			}

}

template<typename T>
void set3DMat(T ***matrix,int box[6], T val)
{
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<box[5];c++)
			{
				matrix[s][r][c]=val;
			}

}

template<typename Tdst, typename Tsrc>
void set3DMat(Tdst ***dst,Tsrc ***src,int box[6])
{
	for(int s=box[0];s<=box[3];s++)
		for(int r=box[1];r<=box[4];r++)
			for(int c=box[2];c<box[5];c++)
			{
				dst[s][r][c]=(Tdst)src[s][r][c];
			}

}

//set a slice of the mat dst from a slice of src
template<typename T>
void set3DMat(T ***dst, T ***src, int dst_z, int src_z, int dims[3])
{
	for(int r=0;r<dims[1];r++)
		for(int c=0;c<dims[2];c++)
		{
			dst[dst_z][r][c]=src[src_z][r][c];
		}

}


#endif