
#include "Mat.h"

template<typename T>
void setValueInsideCube(Mat3D *matrix,int box[],T value)
{
	int *dims=matrix->dims; //?????
	if(matrix==NULL) { printf("error in function setVauleInsideCube, matrix is empty\n"); exit(0); }
	if(box[0]<0 || box[1]<0 || box[2]<0 ||box[3]>dims[0]-1 ||box[4]>dims[1]-1 || box[5]>dims[2]-1)
	{ printf("error in function segVauleInsideCube, cube exceeds the range of data\n"); exit(0); }
	for(int r=box[0];r<=box[3];r++)
		for(int c=box[1];c<=box[4];c++)
			for(int s=box[2];s<=box[5];s++)
				matrix->data[r][c][s]=value;

}

inline float funcExp(float diff, float variance, float labmda)
{
	return lambda*exp(-diff*diff/(variance*variance));
}