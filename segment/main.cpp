#include <stdlib.h>
#include <stdio.h>

template<typename T>
T ***matrix(int dims[],T value)
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

int main()
{
	int dims[3]={50,20,50};
	float ***img=matrix<float>(dims,3.0);
	for(int r=0;r<dims[0];r++)
		for(int c=0;c<dims[1];c++)
			for(int s=0;s<dims[2];s++)
			{
				// if(img[r][c][s]!=3.0)
					printf("%f ",img[r][c][s]);
			}


	return 0;


}