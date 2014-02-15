
#ifndef __MAT_H__
#define __MAT_H__

/*------------------------------------------------------------------------------------------*/
template<typename T>
class Mat2D{
public:
	Mat2D(int rows,int cols,T value=0);

public:
	T data[200][200];
	int dims[2];
};

/*------------------------------------------------------------------------------------------*/
template<typename T>
class Mat3D{
public:
	Mat3D(int rows,int cols,int slices,T value=0);
	void normalize(T min, T max);
public:
	T data[200][200][100];
    int dims[3];
};

#endif