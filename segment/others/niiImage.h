#ifndef __NIIIMAGE_H__
#define __NIIIMAGE_H__

/*------------------------------------------------------------------------------------------*/
// class for nii image, inherited from class Mat3D
template<typename T>
class niiImage: public Mat3D{
public:
	int spacing[3];
};

#endif