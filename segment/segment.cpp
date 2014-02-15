
#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.VC80.CRT' version='8.0.50727.762' processorArchitecture='X86' publicKeyToken='1fc8b3b9a1e18e3b' language='*'\"")
#include <stdio.h>
#include <vector>
#include <math.h>
#include <queue> 
#include "NTadt.h"
#include "graph.h"

using namespace std;

typedef Graph<float, float, float> GraphType;
typedef int node_id;
typedef int bin_id;
#define  HUGE         10000.0
#define SMALL_CONST 0.0001
float STAR_HUGE=1000.0;
// the structure to store pixel after shuffling

int bin_num; 
float bmin, bmax;

//----------------------------------------function declaration--------------------------------------------------------------
void check_error( int boolE, const char *error_message);
inline float computeNLinkCost(float diff,float variance, float lambda);
void gray2float(const GrayImage &grayImage, FloatImage &floatImage);
float  computeImageVariance(FloatImage image);
inline node_id getPixelID(int x, int y, int width);
void normalizeImage(FloatImage &image, float vmin,float vmax);
int getBinID(float x);
void float2gray(const FloatImage &fimage, GrayImage &gimage);
void histogram(vector<float> &imArray, vector<int> &hist, vector<int> &idx);
void computeTWeight(GraphType *g, GrayImage &dataterm, const vector<int> &pixelBinID, const vector<float> &fgProb, const vector<float> &bkProb);
void computeAppearanceModel(const FloatImage &image, const BinaryImage &mask, vector<float> &pixelArray, vector<int> &pixelBinID, vector<float> &fgProb, vector<float> &bkProb);
void computeNWeight(GraphType *g,const FloatImage &image,GrayImage &smooth, float variance,float lambda);
void segmentationToImage(GraphType *g, GrayImage &segmentation);
void segment(const FloatImage &image, GrayImage &segmentation, GrayImage &dataterm, GrayImage &smooth, float lambda, int binNum);
void postProcessing(GrayImage &segmentation);
float computeDiff(const GrayImage &segmentation, const GrayImage &gt);

//--------------------------------------function definition-------------------------------------------------------------------
void check_error( int boolE, const char *error_message)
{ 
   if  (boolE) 
   {
      printf("\n %s \n", error_message);
      exit(1);
    }
}

inline float computeNLinkCost(float diff,float variance, float lambda)
{
	return lambda*exp(-diff*diff/(variance*variance));

}

void gray2float(const GrayImage &grayImage, FloatImage &floatImage)
{
	int width=imGetWidth(grayImage); 
	int height=imGetHeight(grayImage);
	for(int x=0;x<width; x++)
		for(int y=0; y<height; y++) 
			imRef(floatImage,x,y)=(float)imRef(grayImage,x,y);

}

//--------------------------------------------------------------------------------
float  computeImageVariance(FloatImage image)
{
	int width=imGetWidth(image); int height=imGetHeight(image);
	float DEFAULT_VARIANCE = (float) 0.0001;
	float v = (float) 0.0;
	int total = 0;

	for ( int y = 1; y < height; y++ )
		for ( int x = 1; x < width ; x++ ){
			v = v + abs(imRef(image,x,y)-imRef(image,x-1,y))+abs(imRef(image,x,y)-imRef(image,x,y-1));
			total = total + 2;
		}
  //  printf("\nImage variance is %f",v/total); getchar();
	if (v/total < DEFAULT_VARIANCE ) 
		return(DEFAULT_VARIANCE);
	else return( v/total);
}

//------------------------------------------------------------------------------
inline node_id getPixelID(int x, int y, int width) 
{
	return(y*width+x);  // x is col and y is row

}

void normalizeImage(FloatImage &image, float vmin,float vmax)  // vmin and vmax is the min and max value after normalization
{
  float min,max;
  min=max=imRef(image,0,0);
  int width,height;
  width=imGetWidth(image);
  height=imGetHeight(image);
  for(int x=0;x<width;x++)
	  for(int y=0;y<height;y++) {
		  if(imRef(image,x,y)<min) min=imRef(image,x,y);
		  if(imRef(image,x,y)>max) max=imRef(image,x,y);
	  }
  float k=(vmax-vmin)/(max-min);
  for(int x=0;x<width;x++)
	  for(int y=0;y<height;y++) {
		  imRef(image,x,y)=k*(imRef(image,x,y)-min)+vmin;
	  }
}


int getBinID(float x)
{
	bin_id id=(int)(bin_num*(x-bmin)/(bmax-bmin));
	if(id==bin_num) id=bin_num-1;
	return id;
}


void histogram(vector<float> &imArray, vector<int> &hist, vector<int> &idx)
{
	int length=imArray.size();
    int bin_num=hist.size();
	for(int i=0;i<length;i++) {
	bin_id id=getBinID(imArray[i]);
	hist[id]+=1;
	idx[i]=id;
	}

}

//-----------------------------------------------------------------------------
void computeTWeight(GraphType *g, GrayImage &dataterm, const vector<int> &pixelBinID, const vector<float> &fgProb, const vector<float> &bkProb)
{
    int num_pixels=pixelBinID.size();
	int width=imGetWidth(dataterm);
	int height=imGetHeight(dataterm);
	float r; 
	FloatImage likelihood=(FloatImage) imNew(IMAGE_FLOAT,width,height);
	for(int i=0;i<num_pixels;i++) {
		float fgCost=-log(fgProb[pixelBinID[i]]);
		float bkCost=-log(bkProb[pixelBinID[i]]);
		//if(fgCost>10.0) { fgCost=10.0; bkCost=5.0; }
		//if(bkCost>10.0) { bkCost=10.0; fgCost=5.0; }
		g->add_tweights(i,bkCost,fgCost);
	//	r=g->get_trcap(i);
		int y=i/width;
		int x=i-y*width;
		r=bkCost/fgCost;
		imRef(likelihood,x,y)=bkCost-fgCost; // bkCost/fgCost;
	}
	normalizeImage(likelihood,0,255);
	float2gray(likelihood,dataterm);
	imFree(likelihood);
}


void computeAppearanceModel(const FloatImage &image, const BinaryImage &mask, vector<float> &pixelArray, vector<int> &pixelBinID, vector<float> &fgProb, vector<float> &bkProb)
{
	int width=imGetWidth(image);
	int height=imGetHeight(image);
	int num_pixels=width*height;
	vector<bool> maskArray(num_pixels);
	// get bin id for every pixel in the image
    for(int x=0; x<width; x++)
		for(int y=0; y<height; y++){
			node_id pixelID=getPixelID(x,y,width);
			pixelArray[pixelID]=imRef(image,x,y);
			maskArray[pixelID]=imRef(mask,x,y); // set flag for maskArray
			bin_id binID=getBinID(imRef(image,x,y));
			pixelBinID[pixelID]=binID;
		}

	// reset fg and bk histogram to zero
    for(int i=0; i<bin_num; i++) {
		fgProb[i]=0;
		bkProb[i]=0;
	}

	int num_fgpixels=0; int num_bkpixels=0;
	for(int i=0; i<num_pixels; i++) {
		bin_id binID=getBinID(pixelArray[i]);
		if(maskArray[i]==1) {
			fgProb[binID]+=1;
			num_fgpixels+=1;
		}
		else {
			bkProb[binID]+=1;
			num_bkpixels+=1;
		}
	}
	
	for(int i=0;i<bin_num; i++) {
		fgProb[i]=fgProb[i]/num_fgpixels+SMALL_CONST;
		bkProb[i]=bkProb[i]/num_bkpixels+SMALL_CONST;
	}

}




void computeNWeight(GraphType *g,const FloatImage &image,GrayImage &smooth, float variance,float lambda)
{
	
	int width=imGetWidth(image);
	int height=imGetHeight(image);
    FloatImage smoothImage=(FloatImage) imNew(IMAGE_FLOAT,width,height);  
   
	//reset smoothImage
    for(int x=0; x<width; x++)
	   for(int y=0; y<height; y++) imRef(smoothImage,x,y)=0;

	// n weight for direction -
	for(int y=0;y<height;y++) {
		for(int x=0;x<width-1;x++) {
			node_id left=getPixelID(x,y,width);
			node_id right=left+1;
			float diff=abs(imRef(image,x,y)-imRef(image,x+1,y));
			float weight=computeNLinkCost(diff,variance,lambda);
			g->add_edge(left,right,weight,weight);
			imRef(smoothImage,x,y)+=weight;
			imRef(smoothImage,x+1,y)+=weight;
		}
	}
	// n weight for direction |
	for(int x=0;x<width;x++) {
		for(int y=0;y<height-1;y++) {
			node_id up=getPixelID(x,y,width);
			node_id down=up+width;
			float diff=abs(imRef(image,x,y)-imRef(image,x,y+1));
			float weight=computeNLinkCost(diff,variance,lambda);
			g->add_edge(up,down,weight,weight);
			imRef(smoothImage,x,y)+=weight;
			imRef(smoothImage,x,y+1)+=weight;
		}
	}
	/* n weight for direction \ */
	for(int y=0;y<height-1;y++) {
		for(int x=0;x<width-1;x++) {
			node_id left=getPixelID(x,y,width);
			node_id right=left+1+width;
			float diff=abs(imRef(image,x,y)-imRef(image,x+1,y+1));
			float weight=0.707*computeNLinkCost(diff,variance,lambda);
			g->add_edge(left,right,weight,weight);
			imRef(smoothImage,x,y)+=weight;
			imRef(smoothImage,x+1,y+1)+=weight;
		}
	}
	// n weight for direction /
	for(int y=1;y<height;y++) {
		for(int x=0;x<width-1;x++) {
			node_id left=getPixelID(x,y,width);
			node_id right=left+1-width;
			float diff=abs(imRef(image,x,y)-imRef(image,x+1,y-1));
			float weight=0.707*computeNLinkCost(diff,variance,lambda);
			g->add_edge(left,right,weight,weight);
			imRef(smoothImage,x,y)+=weight;
			imRef(smoothImage,x+1,y-1)+=weight;
		}
	}

	normalizeImage(smoothImage,0,255);
	float2gray(smoothImage,smooth);
	imFree(smoothImage);
}

void setMask(BinaryImage &mask, int size_fg)
{
	int width=imGetWidth(mask);
	int height=imGetHeight(mask);
	int cx=width/2; int cy=height/2;
	int left=cx-size_fg/2; int right=cx+size_fg/2;
	int top=cy-size_fg/2; int bottom=cy+size_fg/2;
	for(int x=0; x<width; x++)
		for(int y=0; y<height; y++) {
			imRef(mask,x,y)=0;
			if(x<=right && x>= left && y<=bottom && y>=top) imRef(mask,x,y)=1;
		}

}

void hardConstraint(GraphType *g, int width, int height)
{
	int size_center=(width+height)/8; int size_border=size_center/3;
	//add hard constraints for pixels at border
	for(int y=0; y<size_border; y++) 
		for(int x=0; x<width; x++) {
			node_id pixelID=getPixelID(x,y,width);
			g->add_tweights(pixelID,0,HUGE);  // add hard constraints for top rows
			pixelID=getPixelID(x,height-1-y,width);
			g->add_tweights(pixelID,0,HUGE); // add hard constraints for bottom rows
		}
	for(int x=0; x<size_border; x++)
		for(int y=size_border; y<height-size_border; y++) {
			node_id pixelID=getPixelID(x,y,width);
			g->add_tweights(pixelID,0,HUGE);  // add hard constraints for left colunms
			pixelID=getPixelID(width-1-x,y,width);
			g->add_tweights(pixelID,0,HUGE);  // add hard constraints for right colunms
		}
    // add hard constraints for pixels in the center
	int left=width/2-size_center/2; int top=height/2-size_center/2;
	for(int x=left; x<left+size_center; x++)
		for(int y=top; y<top+size_center; y++) {
			node_id pixelID=getPixelID(x,y,width);
			g->add_tweights(pixelID,HUGE,0);
		}

}

void float2gray(const FloatImage &fimage, GrayImage &gimage)
{
	int w=imGetWidth(fimage);
	int h=imGetHeight(fimage);
	//if(gimage==(GrayImage) NULL) gimage=(GrayImage) imNew(IMAGE_GRAY,w,h);
	for(int x=0;x<w;x++)
		 for(int y=0;y<h;y++)
			 imRef(gimage,x,y)=(unsigned char)imRef(fimage,x,y);

}

void addNextConstraintPixels(GraphType *g, int x,int y,int width,int SeedX,int SeedY,
					         vector<int> &pixels,vector<unsigned char> &markings, 
							 RGB color, RGBImage &im,int height)
							 
{
	float a,b;
	int nextX;
	int numPix = 0;
	bool reverseOrder = false;
	
	
	// solve for line coefficients a and b first
	if ( SeedX != x ){
		a = ((float) (SeedY-y))/(SeedX-x);
		b = y-a*x;
	}

	// next get all pixels on the discretized line. 4 cases
	if ( abs(y-SeedY) > abs(SeedX-x) ){
		int startY,endY;
		if ( y > SeedY) {
			startY = SeedY+1;
			endY = y;
		}
		else  {
			startY = y;
			endY = SeedY-1;
			reverseOrder = true;
		}
			
		for ( int nextY = startY; nextY <= endY; nextY++ ){
				if ( SeedX != x) 
					nextX = (int) ((0.5+nextY-b)/a);
				else nextX = SeedX; 

				pixels[numPix] = nextX+nextY*width;
				numPix++;
			}
	}
	else{
		int startX,endX;
		if ( x > SeedX){
			startX = SeedX+1;
			endX   = x;
		}
		else{
			startX = x;
			endX = SeedX-1;
			reverseOrder = true;
		}
		for ( int nextX = startX; nextX <= endX; nextX++ ){
			int nextY = (int) (a*nextX+b+0.5);
			pixels[numPix] = nextX+nextY*width;
	
			numPix++;
		}
	}
	
	int DISPLAY = 1;
	
	if ( reverseOrder ){
		int i;
		for (  i = 1; i < numPix; i++ ){
			if ( markings[pixels[i-1]] == 0 ){
				markings[pixels[i-1]] = 1;
				g->add_edge(pixels[i-1],pixels[i],STAR_HUGE,0);
				if ( (y+x) == ((y+x)/DISPLAY)*DISPLAY ){
					imRef(im,pixels[i-1]-width*(pixels[i-1]/width),pixels[i-1]/width) = color;
					//imSave(im,"connections.ppm");
				}
			}
			else break;
		}
		
	}
	else {
		int i;
		for (  i  = numPix-1; i > 0; i-- ){
			if ( markings[pixels[i]] == 0 ){
				markings[pixels[i]] = 1;
				g->add_edge(pixels[i],pixels[i-1],STAR_HUGE,0);

				if ( (y+x) == ((y+x)/DISPLAY)*DISPLAY ){
					imRef(im,pixels[i]-width*(pixels[i]/width),pixels[i]/width) = color;
					//imSave(im,"connections.ppm");
				}
			}
			else break;
		}
		
	}
	
	//saveArray(markings,width,height);
	//imSave(im,"connections.ppm");		
	
}

//---------------------------------------------------------
// chooses all image pixels (not just the border) in random order
void setStarTlinksAllPixels(GraphType *g,int width,
				   int height, int seedX, int seedY,int num_pixels){

	vector<unsigned char> markings(num_pixels,0);
	vector<int> pixelOrder(num_pixels);
	vector<int> pixels(num_pixels);
	int numPixels = 0;

	markings[seedX+seedY*width] = 1;

	RGBImage connections = (RGBImage) imNew(IMAGE_RGB,width,height);

	for ( int i = 0; i < num_pixels; i++ )
		pixelOrder[i] = i;

	for ( int i = 0; i < 3*num_pixels; i++ ){
		int first  = rand()*rand()%num_pixels;
		int second = rand()*rand()%num_pixels;
		int temp = pixelOrder[first];
		pixelOrder[first]  = pixelOrder[second];
		pixelOrder[second] = temp;
	}

	for ( int i = 0; i < num_pixels; i++ )
		if ( markings[pixelOrder[i]] == 0 )
		{

			int nextY = pixelOrder[i]/width;
			int nextX = pixelOrder[i]-nextY*width;

			RGB color = {rand()%256,rand()%256,rand()%256};
			addNextConstraintPixels(g,nextX,nextY,width,seedX,seedY,
			                    pixels,markings,color,connections,height);
		}
	int unmarked=0;
	for(int i=0; i<num_pixels; i++) {
		if(markings[i]==0){ unmarked++; }
	}
	imFree(connections);
//	printf("unmarked pixels: %d\n",unmarked);

//	imSave(connections,"connections.pgm");
}

void segmentationToImage(GraphType *g, GrayImage &segmentation)
{
   int width      = imGetWidth(segmentation);
   int height     = imGetHeight(segmentation);
	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++) {
			node_id id=getPixelID(x,y,width);
			if(g->what_segment(id)== GraphType::SOURCE)
			imRef(segmentation,x,y)=255;
			else imRef(segmentation,x,y)=0;
		}

}

float computeDiff(const GrayImage &segmentation, const GrayImage &gt)
{
   int width      = imGetWidth(segmentation);
   int height     = imGetHeight(segmentation);
   float num_diff=0, num_pixels=width*height;
	for(int x=0;x<width;x++)
		for(int y=0;y<height;y++) {
			if(imRef(segmentation,x,y)!=imRef(gt,x,y)) num_diff+=1;
		//	if(imRef(gt,x,y)) num_pixels+=1;
		}
   float error_rate=num_diff/num_pixels;
   return error_rate;

}

void postProcessing(GrayImage &segmentation)
{
   int width      = imGetWidth(segmentation);
   int height     = imGetHeight(segmentation);
   int num_pixels=width*height;
   queue<int> active;
   vector<bool> markings(num_pixels,0);
   int seedID=getPixelID(width/2,height/2,width);
   markings[seedID]=1;
   int seedV=255; //imRef(segmentation,width/2,height/2);
   int backV=0; //imRef(segmentation,0,0);
   active.push(seedID);
   // run region grow for foreground from the seeds in the center
   while (!active.empty()) {
	   int pixelID=active.front();
	//   markings[pixelID]=1; // mark the pixel
	   active.pop();
	   int y=pixelID/width;
	   int x=pixelID-y*width;
	   if(x+1<width) { if(imRef(segmentation,x,y)==imRef(segmentation,x+1,y) && markings[pixelID+1]==0) active.push(pixelID+1); markings[pixelID+1]=1;}
	   if(x-1>=0) { if(imRef(segmentation,x,y)==imRef(segmentation,x-1,y) && markings[pixelID-1]==0) active.push(pixelID-1); markings[pixelID-1]=1;}
	   if(y+1<height) { if(imRef(segmentation,x,y)==imRef(segmentation,x,y+1) && markings[pixelID+width]==0) active.push(pixelID+width); markings[pixelID+width]=1;}
	   if(y-1>=0) { if(imRef(segmentation,x,y)==imRef(segmentation,x,y-1) && markings[pixelID-width]==0) active.push(pixelID-width); markings[pixelID-width]=1; }
	//   imRef(segmentation,x,y)=seedV;
   }

   for(int i=0; i<num_pixels; i++) {
	   if(markings[i]==0) {
		   int y=i/width; int x=i-y*width;
		   imRef(segmentation,x,y)=backV;
	   }
	   markings[i]=1;  //reset to forground for later use (region grow from background)
   }

   active.push(0);
   while (!active.empty()) {
	   int pixelID=active.front();
	   markings[pixelID]=0; // mark the pixel
	   active.pop();
	   int y=pixelID/width;
	   int x=pixelID-y*width;
	   if(x+1<width) { if(imRef(segmentation,x,y)==imRef(segmentation,x+1,y) && markings[pixelID+1]==1) active.push(pixelID+1); markings[pixelID+1]=0;}
	   if(x-1>=0) { if(imRef(segmentation,x,y)==imRef(segmentation,x-1,y) && markings[pixelID-1]==1) active.push(pixelID-1); markings[pixelID-1]=0;}
	   if(y+1<height) { if(imRef(segmentation,x,y)==imRef(segmentation,x,y+1) && markings[pixelID+width]==1) active.push(pixelID+width); markings[pixelID+width]=0;}
	   if(y-1>=0) { if(imRef(segmentation,x,y)==imRef(segmentation,x,y-1) && markings[pixelID-width]==1) active.push(pixelID-width); markings[pixelID-width]=0; }
   }

   for(int i=0; i<num_pixels; i++) {
	   if(markings[i]==1) {
		   int y=i/width; int x=i-y*width;
		   imRef(segmentation,x,y)=seedV;
	   }
   }

}

GrayImage saveToResult(const GrayImage &image,const GrayImage &gt, const GrayImage &segmentation, const GrayImage &dataterm, const GrayImage &smooth)
{
   int width      = imGetWidth(image);
   int height     = imGetHeight(image);
   GrayImage hd=(GrayImage) imNew(IMAGE_GRAY,width,height);
   int size_center=(width+height)/8; int size_border=size_center/3;
   for(int x=0;x<width;x++)
	   for(int y=0;y<height;y++) {
		   imRef(hd,x,y)=imRef(image,x,y);
           if(y==height/4 && x>=width/4 && x<=3*width/4) imRef(hd,x,y)=0;
		   if(y==3*height/4 && x>=width/4 && x<=3*width/4)  imRef(hd,x,y)=0;
		   if(x==width/4 && y>=height/4 && y<=3*height/4) imRef(hd,x,y)=0;
		   if(x==3*width/4 && y>=height/4 && y<=3*height/4) imRef(hd,x,y)=0;
	   }
		  
   GrayImage result=(GrayImage) imNew(IMAGE_GRAY,width*3+6,height*2+3);
   
   for(int x=0;x<width*3+6;x++)
	   for(int y=0;y<height*2+3;y++) imRef(result,x,y)=255;
   for(int x=0;x<width;x++)
	   for(int y=0;y<height;y++) {
		   imRef(result,x,y)=imRef(image,x,y);
		   imRef(result,x+width+3,y)=imRef(smooth,x,y);
		   imRef(result,x+2*width+6,y)=imRef(dataterm,x,y);
		   imRef(result,x,y+height+3)=imRef(hd,x,y);
		   imRef(result,x+width+3,y+height+3)=imRef(gt,x,y);
		   imRef(result,x+2*width+6,y+height+3)=imRef(segmentation,x,y);
	   }

	   return result;


}

void segment(const FloatImage &image, GrayImage &segmentation, GrayImage &dataterm, GrayImage &smooth, float lambda, int binNum)
{
   int width      = imGetWidth(image);
   int height     = imGetHeight(image);
   int num_pixels = (width)*(height);
   bin_num=binNum; 
   int size_fg=(width+height)/4;
   int size_border=5; int size_center=2;

   float variance = computeImageVariance(image); // image variance used in neigborhood edge computation
   
   BinaryImage mask=(BinaryImage) imNew(IMAGE_BINARY,width,height);
   setMask(mask,size_fg);  //set the initial foreground and background

   // some vectors 
   vector<float> pixelArray(num_pixels,0);
   vector<int> pixelBinID(num_pixels,0);
   vector<float> fgProb(bin_num,0);
   vector<float> bkProb(bin_num,0);
//   vector<int> pixelOrder(num_pixels);

 //  for(int i=0;i<num_pixels; i++) pixelOrder[i]=i;   //initialize pixel order

   GraphType *g=new GraphType(num_pixels,8*num_pixels);
   g->add_node(num_pixels);

   smooth=(GrayImage) imNew(IMAGE_GRAY,width,height);
   dataterm=(GrayImage) imNew(IMAGE_GRAY,width,height);
   segmentation=(GrayImage) imNew(IMAGE_GRAY,width,height);

   hardConstraint(g,width,height);
   computeAppearanceModel(image,mask,pixelArray,pixelBinID,fgProb,bkProb);
   computeTWeight(g,dataterm,pixelBinID,fgProb,bkProb);
   computeNWeight(g,image,smooth,variance,lambda);
//   setStarTlinksAllPixels(g,width,height,width/2,height/2,num_pixels);
   
   float flow = g -> maxflow();
   segmentationToImage(g,segmentation);
//   postProcessing(segmentation);

   delete g;
   imFree(mask);
   //imSave(smooth,"smooth.pgm");
   //imSave(dataterm,"dataterm.pgm");

}



int main(int argc, char **argv) {

	char file_name[200]="D:\\workspace\\project\\segment\\data\\image_name.txt";
	char gt_path[200]="D:\\workspace\\project\\segment\\data\\crop_segmentations\\";
	char image_path[200]="D:\\workspace\\project\\segment\\data\\crop_images\\";
	char result_path[200]="D:\\workspace\\project\\segment\\results\\result_no_starshape\\segmentation_pgm\\";
    FILE *fimage=fopen(file_name,"r");
	FloatImage error=(FloatImage) imNew(IMAGE_FLOAT,50,30); //[50][30]={0.0};
	int lambda; int binNum;
	int min_lambda,min_binNum; float minerr=100;

    char image_name[29][100];
	for(int i=0;i<29;i++) fscanf(fimage,"%s\n",image_name[i]);
	fclose(fimage);


		//STAR_HUGE=10;
	int l=18; int b=16;
	//for(int l=0; l<50; l++){
	//for(int b=0; b<30; b++) {
	lambda=l+1; binNum=(b+1)*10;
	imRef(error,l,b)=0;
	for(int i=0; i<29; i++) {
        char image_full_name[200], gt_full_name[200], output_name[200];
	    strcpy(image_full_name,image_path);
	    strcpy(gt_full_name,gt_path);
		strcpy(output_name,result_path);

//		fscanf(fimage,"%s\n",image_name);
		strcat(image_full_name,image_name[i]);
		strcat(image_full_name,"roi.pgm"); // get the full path + name for roi image
		strcat(gt_full_name,image_name[i]);
		strcat(gt_full_name,"seg.pgm"); // get the full path + name for roi image
		strcat(output_name,image_name[i]);
		strcat(output_name,"_rlt.pgm");

		GrayImage image=(GrayImage) imLoad(IMAGE_GRAY, image_full_name);
		GrayImage gt=(GrayImage) imLoad(IMAGE_GRAY,gt_full_name);
	//	printf("\n%s",gt_full_name);
        check_error( ( image == (GrayImage) NULL ), "Cannot load input image");
		check_error( ( gt == (GrayImage) NULL ), "Cannot load ground truth");
	    int width=imGetWidth(image); int height=imGetHeight(image);
	    FloatImage floatImage=(FloatImage) imNew(IMAGE_FLOAT,width,height);
	    gray2float(image,floatImage);
	    normalizeImage(floatImage,0,1);
	    bmin=0; bmax=1;
		GrayImage dataterm, smooth;
		GrayImage segmentation;
        segment(floatImage,segmentation,dataterm,smooth,lambda,binNum);
		imRef(error,l,b)+=computeDiff(segmentation,gt);
		GrayImage result=saveToResult(image,gt,segmentation,dataterm,smooth);
//		printf("\n%s\n",output_name);
		imSave(result,output_name);
		imFree(segmentation); imFree(floatImage); imFree(dataterm); imFree(smooth); imFree(image); imFree(gt);
	}
	imRef(error,l,b)=imRef(error,l,b)/29;
	if(imRef(error,l,b)<minerr) { minerr=imRef(error,l,b); min_lambda=lambda; min_binNum=binNum; }
	printf("lambda %d #bin %d error rate %3f\n",lambda,binNum,imRef(error,l,b));
	//}
	//}

	//printf("min error is %5f, best paramters are %d %d \n",minerr, min_lambda, min_binNum);
	//normalizeImage(error,0,255);
	//GrayImage error_image=(GrayImage) imNew(IMAGE_GRAY,50,30);
	//float2gray(error,error_image);
	//char *fileE="D:\\workspace\\project\\segment\\results\\error_starshape.pgm";
	//imSave(error_image,fileE);

	return 1;

}