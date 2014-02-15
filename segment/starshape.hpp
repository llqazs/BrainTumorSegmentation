#include <vector>
#include "graph.h"
#include "segment_3d.hpp"
using namespace std;

// discretize lines in an image,  x2-x1>abs(y2-y1)
void bresenham(int x1,int y1, int x2, int y2, vector<int> &py)
{
	 if((x2-x1)<abs(y2-y1)) { printf("x2-x1 should be larger than abs(y2-y1)\n"); exit(0); }
	 int numPix=x2-x1+1;
	 if(y2==y1)
	 {
		 for(int i=0;i<numPix;i++)
			 py[i]=y1;
		 return;
	 }

     int deltax=x2-x1;
     int deltay =abs(y2-y1);
     int error=deltax/2;
     int ystep;
     int y=y1;
	 if(y1<y2) ystep=1;
	 else ystep=-1;
	 for(int i=0;i<numPix;i++)
	 {
		 py[i]=y;
         error=error-deltay;
		 if(error<0)
		 {
			 y=y+ystep;
		     error=error+deltax;
	     }
	 }
}

void addNextConstraintPixels2D(GraphType *g,int dims[],int seed[], int point[],
					         vector<int> &pixelID,vector<unsigned char> &markings)
							 
{
	int width=dims[1];
	int height=dims[0];
	bool reverseOrder = false;

	int index[2]; int max=0;
	for(int i=0;i<2;i++)
	{
		if(abs(seed[i]-point[i])>max) { max=abs(seed[i]-point[i]); index[0]=i; }
	}
	int j=1;
	for(int i=0;i<2;i++)
	{
		if(i!=index[0]) { index[j]=i;}
	}
	
	int x1,x2,y1,y2;
	if(point[index[0]]<seed[index[0]])
	{
		reverseOrder=true;
		x1=point[index[0]];
		x2=seed[index[0]];
		y1=point[index[1]];
		y2=seed[index[1]];
	}
	else {
		x1=seed[index[0]];
		x2=point[index[0]];
		y1=seed[index[1]];
		y2=point[index[1]];
	}

	int numPix=x2-x1+1;
	vector<vector<int> > pt(2,vector<int>(numPix,0));
	vector<int> pixels(numPix);  //pixel ID for current dimension

	for(int i=0;i<numPix;i++)
	{
		pt[index[0]][i]=x1+i;
	}
	bresenham(x1,y1,x2,y2,pt[index[1]]);

	for(int i=0;i<numPix;i++)
	{
		pixels[i]=pt[0][i]*width+pt[1][i];
	}
	
	if ( reverseOrder ){
		int i;
		for (  i = 1; i < numPix; i++ ){
		//	if ( markings[pixels[i-1]] == 0 ){
				markings[pixels[i-1]] = 1;
				node_id id1=pixelID[pixels[i-1]];
				node_id id2=pixelID[pixels[i]];
				g->add_edge(id1,id2,STAR_HUGE,0);
		//	}
		//	else break; // break if current pixel has been added to the constraint
		}
		
	}
	else {
		int i;
		for (  i  = numPix-1; i > 0; i-- ){
	//		if ( markings[pixels[i]] == 0 ){
				markings[pixels[i]] = 1;
				node_id id1=pixelID[pixels[i]];
				node_id id2=pixelID[pixels[i-1]];
				g->add_edge(id1,id2,STAR_HUGE,0);
		//	}
		//	else break; // break if current pixel has been added to the constraint
		}
		
	}	
	
}

//---------------------------------------------------------
// chooses all image pixels (not just the border) in random order
void setStarTlinksAllPixels2D(GraphType *g,vector<int> &pixelID,int dims[3],int seed[3]){ //vector pixels is the pixelID for pixels on this slice

	int num_pixels=pixelID.size();
	vector<unsigned char> markings(num_pixels,0);
	vector<int> pixelOrder(num_pixels);

	int numPixels = 0;
	int seedX=seed[1], seedY=seed[0];
	int width=dims[1], height=dims[0];
	markings[seedX+seedY*width] = 1;


	for ( int i = 0; i < num_pixels; i++ )
		pixelOrder[i] = i;

	for ( int i = 0; i < 3*num_pixels; i++ ){
		int first  = rand()*rand()%num_pixels;
		int second = rand()*rand()%num_pixels;
		int temp = pixelOrder[first];
		pixelOrder[first]  = pixelOrder[second];
		pixelOrder[second] = temp;
	}

	int unmarked=0;
	for ( int i = 0; i < num_pixels; i++ )
		if ( markings[pixelOrder[i]] == 0 )
		{

			int nextY = pixelOrder[i]/width;
			int nextX = pixelOrder[i]-nextY*width;
			int point[2]={nextY,nextX};
			addNextConstraintPixels2D(g,dims,seed,point,pixelID,markings);
	//		if(markings[pixelOrder[i]]==0) printf("%d %d %d\n",i,nextY,nextX);
			if(markings[pixelOrder[i]]==0){ unmarked++; }
		}
	 
	if(unmarked>0) printf("2D starshape unmarked pixels: %d\n",unmarked);

}



void addNextConstraintPixels3D(GraphType *g, int dims[3],int seed[3], int point[3],
					         vector<int> &pixelID,vector<unsigned char> &markings)
							 
{
	bool reverseOrder = false;
	int index[3]; int max=0;
	for(int i=0;i<3;i++)
	{
		if(abs(seed[i]-point[i])>max) { max=abs(seed[i]-point[i]); index[0]=i; }
	}
	int j=1;
	for(int i=0;i<3;i++)
	{
		if(i!=index[0]) { index[j]=i; j++;}
	}
	
	int x1,y1,x2,y2;
	// ensure that x2>=x1
	x1=seed[index[0]]; x2=point[index[0]];
	if(point[index[0]]<seed[index[0]]) {
		reverseOrder=true;
		x1=point[index[0]];
		x2=seed[index[0]];
	}

	int numPix=x2-x1+1;
	vector<vector<int> > pt(3,numPix);
	//initial values for direction index[0]
	for(int i=0;i<numPix;i++)
		pt[index[0]][i]=x1+i;

	// get coordinates for other directions
	for(int i=1;i<3;i++)
	{
	    
		if(reverseOrder) {
			y1=point[index[i]];
			y2=seed[index[i]];
		}
		else{
			y1=seed[index[i]];
			y2=point[index[i]];
		}

		bresenham(x1,y1,x2,y2,pt[index[i]]);
	}

	// get current pixel ID 
	vector<int> pixels(numPix);
	for(int i=0;i<numPix;i++)
	{
		pixels[i]=pt[0][i]*dims[1]*dims[2]+pt[1][i]*dims[2]+pt[2][i]; 
	}

	
	if ( reverseOrder ){
		int i;
		for (  i = 0; i < numPix-1; i++ ){
	//		if ( markings[pixels[i]] == 0 ){
				markings[pixels[i]] = 1;
				node_id id1=pixelID[pixels[i]];
				node_id id2=pixelID[pixels[i+1]];
				g->add_edge(id1,id2,STAR_HUGE,0);
		//	}
		//	else break; // break if current pixel has been added to the constraint
		}
	}
	
	else {
		int i;
		for (  i  = numPix-1; i > 0; i-- ){
		//	if ( markings[pixels[i]] == 0 ){
				markings[pixels[i]] = 1;
				node_id id1=pixelID[pixels[i]];
				node_id id2=pixelID[pixels[i-1]];
				g->add_edge(id1,id2,STAR_HUGE,0);

	//		}
		//	else break; // break if current pixel has been added to the constraint
		}
		
	}
	
}

//---------------------------------------------------------
// chooses all image pixels (not just the border) in random order
void setStarTlinksAllPixels3D(GraphType *g,int dims[3], int seed[3],vector<int> &pixelID){

	int num_pixels=dims[0]*dims[1]*dims[2];
	vector<unsigned char> markings(num_pixels,0);
	vector<int> pixelOrder(num_pixels);
	int numPixels = 0;

	node_id seed_id=seed[0]*dims[1]*dims[2]+seed[1]*dims[2]+seed[2];
	markings[seed_id] = 1;

	for ( int i = 0; i < num_pixels; i++ )
		pixelOrder[i] = i;

	for ( int i = 0; i < 3*num_pixels; i++ ){
		int first  = rand()*rand()%num_pixels;
		int second = rand()*rand()%num_pixels;
		int temp = pixelOrder[first];
		pixelOrder[first]  = pixelOrder[second];
		pixelOrder[second] = temp;
	}

	int next[3]; int unmarked=0;
	for ( int i = 0; i < num_pixels; i++ )
		if ( markings[pixelOrder[i]] == 0 )
		{

			int nextS = pixelOrder[i]/(dims[1]*dims[2]);
			int nextR = (pixelOrder[i]-nextS*dims[1]*dims[2])/dims[2];
			int nextC=pixelOrder[i]-nextS*dims[1]*dims[2]-nextR*dims[2];

			next[0]=nextS; next[1]=nextR; next[2]=nextC;
		//	RGB color = {rand()%256,rand()%256,rand()%256};
			addNextConstraintPixels3D(g,dims,seed,next,
			                    pixelID,markings);
			if(markings[pixelOrder[i]]==0) unmarked++; //printf("%d %d %d %d\n",i,nextR,nextC,nextS);
		}
	if(unmarked>0) printf("3D starshape unmarked pixels: %d\n",unmarked);
	//int unmarked=0;
 //   for(int r=0;r<dims[0];r++)
	//	for(int c=0;c<dims[1];c++)
	//		for(int s=0;s<dims[2];s++)
	//		{
	//			int pos[3]={r,c,s};
	//			node_id id=getPixelID(pos,dims);
	//			if (markings[id]==0) {
	//				printf("%d %d %d %d\n",id,r,c,s);
	//				unmarked++;
	//			}
	//		}

	//printf("unmarked pixels: %d\n",unmarked);

}