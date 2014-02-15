#include <vector>
#include "graph.h"
#include "segment_3d.h"
using namespace std;

void addNextConstraintPixels2D(GraphType *g,int dims[],int seed[], int point[],
					         vector<int> &pixelID,vector<unsigned char> &markings)
							 
{
	int width=dims[1];
	int height=dims[0];
	int num_pixels=width*height;
	vector<int> pixels(num_pixels);
	int SeedX=seed[1]; int SeedY=seed[0];
	int x=point[1]; int y=point[0];
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
				node_id id1=pixelID[pixels[i-1]];
				node_id id2=pixelID[pixels[i]];
				g->add_edge(id1,id2,STAR_HUGE,0);
				if ( (y+x) == ((y+x)/DISPLAY)*DISPLAY ){

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
				node_id id1=pixelID[pixels[i]];
				node_id id2=pixelID[pixels[i-1]];
				g->add_edge(id1,id2,STAR_HUGE,0);

				if ( (y+x) == ((y+x)/DISPLAY)*DISPLAY ){
				}
			}
			else break;
		}
		
	}	
	
}

//---------------------------------------------------------
// chooses all image pixels (not just the border) in random order
void setStarTlinksAllPixels2D(GraphType *g,vector<int> &pixelID,int dims[],int seed[]){ //vector pixels is the pixelID for pixels on this slice

	int num_pixels=dims[0]*dims[1];
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

	for ( int i = 0; i < num_pixels; i++ )
		if ( markings[pixelOrder[i]] == 0 )
		{

			int nextY = pixelOrder[i]/width;
			int nextX = pixelOrder[i]-nextY*width;
			int point[2]={nextY,nextX};
			addNextConstraintPixels2D(g,dims,seed,point,pixelID,markings);
		}
	//int unmarked=0;
	//for(int i=0; i<num_pixels; i++) {
	//	if(markings[i]==0){ unmarked++; }
	//}
}



void addNextConstraintPixels3D(GraphType *g, int dims[],int seed[], int point[],
					         vector<int> &pixels,vector<unsigned char> &markings)
							 
{
	float ay,by,az,bz;
	int nextX;
	int numPix = 0;
	bool reverseOrder = false;
	float a[3]; float b[3];
	int index[3], index_else[2]; int max=0;
	for(int i=0;i<3;i++)
	{
		if(abs(seed[i]-point[i])>max) { max=abs(seed[i]-point[i]); index[0]=i; }
	}
	int j=1;
	for(int i=0;i<3;i++)
	{
		if(i!=index[0]) { index[j]=i; j++;}
	}

	// solve for line coefficients a and b first
	if ( seed[index[0]]!=point[index[0]] ){
		a[index[1]]=((float) (seed[index[1]]-point[index[1]]))/(seed[index[0]]-point[index[0]]);
		b[index[1]] = point[index[1]]-a[index[1]]*point[index[0]];
		a[index[2]]=((float) (seed[index[2]]-point[index[2]]))/(seed[index[0]]-point[index[0]]);
		b[index[2]] = point[index[2]]-a[index[2]]*point[index[0]];
	}

	// next get all pixels on the discretized line. 4 cases
	
		int start,end;
		if ( point[index[0]] > seed[index[0]]) {
			start = seed[index[0]]+1;
			end = point[index[0]];
		}
		else  {
			start = point[index[0]];
			end = seed[index[0]]-1;
			reverseOrder = true;
		}
			
		int next[3];
		for (  next[index[0]] = start; next[index[0]] <= end; next[index[0]]++ ){
			for(int i=1;i<3;i++) {
				if ( seed[index[i]] != point[index[i]]) 
					next[index[i]] = (int) (0.5+a[index[i]]*next[index[0]]+b[index[i]]);
				else next[index[i]] = seed[index[i]]; 
			}
				pixels[numPix] = getPixelID(next,dims);
				numPix++;
			}
		
	
	int DISPLAY = 1;
	
	if ( reverseOrder ){
		int i;
		for (  i = 1; i < numPix; i++ ){
			if ( markings[pixels[i-1]] == 0 ){
				markings[pixels[i-1]] = 1;
				g->add_edge(pixels[i-1],pixels[i],STAR_HUGE,0);
				//if ( (y+x) == ((y+x)/DISPLAY)*DISPLAY ){
				//	imRef(im,pixels[i-1]-width*(pixels[i-1]/width),pixels[i-1]/width) = color;
					//imSave(im,"connections.ppm");
				//}
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

			}
			else break;
		}
		
	}
	
	//saveArray(markings,width,height);
	//imSave(im,"connections.ppm");		
	
}

//---------------------------------------------------------
// chooses all image pixels (not just the border) in random order
void setStarTlinksAllPixels3D(GraphType *g,int dims[], int seed[]){

	int num_pixels=dims[0]*dims[1]*dims[2];
	vector<unsigned char> markings(num_pixels,0);
	vector<int> pixelOrder(num_pixels);
	vector<int> pixels(num_pixels);
	int numPixels = 0;

	node_id seed_id=getPixelID(seed,dims);
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

	int next[3];
	for ( int i = 0; i < num_pixels; i++ )
		if ( markings[pixelOrder[i]] == 0 )
		{

			int nextS = pixelOrder[i]/(dims[0]*dims[1]);
			int nextR = (pixelOrder[i]-nextS*dims[0]*dims[1])/dims[1];
			int nextC=pixelOrder[i]-nextS*dims[0]*dims[1]-nextR*dims[1];

			next[0]=nextR; next[1]=nextC; next[2]=nextS;
		//	RGB color = {rand()%256,rand()%256,rand()%256};
			addNextConstraintPixels3D(g,dims,seed,next,
			                    pixels,markings);
		}
	//int unmarked=0;
	//for(int i=0; i<num_pixels; i++) {
	//	if(markings[i]==0){ unmarked++; }
	//}
//	printf("unmarked pixels: %d\n",unmarked);

}