#include<vector>
using namespace std;
//kernel density estimation 
// hist: histogram of sample
// n:   number of data
// bin_value: corresponding value of each bin
// sigma: parameter that controls density estimation
// bin with zero value will not be counted
// 
void parzenwin(vector<float> &hist,vector<float> &bin_value, vector<float> &flag,float sigma)
{
	float SMALLCONST=1.0e-20;
	vector<float> temp(hist);
	int bin_num=hist.size();
	float sigma2=sigma*sigma;
	float v; float n;
	float min=100;
	for(int i=0;i<bin_num;i++)
	{
		hist[i]=0;
		if(flag[i]==0) continue;
		n=0;
		for(int j=0;j<bin_num;j++)
		{
			if(flag[i]==0)continue;
			n+=temp[j];
			v=bin_value[i]-bin_value[j];
			v=v*v;
			hist[i]+=temp[j]*exp(-0.5*v/sigma2);
		}
		hist[i]=hist[i]/2.5067;
		hist[i]=hist[i]/sigma;
		hist[i]=hist[i]/n;
	}
    for(int i=0;i<bin_num;i++)
	{
		hist[i]=hist[i]+SMALLCONST;
	}

}


float computeSigmaFromHist(vector<float> &hist,vector<float> &bin_value)
{
	int bin_num=hist.size();
	float sum=0;
	float num=0;
	for(int i=0;i<bin_num;i++)
	{
		sum+=hist[i]*bin_value[i];
		num+=hist[i];
	}
	if(num==1) return 1;
	float mean=sum/num;
	float std=0;
	for(int i=0;i<bin_num;i++)
	{
		if(hist[i]==0) continue;
		std+=hist[i]*(bin_value[i]-mean)*(bin_value[i]-mean);
	}
	return sqrt(std/num);

}

void testZero(vector<vector<float> > &Prob,vector<vector<float> > &flagHist)
{
	int bin_num=Prob[0].size();
	int slices=Prob.size();
//	printf("foreground prob:\n");
	for(int s=1;s<slices-1;s++)
	{
		for(int j=0;j<bin_num;j++)
		{
			if(flagHist[s][j]==0) continue;
			if(Prob[s][j]==0) printf("fg slice #bin %d %d\n",s,j);
		}
	}

}

float computeSum(vector<float> &arry)
{
	int num=arry.size();
	float sum=0;
	for(int i=0;i<num;i++)
		sum+=arry[i];

	return sum;

}

void kernelDensityEstLocal(vector<vector<float> > &fgProb,vector<vector<float> > &bkProb, vector<float> &bin_value,vector<float> &flagHist,float ratio)
{
	int slices=fgProb.size();
	int bin_num=fgProb[0].size();
	for(int s=0;s<slices;s++)
	{
		float fgnum=computeSum(fgProb[s]);
		if(fgnum>0)
		{
		float fgSigma=computeSigmaFromHist(fgProb[s],bin_value);
		fgSigma=1.06*fgSigma/pow(fgnum,0.2f);
		fgSigma=ratio*fgSigma;
	//	printf("fgsigma %d %f\n",s,fgSigma);
		parzenwin(fgProb[s],bin_value,flagHist,fgSigma);
		}
		float bknum=computeSum(bkProb[s]);
		if(bknum>0)
		{
		float bkSigma=computeSigmaFromHist(bkProb[s],bin_value);
		bkSigma=1.06*bkSigma/pow(bknum,0.2f);
		bkSigma=ratio*bkSigma;
	//	printf("bksigma %d %f\n",s,bkSigma);
		parzenwin(bkProb[s],bin_value,flagHist,bkSigma);
		}
	}
}

void kernelDensityEst(vector<float> &hist, vector<float> &bin_value, vector<float> &flagHist,float ratio)
{
	int bin_num=bin_value.size();
	float fgnum=computeSum(hist);
	if(fgnum==0) return;

		float fgSigma=computeSigmaFromHist(hist,bin_value);
		fgSigma=1.06*fgSigma/pow(fgnum,0.2f);
		fgSigma=ratio*fgSigma;
		if(fgSigma==0 || fgnum==1) fgSigma=1;
		parzenwin(hist,bin_value,flagHist,fgSigma);

	float sum=0;
	for(int i=0;i<bin_num;i++)
        sum+=hist[i];
	for(int i=0;i<bin_num;i++)
		hist[i]=hist[i]/sum;
}



//----------------------------------------------------------------------------------------
void parzenwin(vector<float> &hist,vector<float> &bin_value, float sigma)
{
	float SMALLCONST=1.0e-20;
	vector<float> temp(hist);
	int bin_num=hist.size();
	float sigma2=sigma*sigma;
	float v; float n;
	float min=100;
	for(int i=0;i<bin_num;i++)
	{
		hist[i]=0;
	//	if(flag[i]==0) continue;
		n=0;
		for(int j=0;j<bin_num;j++)
		{
		//	if(flag[i]==0)continue;
			n+=temp[j];
			v=bin_value[i]-bin_value[j];
			v=v*v;
			hist[i]+=temp[j]*exp(-0.5*v/sigma2);
		}
		hist[i]=hist[i]/2.5067;
		hist[i]=hist[i]/sigma;
		hist[i]=hist[i]/n;
		if(hist[i]>0 && hist[i]<min) min=hist[i];
	}
    for(int i=0;i<bin_num;i++)
	{
	//	if(hist[i]==0) hist[i]=min;
		hist[i]=hist[i]+SMALLCONST;
	}

}





void kernelDensityEstLocal(vector<vector<float> > &fgProb,vector<vector<float> > &bkProb, vector<float> &bin_value,float ratio)
{
	int slices=fgProb.size();
	int bin_num=fgProb[0].size();
	for(int s=0;s<slices;s++)
	{
		float fgnum=computeSum(fgProb[s]);
		if(fgnum>0)
		{
		float fgSigma=computeSigmaFromHist(fgProb[s],bin_value);
		fgSigma=1.06*fgSigma/pow(fgnum,0.2f);
		fgSigma=ratio*fgSigma;
	//	printf("fgsigma %d %f\n",s,fgSigma);
		parzenwin(fgProb[s],bin_value,fgSigma);
		}
		float bknum=computeSum(bkProb[s]);
		if(bknum>0)
		{
		float bkSigma=computeSigmaFromHist(bkProb[s],bin_value);
		bkSigma=1.06*bkSigma/pow(bknum,0.2f);
		bkSigma=ratio*bkSigma;
	//	printf("bksigma %d %f\n",s,bkSigma);
		parzenwin(bkProb[s],bin_value,bkSigma);
		}
	}
}

void kernelDensityEst(vector<float> &hist, vector<float> &bin_value, float ratio)
{
	int bin_num=bin_value.size();
	float fgnum=computeSum(hist);
	if(fgnum==0) return;

		float fgSigma=computeSigmaFromHist(hist,bin_value);
		fgSigma=1.06*fgSigma/pow(fgnum,0.2f);
		fgSigma=ratio*fgSigma;
		if(fgSigma==0 || fgnum==1) fgSigma=1;
		parzenwin(hist,bin_value,fgSigma);

	float sum=0;
	for(int i=0;i<bin_num;i++)
        sum+=hist[i];
	for(int i=0;i<bin_num;i++)
		hist[i]=hist[i]/sum;
}