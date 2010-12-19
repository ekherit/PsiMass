#ifndef   SumPoints
#define   SumPoints
#include "TMath.h"


void SumPointsBySigma(Int_t npini,Int_t npfin,Int_t* UsePoints,Double_t* v,Double_t* s,Double_t* vn,Double_t* sn,bool ZeroOpt)
{
    Double_t w;
    Double_t a;
    for(Int_t i=0;i<npfin;i++)
    {
	w=0;
	a=0;
	for(Int_t j=0;j<npini;j++)
    	{
		if(UsePoints[j]==i)
		{
        		w+=(1./sq(s[j]));
			a+=(v[j]/sq(s[j]));
		}
	}
	vn[i]=a/w;
	if(ZeroOpt) sn[i]=0;
	else sn[i]=1./sqrt(w);

    }
};
void SumPointsByQuant(Int_t npini,Int_t npfin,Int_t* UsePoints,Double_t* v,Double_t* q,Double_t*s,Double_t* vn,Double_t* sn,bool ZeroOpt)
{
    Double_t w;
    Double_t a;
    Double_t norma;

    for(Int_t i=0;i<npfin;i++)
     {
	w=0;
	a=0;
	norma=0;
//	normw=0;
	for(Int_t j=0;j<npini;j++)
    	{
		if(UsePoints[j]==i)
		{
		    w+=sq(q[j]*s[j]);
//		    normw+=sq(q[j]);
		    norma+=q[j];
		    a+=(v[j]*q[j]);
		}
	}
	vn[i]=a/norma;
	if(ZeroOpt) sn[i]=0;
	else     sn[i]=sqrt(w/sq(norma));
    }
};
void SumPointsSimple(Int_t npini,Int_t npfin,Int_t* UsePoints,Double_t* v,Double_t* vn)
{    
    Double_t a;  
    for(Int_t i=0;i<npfin;i++)
    {       
	a=0;
	for(Int_t j=0;j<npini;j++)
    	{
		if(UsePoints[j]==i)
		{                  
		    a+=v[j];
		}
	}
	vn[i]=a;
    }
};

void SeparatePointsPart(Int_t npini,Int_t* npfin,Int_t* NpUse,Double_t* SortArray,Double_t Diff)
{

//    Double_t  ReperPoint=SortArray[0];
    Double_t  ReperPoint=0;
    Int_t     npf=0,nuse=1;
    NpUse[0]=npf;

    for(Int_t i=0;i<npini;i++)
    {
	if(ReperPoint<SortArray[i])
	{
	    ReperPoint=SortArray[i];
	}
    }

    for(Int_t i=1;i<npini;i++)
    {
	if(fabs(ReperPoint-SortArray[i])<Diff)
	{
	      NpUse[i]=npf;
	      ReperPoint=(ReperPoint*nuse+SortArray[i])/(nuse+1);
	      nuse++;

	}
	else
	{
	    npf++;
	    NpUse[i]=npf;
	    ReperPoint=SortArray[i];
	    nuse=1;
	}
    }
    *npfin=(npf+1);//for dimension of arrays
};

void SeparatePointsPartNew(Int_t npini,Double_t* q,Int_t* npfin,Int_t* NpUse,Double_t* SortArray,Double_t Diff)
{

//    Double_t  ReperPoint=SortArray[0];
    Double_t  ReperPoint=0;
    Int_t     npf=0;
//    Int_t     nuse=1;
    Double_t norma;
    NpUse[0]=npf;

    /*
    for(Int_t i=0;i<npini;i++)
    {
	if(ReperPoint<SortArray[i])
	{
	    ReperPoint=SortArray[i];
	}
    }*/
    norma=q[0];
    ReperPoint=SortArray[0];
    for(Int_t i=1;i<npini;i++)
    {
	if(fabs(ReperPoint-SortArray[i])<Diff)
	{
	      NpUse[i]=npf;
	      ReperPoint=(ReperPoint*norma+SortArray[i]*q[i])/(norma+q[i]);
	      norma+=q[i];
	}
	else
	{
	    npf++;
	    NpUse[i]=npf;
	    ReperPoint=SortArray[i];
	    norma=q[i];
	}
    }
    *npfin=(npf+1);//for dimension of arrays
};

void SeparatePointsPartPsiPP(Int_t npini,Double_t* q,Int_t* npfin,Int_t* NpUse,Double_t* SortArray,Double_t Diff1,Double_t Diff2,Double_t EEdge)
{
    Double_t  ReperPoint=0;
    Double_t  Diff=0;
    Int_t     npf=0;
    Double_t norma;
    NpUse[0]=npf;
    norma=q[0];
    ReperPoint=SortArray[0];
    
    for(Int_t i=1;i<npini;i++)
    {
      if(ReperPoint<EEdge) Diff=Diff1;
      else  Diff=Diff2;
	if(fabs(ReperPoint-SortArray[i])<Diff)
	{
	      NpUse[i]=npf;
	      ReperPoint=(ReperPoint*norma+SortArray[i]*q[i])/(norma+q[i]);
	      norma+=q[i];
	}
	else
	{
	    npf++;
	    NpUse[i]=npf;
	    ReperPoint=SortArray[i];
	    norma=q[i];
	}
    }
    *npfin=(npf+1);//for dimension of arrays
};

Double_t AproxCorrBB(Int_t dim,Double_t place, Double_t** Array)
{
      Double_t Corr=1;
      for(Int_t ic=0;ic<(dim-1);ic++)
      {
      	if(Array[ic+1][0]>place&&Array[ic][0]<=place)
      	{
      		Corr=(1.+
		Array[ic][1]+
		(Array[ic+1][1]-Array[ic][1])
		/(Array[ic+1][0]-Array[ic][0])
		*(place-Array[ic][0]));
      	}
      }
      return Corr;
}

Double_t AproxCorr(Double_t add,Int_t dim,Double_t place,Double_t sigma,Double_t Array[][2])
{
      Double_t Corr=1;
      for(Int_t ic=0;ic<(dim-1);ic++)
      {
      	if(Array[ic+1][0]>place&&Array[ic][0]<=place)
      	{
      		Corr=add+(Array[ic][1]+
		(Array[ic+1][1]-Array[ic][1])
		/(Array[ic+1][0]-Array[ic][0])
		*(place-Array[ic][0]))/sigma;
      	}
      }
      return Corr;
}
#endif
