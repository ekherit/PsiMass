#include<iostream.h>
#include <math.h>
#include <TFormula.h>
#include <TF1.h>
#include <TGraphErrors.h>

#include"FitTools/MathLibrary.h"

const Double_t sqrtln4=1.177410022515475;
Double_t AGauss(Double_t *x, Double_t *par)
{
  Double_t A=par[0];// const
  Double_t mean=par[1];// mean
  Double_t sigma=par[2];//FWHM/2.35
  Double_t W=par[3];// assymetry
  Double_t y=x[0];
  if (TMath::Abs(W)<=1e-6) W=TMath::Sign(1.e-6,W);
  if (sigma<=0.1) return 0;
  Double_t e=1+TMath::SinH(W*sqrtln4)/sqrtln4*(y-mean)/sigma;
  if (e<1.e-07) return 0;
  Double_t delta=TMath::Log(e);
  Double_t result=A*TMath::Exp(-0.5*(delta*delta/W/W+W*W));
  return result;
}
Double_t TMFCut(Double_t* x, Double_t* par)
{
	Double_t E=x[0];
	Double_t TMup=0;
	TMup=par[0]+(1.-par[0])*exp(-E/par[1]);
	return TMup;
}
Double_t BetheBloch(Double_t* x,Double_t* par){
    Double_t  y=pow(10,x[0]); //beta*gamma
    Double_t  bb;
    Double_t  beta,gamma,beta2;
    gamma=sqrt(1.+sq(y));
    beta=y/gamma;
    beta2= pow(beta,par[4]);
    bb=par[0]*(par[1]-beta2-log(par[2]+pow(y,-par[3])))/beta2;
    return bb;
};
Double_t BetheBloch2(Double_t* x,Double_t* par){
    Double_t  y=pow(10,x[0]); //beta*gamma
    Double_t  bb;
    Double_t  beta,gamma,beta2;
    gamma=sqrt(1.+sq(y));
    beta=y/gamma;
    beta2= pow(beta,par[3]);
    bb=par[0]*(par[1]-2*beta2-log(pow(y,-par[2])))/beta2;
    return bb;
};

Double_t BetheBloch_new(Double_t* x,Double_t* par){
    Double_t  y=pow(10,x[0]); //beta*gamma
    Double_t  bb;
    Double_t  beta,gamma,beta2;
    gamma=sqrt(1.+sq(y));
    beta=y/gamma;
    beta2=pow(beta,par[4]);
    bb=par[0]*(par[1]-beta2-par[6]*log((1+par[2]*pow(y,-par[3]))/(1+par[5]*pow(y,-par[3]))))/beta2;
    return bb;
};
Double_t BetheBloch_new1(Double_t* x,Double_t* par){
    Double_t  y=pow(10,x[0]); //beta*gamma
    Double_t  bb;
    Double_t  beta,gamma,beta2,betaP;
    gamma=sqrt(1.+sq(y));
    beta=y/gamma;
    beta2=pow(beta,par[4]+par[7]*pow(beta,par[8]));
    //    betaP=pow(beta,par[4]-fabs(par[7]));    
    bb=par[0]*(par[1]-beta2-par[6]*log((1+par[2]*pow(y,-par[3]))/(1+par[5]*pow(y,-par[3]))))/beta2;
    return bb;
};

Double_t dEdxResolutionN(Double_t* x,Double_t* p){
    Double_t  y=x[0]/42.; 
    Double_t  R=p[0]*pow(y,p[1]);
    return R;
};

Double_t dEdxResolutionL(Double_t* x,Double_t* p){
    Double_t  y=x[0]/10.; 
    Double_t  R=p[0]*pow(y,p[1]);
    return R;
};

Double_t dEdxResolutionM(Double_t* x,Double_t* par){
    Double_t Num=x[0];
    Double_t L=x[1];
    Double_t BB=x[2];
    Double_t pt=x[3];            
    Double_t  R=0;
    R=par[0]*pow(Num/42.,par[1])*pow(L/21.,par[2])*(1.+par[3]*sq(BB-par[4]))*(1.+par[5]*sq(pt-par[6]));
    return R;
};

/*
Double_t DistExp(Double_t* x,Double_t* par)
{
  Double_t  d=x[0];
  Double_t  dmult=1;
  
  if(d<=par[3])
  {    
    dmult=par[0]*(1+par[2]*(sq(d-par[1])-sq(par[1])));
  }
  else {
    //dmult=par[0]*(1-pow(par[1],2)+pow(par[2]-par[1],2));
    dmult=par[0]*(1+par[2]*(sq(par[3]-par[1])-sq(par[1])));
    dmult*=(exp(-(d-par[3])/par[4])*(1+par[5]*(d-par[3])));
  }
  return dmult;
}*/
Double_t DistExp(Double_t* x,Double_t* p)
{
  Double_t  d=x[0];
  Double_t  dmult=0;  
  if(d<p[1]){
    dmult=1+p[0]*d;
  }
  if(d>=p[1]) dmult=1+p[0]*p[1]+p[2]*(d-p[1]);
  if(d>=p[3]) dmult=1+p[0]*p[1]+p[2]*(p[3]-p[1]);
  if(d>=p[4]) dmult=1+p[0]*p[1]+p[2]*(p[3]-p[1])+p[5]*(d-p[4]);  
  return dmult;
}
Double_t DistF(Double_t* x,Double_t* p)
{
  Double_t  d=x[0];
  Double_t  f=0;
  if(d<p[0])
  {
    f=p[3]+p[2]*sq(d-p[1]);
  }
  else if(d>=p[0]&&d<p[4])
  {
    f=p[3]+p[2]*sq(p[0]-p[1])+p[5]*(d-p[0]); //+p[6]*sq(d-p[0]);
  }
  else if(d>=p[4]&&d<4)
  {
    f=p[3]+p[2]*sq(p[0]-p[1])+p[5]*(p[4]-p[0])+p[6]*(d-p[4])+p[7]*sq(d-p[4])+
      p[8]*cub(d-p[4]);
  }
  return f;
}

Double_t DistFtest(Double_t* x,Double_t* p)
{
  Double_t  d=x[0];
  Double_t  f=0;
  if(d<p[4])
  {
    f=p[0]+p[1]*d+p[2]*d*d+p[3]*cub(d);
  }
  else if(d>=p[4]&&d<4)
  {
    f=p[0]+p[1]*p[4]+p[2]*p[4]*p[4]+p[3]*cub(p[4])+p[5]*(d-p[4])+p[6]*sq(d-p[4])+p[7]*cub(d-p[4])+p[8]*sq4(d-p[4]);    //+p[9]*cub(d)*sq(d);
  }
  return f;
}

Double_t DistFN(Double_t* x,Double_t* p)
{
  Double_t  d=x[0];
  Double_t  f=0;
  if(d<p[0])
  {
    f=p[2]+p[1]*d;
  }
  else if(d>=p[0]&&d<p[3])
  {
    f=p[2]+p[0]*p[1]+p[4]*(d-p[0])+p[5]*sq(d-p[0]);
  }
  else if(d>=p[3]&&d<4)
  {
    f=p[2]+p[0]*p[1]+p[4]*(p[3]-p[0])+p[5]*sq(p[3]-p[0])+p[6]*sq(d-p[3])+
      p[7]*cub(d-p[3]);
  }
  return f;
}

Double_t DistF9(Double_t* x,Double_t* p)
{
  Double_t  t=x[0];
  Double_t  f=0;
  
  f=p[0]*T0()+p[1]*T1(t)+p[2]*T2(t)+p[3]*T3(t)+
    p[4]*T4(t)+p[5]*T5(t)+p[6]*T6(t)+p[7]*T7(t)+p[8]*T8(8);
  return f;
}
Double_t ThetaF(Double_t* x,Double_t* p)
{
  Double_t  d=x[0];
  Double_t  dmult=0;  
  if(d<p[1])
  {
    dmult=p[0]+p[2]*sq(d);
  }
  else if(d>=p[1]&&d<p[6])
  {
    dmult+=(p[0]+p[2]*sq(p[1]));
    dmult+=p[3]*(d-p[4])+p[5]*cub(d-p[4])-p[3]*(p[1]-p[4])-p[5]*cub(p[1]-p[4]);
  }    
  else if(d>=p[6])
  {
    dmult+=(p[0]+p[2]*sq(p[1]));
    dmult+=p[3]*(p[6]-p[4])+p[5]*cub(p[6]-p[4])-p[3]*(p[1]-p[4])-p[5]*cub(p[1]-p[4]);  
    dmult+=p[7]*(1./(p[8]-d)-1./(p[8]-p[6]));
  }    
  return dmult;
}
Double_t TimeF(Double_t* x,Double_t* p)
{  
  Double_t  d=x[0];
  Double_t  f=0;
  if(d<p[2])
  {
    f=p[0];
  }
  else if(d>=p[2]&&d<p[4])
  {
    f=p[0]+p[1]*(d-p[2])+p[3]*(cub(d-p[6])-cub(p[2]-p[6]));
  }
  else if(d>=p[4])
  {
    f=p[0]+p[1]*(p[4]-p[2])+p[3]*(cub(p[4]-p[6])-cub(p[2]-p[6]))+p[5]*(d-p[4]);
  }
  return f;  
 
}

Double_t TimeFtest(Double_t* x,Double_t* p)
{  
  Double_t  d=x[0];
  Double_t  f=0;
  /*    if(d<p[2])
  {
    f=p[0]+p[1]*d+p[6]*sq(d);
  }
  else if(d>=p[2])
  {
    f=p[0]+p[1]*p[2]+p[6]*sq(p[2])-p[3]*(d-p[2])+p[4]*atan(p[5]*(d-p[2]));
    }*/

  //  f=p[0]+p[1]*p[4]+p[2]*sq(p[4])+p[3]*cub(p[4])+p[5]*(d-p[4])+p[6]*sq(d-p[4])+p[7]*cub(d-p[4])+p[8]*sq4(d-p[4]);
  
    if(d<=p[2])
    {
      f=p[0]+p[1]*sq(d);
    }
  if((d>p[2])&&(d<p[5]))
    {
      f=p[0]+p[1]*sq(p[2])+p[3]*(d-p[2])+p[4]*sq(d-p[2]);
    }
  else if(d>=p[5])
    {
      f=p[0]+p[1]*sq(p[2])+p[3]*(p[5]-p[2])-p[6]*(d-p[5])+p[7]*atan(p[6]*(d-p[5]))+p[4]*sq(p[5]-p[2])+p[8]*sq(d-p[5]);
      }
  return f;  
}

Double_t Pol3(Double_t* x,Double_t* p)
{
  Double_t t=x[0];
 Double_t pol=0;
 pol=p[0]+p[1]*t+p[2]*sq(t)+p[3]*cub(t);
 return pol;
}
Double_t PolN(Double_t* x, Double_t* p,Int_t N)
{
  Double_t pol=0;
  Double_t t = x[0];
  if(N==1) pol=p[0]+p[1]*t;
  if(N==2) pol=p[0]+p[1]*t+p[2]*sq(t);
  if(N==3) pol=p[0]+p[1]*t+p[2]*sq(t)+p[3]*cub(t);
  if(N==4) pol=p[0]+p[1]*t+p[2]*sq(t)+p[3]*cub(t)+p[4]*sq4(t);
  if(N==8) pol=p[0]+p[1]*t+p[2]*sq(t)+p[3]*cub(t)+p[4]*sq4(t)+p[5]*sq5(t)+p[6]*sq6(t)+p[7]*sq7(t)+p[8]*sq8(t);
  return pol;
}

void FitPoln(TF1** poln,Double_t* param, TGraphErrors* dEvsDist,Double_t* x,Double_t* y, Double_t* sigmay,Int_t npoints,Int_t* degr,Int_t NumPar)
{
  Double_t chisq1=0;
  Double_t chisq2=0;
  Double_t chisq3=0;
  Double_t delta1,delta2,delta3;
  Double_t* par1=new Double_t [NumPar];
  Double_t* par2=new Double_t [NumPar];
  Double_t* par3=new Double_t [NumPar];
  poln[0]->SetParameters(param);
  poln[1]->SetParameters(param);
  poln[2]->SetParameters(param);
  dEvsDist->Fit(poln[0],"LIVMR");
  poln[0]->GetParameters(par1);
  dEvsDist->Fit(poln[1],"LIVMR");
  poln[1]->GetParameters(par2);
  dEvsDist->Fit(poln[2],"LIVMR");
  poln[2]->GetParameters(par3);
    for(int j=0;j<npoints;j++){
      delta1=(y[j]-PolN(&x[j],par1,degr[0]))/sigmay[j];
      chisq1+=delta1*delta1;

      delta2=(y[j]-PolN(&x[j],par2,degr[1]))/sigmay[j];
      chisq2+=delta2*delta2;

      delta3=(y[j]-PolN(&x[j],par3,degr[2]))/sigmay[j];
      chisq3+=delta3*delta3;
    }
    if(chisq1<=chisq2 && chisq1<=chisq3){
      dEvsDist->Fit(poln[0],"LIVMR");
      for(int j=0;j<NumPar-1;j++) param[j+1]=par1[j];
      param[0]=1;
    }
    if(chisq2<=chisq1 && chisq2<=chisq3){
      dEvsDist->Fit(poln[1],"LIVMR");
      for(int j=0;j<NumPar-1;j++) param[j+1]=par2[j];
      param[0]=4;
    }
    if(chisq3<=chisq1 && chisq3<=chisq2){
      dEvsDist->Fit(poln[2],"LIVMR");
      for(int j=0;j<NumPar-1;j++) param[j+1]=par3[j];
      param[0]=2;
    }

    delete [] par1;
    delete [] par2;
    delete [] par3; 
}


Double_t ThetaF9(Double_t* x,Double_t* p)
{
  Double_t  t=x[0];
  Double_t  dmult=0;  
  //  p[5]=0;
  // p[6]=0;
  //  p[7]=0;
  //  p[8]=0;
  //  dmult=1+p[1]/sqrt(sq(p[0])-sq(d))+p[2]*(d-0.5)+p[3]*sq4(d-0.5);
  //  dmult=(1+p[1]/sqrt(sq(p[0])-sq(d)))*(1+p[2]*(d-0.5)+p[3]*cub(d-0.5));
  dmult=p[0]*T0()+p[1]*T1(t)+p[2]*T2(t)+p[3]*T3(t)+
    p[4]*T4(t)+p[5]*T5(t)+p[6]*T6(t)+p[7]*T7(t)+p[8]*T8(t);

  return dmult;
}

Double_t ThetaF9test(Double_t* x, Double_t* p)
{
  Double_t  t=x[0];
  Double_t  f=0;
  if(t<p[2]){
    //  f=p[0]+p[1]*t+p[2]*sq(t);  //+p[3]*cub(t);
    f=p[0]+p[1]*sq(t);
  }
  else if(t>=p[2]){
    // f=p[0]+p[1]*p[4]+p[2]*sq(p[4])+p[3]*cub(p[4])+p[5]*(t-p[4])+p[6]*sq(t-p[4])+p[7]*cub(t-p[4])+p[8]*sq4(t-p[4])+p[3]*sq6(t-p[4]);
    f=p[0]+p[1]*sq(p[2])+p[3]*T3(t-p[2])+p[4]*T4(t-p[2])+p[5]*T5(t-p[2])+p[6]*T6(t-p[2])+p[7]*T7(t-p[2])+p[8]*T8(t-p[2]);
  }
  //  else if(t>=p[6]) f=p[0]+p[1]*sq(p[3]-p[2])+p[4]*(t-p[3])+p[5]*sq(t-p[3])+p//[7]*pow((t-p[6]),-p[8]);
  return f;
}

Double_t AlphaF(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]*T0()+p[1]*T1(a)+p[2]*T2(a)+p[3]*T3(a)+
    p[4]*T4(a)+p[5]*T5(a)+p[6]*T6(a)+p[7]*T7(a)+p[8]*T8(a);
  return dmult;  
}
Double_t RF6(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]*T0()+p[1]*T1(a)+p[2]*T2(a)+p[3]*T3(a)+
        p[4]*T4(a)+p[5]*T5(a);
  return dmult;  
}
Double_t RF4(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]+p[1]*(a-p[2])+p[3]*sq(a-p[2]);
  return dmult;  
}
Double_t RF5(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]+p[1]*(a-p[2])+p[3]*sq(a-p[2])+p[4]*cub(a-p[2]);
  return dmult;  
}
Double_t AlphaF3(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]+p[1]*T1(a)+p[2]*T2(a)+p[3]*T3(a);
  return dmult;  
}
Double_t AlphaF5(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  dmult=p[0]+p[1]*T1(a)+p[2]*T2(a)+p[3]*T3(a)+p[4]*T4(a);
  return dmult;  
}
Double_t AlphaFtest(Double_t* x,Double_t* p)
{
  Double_t  dmult=0;  
  Double_t  a=x[0];
  // p[2]=0;
  // p[6]=0;
  // p[7]=0;
  if(a<=p[3]) dmult=p[0]+p[1]*T1(a)+p[2]*T2(a);
  if (a>p[3]) dmult=p[0]+p[1]*T1(p[3])+p[2]*T2(p[3])+p[4]*(a-p[3])+p[5]*sq(a-p[3])+p[6]*cub(a-p[3])+p[7]*sq4(a-p[3]);
  /* if(a<=p[2]) dmult=p[0]+p[1]*a;
  if (a>p[2]) dmult=p[0]+p[1]*p[2]+p[3]*(a-p[2])+p[4]*sq(a-p[2])+p[5]*cub(a-p[2])+p[6]*sq4(a-p[2]);*/
  return dmult;  
}

Double_t test1(Double_t x,Double_t* par)
{
    return (par[0]*pow(x,4)+par[1]*pow(x,2)+par[2])*exp(sin(x));
}
Double_t test(Double_t x,Double_t* par)
{
    return par[0]*pow(x,3)-2*par[1]*x-5*par[2];
}
Int_t compar(const void* a,const void* b)
{
        Double_t x=*(Double_t*)a-*(Double_t *)b;
        return (x>0 ? 1:(x<0 ? -1:0));
};
Int_t comparD(Int_t UpDown,Double_t a,Double_t b)
{
    Double_t x=UpDown*(a-b);
    return (x>0 ? 1:(x<0 ? -1:0));
}
Int_t comparDRows(Int_t UpDown,Double_t* a,Double_t* b,Int_t n)
{
    Int_t x=comparD(UpDown,a[n],b[n]);
    return x;
}
void swapD(Double_t& a,Double_t& b)
{
        Double_t temp=a;
        a=b;
        b=temp;
};
void swapI(Int_t& a,Int_t& b)
{
        Int_t temp=a;
        a=b;
        b=temp;
};
void swapDRows(Double_t* a,Double_t* b,Int_t n)
{
    for(Int_t i=0;i<n;i++)
    {
	swapD(a[i],b[i]);
    }
}
/*
inline Double_t sq(Double_t x);{return (x*x);}
inline Double_t cub(Double_t x){return (x*x*x);}
inline Double_t sq4(Double_t x){return (x*x*x*x);}
inline Double_t sq5(Double_t x){return sq(x)*cub(x);}
inline Double_t sq6(Double_t x){return cub(x)*cub(x);}
inline Double_t sq7(Double_t x){return sq4(x)*cub(x);}
inline Double_t sq8(Double_t x){return sq4(x)*sq4(x);}
inline Double_t sq9(Double_t x){return sq4(x)*sq5(x);}
inline Double_t sq10(Double_t x){return sq5(x)*sq5(x);}
inline Double_t sq11(Double_t x){return sq5(x)*sq6(x);}

inline Double_t Min(Double_t a,Double_t b){return TMath::Min(a,b);}
inline Double_t Max(Double_t a,Double_t b){return TMath::Max(a,b);}

inline Double_t T0(){return 1;}
inline Double_t T1(Double_t x){return x;}
inline Double_t T2(Double_t x){return (2*sq(x)-1);}
inline Double_t T3(Double_t x){return (4*cub(x)-3*x);}
inline Double_t T4(Double_t x){return (8*sq4(x)-8*sq(x)+1);}
inline Double_t T5(Double_t x){return (16*sq5(x)-20*cub(x)+5*x);}
inline Double_t T6(Double_t x){return (32*sq6(x)-48*sq4(x)+18*sq(x)-1);}
inline Double_t T7(Double_t x){return (64*sq7(x)-112*sq5(x)+56*cub(x)-7*x);}
inline Double_t T8(Double_t x){return (128*sq8(x)-256*sq6(x)+160*sq4(x)-32*sq(x)+1);}
inline Double_t T9(Double_t x){return (256*sq9(x)-576*sq7(x)+432*sq5(x)-120*cub(x)+9*x);}
inline Double_t T10(Double_t x){return (512*sq10(x)-1280*sq8(x)+1120*sq6(x)-400*sq4(x)+50*sq(x)-1);}
inline Double_t T11(Double_t x){return (1024*sq11(x)-2816*sq9(x)+2816*sq7(x)-1232*sq5(x)+220*cub(x)-11*x);}
inline Double_t T12(Double_t x){return (2*x*T11(x)-T10(x));}
inline Double_t T13(Double_t x){return (2*x*T12(x)-T11(x));}
inline Double_t T14(Double_t x){return (2*x*T13(x)-T12(x));}
inline Double_t T15(Double_t x){return (2*x*T14(x)-T13(x));}
inline Double_t T16(Double_t x){return (2*x*T15(x)-T14(x));}
inline Double_t T17(Double_t x){return (2*x*T16(x)-T15(x));}
inline Double_t T18(Double_t x){return (2*x*T17(x)-T16(x));}
inline Double_t T19(Double_t x){return (2*x*T18(x)-T17(x));}
inline Double_t T20(Double_t x){return (2*x*T19(x)-T18(x));}
  */
Double_t factorial(Double_t* N)
{
    Double_t f=1.0;
    if(N) for(int i=1;i<=*N;i++) f*=(Double_t)i;
    return f;
};
Double_t PartPuasson(Double_t Nexpect,Double_t Nobserv)
{
    Double_t pp;
    pp=Nexpect-Nobserv;
    if(Nobserv>0&&Nexpect>0) pp+=Nobserv*log(Nobserv/Nexpect);
    return 2.0*pp;
};
Double_t PartMultinom(Double_t Nexpect,Double_t Nobserv)
{
    Double_t pp=0;
    if(Nobserv>0&&Nexpect>0) pp+=Nobserv*log(Nobserv/Nexpect);
    return 2.0*pp;

};
Double_t chebysh_exp(Double_t* x,Double_t* par){
	Double_t y;
	y=x[0];
	Double_t z;
	Double_t cut=par[12];
	if(x[0]<=cut) z=0.5*par[0]+par[1]*T1(y)+par[2]*T2(y)+par[3]*T3(y)+par[4]*T4(y)+par[5]*T5(y)+par[6]*T6(y)+par[7]*T7(y)+par[8]*T8(y)+par[9]*T9(y)+par[10]*T10(y)+par[11]*T11(y);
       else z=(0.5*par[0]+par[1]*T1(cut)+par[2]*T2(cut)+par[3]*T3(cut)+par[4]*T4(cut)+par[5]*T5(cut)+par[6]*T6(cut)+par[7]*T7(cut)+par[8]*T8(cut)+par[9]*T9(cut)+par[10]*T10(cut)+par[11]*T11(cut))*exp(-(x[0]-cut)/par[13]);
//	+par[4]*T4(y));
//   else z=(0.5*par[0]+par[1]*T1(par[6])+par[2]*T2(par[6])+par[3]*T3(par[6])+par[4]*T4(par[6]))*exp(-(1.0/y-par[6])/par[5]);
	return z;
};
Double_t power(Double_t* x, Double_t* par){
	Double_t y=x[0];
	return   par[0]*pow(y,par[1]);
};

Double_t powerN(Double_t* x, Double_t* p)
{
	Double_t y=x[0];
	return   p[0]*pow(y/42,p[1]);
};

Double_t powerL(Double_t* x, Double_t* p)
{
	Double_t y=x[0];	
	return   p[0]*pow(y/21.,p[1]);
};
Double_t powerBB(Double_t* x, Double_t* p)
{
	Double_t y=x[0];
	return   p[0]*(1+p[1]*sq(y-p[2]));
};
Double_t powerP(Double_t* x, Double_t* p)
{
	Double_t y=x[0];
	return   p[0]*(1+p[1]*sq(y-p[2]));
};


Double_t linear(Double_t* x, Double_t* par)
{
	Double_t y=x[0];
	return   par[0]+par[1]*y;
};
Double_t linear_(Double_t x, Double_t* par)
{
	Double_t y=x;
	return   par[0]+par[1]*exp(-sq(y)/2)/sqrt(2*TMath::Pi());
};

Double_t linear1(Double_t* x, Double_t* par)
{
	Double_t y=x[0];
	return   par[0]+par[1]*(y-250.);
};
Double_t linear2(Double_t* x, Double_t* par)
{
	Double_t y=x[0];
	return   par[0]+par[1]*(y-11);
};
Double_t linear0(Double_t* x, Double_t* par)
{
	return   par[0];
};
Double_t along(Double_t* x, Double_t* par){
	Double_t y=x[0];
	return   par[0]*(1.+par[1]*cos(2*TMath::Pi()*(y+par[3])/par[2]));
};
Double_t alongz(Double_t* x, Double_t* par){
	Double_t y=x[0];
	return   (par[0]+par[1]*y+par[2]*y*y);
};
Double_t ZF(Double_t* x, Double_t* par){
	Double_t y=x[0];
	return   par[0]*(1+par[1]*y+par[2]*sq(y)+par[3]*cub(y)+par[4]*sq4(y));
};
Double_t alongzA(Double_t* x, Double_t* par){
    Double_t y=x[0];
    Double_t z;
    if(fabs(y)<=par[2]) z=par[0];
    else z=par[0]+par[1]*sq4(y-par[2]);
    return  z;
};
Double_t standard_landau(Double_t* x, Double_t* par){
	Double_t y=x[0];
//  Emp=par[0],ksi=par[1],;
	Double_t v;
        v=(y-par[0])/par[1];
	return  par[2]*exp(-0.5*(v+exp(-v)));
};
Double_t standard_Landau(Double_t* x, Double_t* par){
	Double_t y=x[0];
//  Emp=par[0],ksi=par[1],;
	Double_t v;
        v=(y-par[1])/par[2];
	return  par[0]*exp(-0.5*(v+exp(-v)));
};

Double_t landauFF(Double_t* x, Double_t* par){
	Double_t y=x[0];
//  Emp=par[0],ksi=par[1],;
	Double_t v;
        v=(y-par[0])/par[1]-0.225;
	return  par[2]*exp(-0.5*(v+exp(-v)));
};
Double_t de_loss(Double_t* x, Double_t* par)
{
	Double_t  y=x[0];
	Double_t  b=y/sqrt(sq(par[5])+sq(y));
	return par[0]*(par[1]-pow(b,2)+log(par[2]+1/pow(b*(1-b*b),par[3])))/pow(b,par[4]) ;
};
Double_t de_loss_fit(Double_t x, Double_t* par)
{
	Double_t  y=x;
	return par[0]*(par[1]-pow(y,2)+log(par[2]+1/pow(y*(1-y*y),par[3])))/pow(y,par[4]) ;
};
Double_t ChargeR(Double_t ChargeI, Double_t* par)
{
	Double_t  y=ChargeI;
        Double_t  Const=3.125*(10e+3) ;// 5*(10e-16)/(1.6*(10e-19))
	return  4.*(y-par[1])*Const*exp(par[2]/par[3])*par[0]*par[5]/par[4];
};
Double_t dE_Real(Double_t ChargeI, Double_t* par)
{
	Double_t  y=ChargeI;
        Double_t  Const=3.125*(10e+3) ;// 5*(10e-16)/(1.6*(10e-19))
	return  4.*(y-par[1])*Const*exp(par[2]/par[3])*par[0]*15.*par[5]/par[4];
};
Double_t Charge_R(Double_t ChargeI, Double_t* par)
{
	Double_t  y=ChargeI;
        Double_t  Const=3.125*(10e+3) ;// 5*(10e-16)/(1.6*(10e-19))
	return  par[0]*((y-par[1])*Const/par[4])*60.*4;
};

Double_t signum(Double_t x)
{
	Double_t s=1;
	if(x<0) s=-1;
	return s;
};
void  ScalePoints(Double_t* a,Double_t xmin,Double_t xmax,Double_t xrmin,Double_t xrmax,Double_t xshift,Int_t nleft,Int_t nr,Int_t nright){
    Double_t dlleft,dlright,dlr;
    dlleft=(xrmin-xmin)/nleft;
    dlr=(xrmax-xrmin)/(nr-1);
    dlright=(xmax-xrmax)/nright;
    for(Int_t i=0;i<nleft;i++){
        a[i]=dlleft*i+xmin+xshift;
    }
    for(Int_t i=0;i<nr;i++){
        a[i+nleft]=xrmin+dlr*i+xshift;
    }
    for(Int_t i=1;i<=nright;i++){
        a[i+nleft+nr-1]=xrmax+dlright*i+xshift;
    }
};
void devide_l(Int_t n,Double_t* a,Double_t xmin,Double_t xmax,Double_t xbet,Int_t nmin){
	Double_t dlmin,dlmax;
	dlmin=(xbet-xmin)/nmin;
        dlmax=(xmax-xbet)/(n-nmin-1);
	for(Int_t i=0;i<=nmin;i++){
		a[i]=dlmin*i;
	}
        for(Int_t i=(nmin+1);i<n;i++){
		a[i]=dlmin*nmin+dlmax*(i-nmin)+xmin;
 	}
};
Double_t out_l(Int_t i,Double_t* a){
	Double_t sa;
	sa=(a[i]+a[i+1])/2.0;
	return sa;
};

void aver_value(Int_t dim,Int_t *dnp,Double_t* v,Double_t* s,Double_t* vn,Double_t* sn){
    Int_t ij=0;
    Double_t w;
    for(Int_t i=0;i<dim;i++)
    {
        w=0.;
        for(Int_t j=ij;j<(ij+dnp[i]);j++){
            w+=(1/(s[j]*s[j]));
        }
            sn[i]=1./sqrt(w);
        for(Int_t j=ij;j<(ij+dnp[i]);j++){
            vn[i]+=(v[j]/(s[j]*s[j]));
        }
            vn[i]/=w;
            ij+=dnp[i];
    }
};
void aver_value_scan(Int_t npe,Int_t isc,Double_t** v,Double_t** s,Double_t* vn,Double_t* sn){
//    Int_t ij=0;
    Double_t w;
    for(Int_t i=0;i<npe;i++){
        w=0.;
        for(Int_t j=0;j<isc;j++){
            w+=(1./(s[i][j]*s[i][j]));
        }
        sn[i]=1./sqrt(w);
        for(Int_t j=0;j<isc;j++){
            vn[i]+=(v[i][j]/(s[i][j]*s[i][j]));
        }
        vn[i]/=w;
    }
};
void aver_value_cut(Int_t nps,Int_t ns,Double_t* v,Double_t* s,Double_t* vn,Double_t* sn){
    Int_t ij=0;
    Double_t w;
    for(Int_t i=0;i<nps;i++){
        w=0.;
        for(Int_t j=0;j<ns;j++){
            ij=i+j*nps;
            w+=(1./(s[ij]*s[ij]));
        }
           sn[i]=1./sqrt(w);
        for(Int_t j=0;j<ns;j++){
            ij=i+j*nps;
            vn[i]+=(v[ij]/(s[ij]*s[ij]));
        }
            vn[i]/=w;
    }
};
////// for dE/dx
Double_t dist_exp(Double_t* x, Double_t* par)
{
	Double_t y=x[0];
	Double_t z=par[0];
	if(y>par[2]) z=par[0]*exp(-(y-par[2])/par[1]);
        return z;
};
Float_t define_bins(Int_t nslices,const Int_t* scaleweight,
	Float_t xmin,Float_t xmax,Float_t* bincenter,Int_t* defplace)
{
	Float_t step=0;
	Int_t dev=0;
	Int_t ind=0;
	Float_t* leftedge=new Float_t [nslices];
	for(Int_t i=0;i<nslices;i++)
	{
		dev+=scaleweight[i];
	}
	
	if(dev>0&&xmax>xmin)
	{
		step=(xmax-xmin)/dev;
        for(Int_t i=0;i<nslices;i++)
		{
			if(i!=0)
			{
         		leftedge[i]=leftedge[i-1]+step*scaleweight[i-1];
            	bincenter[i-1]=(leftedge[i]+leftedge[i-1])/2.;
			}
			else
			{
				leftedge[0]=xmin;
			}
			if(i==nslices-1)
			{
				bincenter[i]=xmax-step*scaleweight[i]/2.;
			}
			for(Int_t j=0;j<scaleweight[i];j++)
			{
				defplace[ind+j]=i;
			}
			ind+=scaleweight[i];
		}
		
	}
	else
	{
		cout<<"\n mistake in bins definition ! step=0"<<endl;
	}
		
    return step;
};
Float_t define_bins_(Int_t nslices,Float_t xmin,Float_t xmax,Float_t* bincenter)
{
	Float_t step=0;
	Int_t ind=0;
	Float_t* leftedge=new Float_t [nslices];
	if(nslices>0&&xmax>xmin)
	{
		step=(xmax-xmin)/nslices;
        for(Int_t i=0;i<nslices;i++)
		{
			if(i!=0)
			{
         		leftedge[i]=leftedge[i-1]+step;
            	bincenter[i-1]=(leftedge[i]+leftedge[i-1])/2.;
			}
			else
			{
				leftedge[0]=xmin;
			}
			if(i==nslices-1)
			{
				bincenter[i]=xmax-step/2.;
			}
		}
		
	}
	else
	{
		cout<<"\n mistake in bins definition ! step=0"<<endl;
	}
		
    return step;
};
Double_t BB(Double_t* x,Double_t* par){
	Double_t  y=x[0];
	Double_t  bb;
	Double_t  gamma;
	Double_t  beta;
	gamma=sqrt(1+sq(y));
	beta=y/gamma;
//	PARAMETER (p1=31.911,p2=12.257,p3=6.9468e-05,p4=2.0851,p5=2.3047)
	bb=par[0]*(par[1]-pow(beta,par[4])-log(par[2]+pow((1./y),par[3])))/pow(beta,par[4]);
    return bb;
};
Double_t BBfix(Double_t x)
{
	Double_t y=x;
	Double_t  par[5];
	Double_t  bb;
	Double_t  gamma;
	Double_t  beta;
	gamma=sqrt(1+sq(y));
	beta=y/gamma;
    /*
	par[0]=31.911;
	par[1]=12.257;
	par[2]=6.9468e-05;
	par[3]=2.0851;
	par[4]=2.3047;
    */
	par[0]=15.48;
	par[1]=23.28;
	par[2]=0.1232e-8;
	par[3]=4.241;
	par[4]=2.427;
	bb=par[0]*(par[1]-pow(beta,par[4])-log(par[2]+pow((1./y),par[3])))/pow(beta,par[4]);
    return bb;
};
Double_t BB3(Double_t* x,Double_t* par)
{
	Double_t  y=x[0];
	Double_t  bb;
//	PARAMETER (p1=31.911,p2=12.257,p3=6.9468e-05,p4=2.0851,p5=2.3047)
	bb=par[0]/pow(y,par[1])+par[2];
    return bb;
};
// BaBar functions
Double_t ChargeCorrCoupal(Double_t chargein,Double_t sindip, Double_t* CoupalPar)
{
	Double_t y=sindip;
	Double_t corrT[2];
	Double_t chargeout;
	corrT[0] = CoupalPar[0]/sqrt(sq(y) + sq(CoupalPar[1]));
	corrT[1] = 0.5*sq(CoupalPar[0])/(sq(y) + sq(CoupalPar[1]));
	chargeout =chargein*(1.+corrT[0]*chargein+corrT[1]*sq(chargein));//-1/3*corrT[2]*cub(chargein));
    return chargeout;
};
Double_t ChargeCorrCheby(Double_t chargein,Double_t sindip, Double_t* ChebyPar)
{
	Double_t y;
	Double_t chargeout;
	Double_t CorrCheby;
	y=2*fabs(sindip)-1.;
	CorrCheby = ChebyPar[0]+ChebyPar[1]*cos(acos(y))
	+ChebyPar[2]*cos(2.*acos(y))
	+ChebyPar[3]*cos(3.*acos(y))
	+ChebyPar[4]*cos(4.*acos(y));
    chargeout=chargein/CorrCheby;
	return chargeout;
};
Double_t ChebyP4(Double_t* sindip, Double_t* par)
{
	Double_t y=2*fabs(sindip[0])-1;
	Double_t normch;
	normch=par[0]+par[1]*cos(acos(y))
	+par[2]*cos(2.*acos(y))
	+par[3]*cos(3.*acos(y))
	+par[4]*cos(4.*acos(y));
	return normch;
};
Double_t ChebyP3(Double_t* z, Double_t* par)
{
	Double_t y=z[0];
	Double_t normz;
	normz=par[0]+par[1]*T1(y)+par[2]*T2(y)+par[3]*T3(y);
	return normz;
};
Float_t ChebyP3Fix(Float_t z, Float_t* par)
{
	Float_t y=z;
	Float_t normz;
	normz=par[0]+par[1]*T1(y)+par[2]*T2(y)+par[3]*T3(y);
	return normz;
};
Float_t ChebyP4Fix(Float_t sindip, Double_t* par)
{
	Float_t y=2*fabs(sindip)-1;
	Float_t normch;
	normch=par[0]+par[1]*cos(acos(y))
	+par[2]*cos(2.*acos(y))
	+par[3]*cos(3.*acos(y))
	+par[4]*cos(4.*acos(y));
	return normch;
};
Double_t ChargeAfterSat(Double_t sindip,Double_t rhoraw, Double_t* par)
{
	Double_t mult=par[0];
	Double_t k=par[1];
	Double_t d=par[2];
	Double_t rhomeas;
//	qmeas=mult*exp(-(k*rhoraw))*(T0()+p1*T1(smooth)+p2*T2(smooth)+p3*T3(smooth)+p4*T4(smooth));
//  qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)/pow(sq(a*rhoraw)+1,c)));
//  qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)*pow(sq(a*rhoraw)+1,c)));
//  rhomeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)*pow(sq(a*rhoraw)+1,c)));
    rhomeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)));
	return rhomeas;
};
Double_t ChargeAfterSat_(Double_t* x, Double_t* par)
{
	Double_t sindip=x[0];
	Double_t rhoraw=x[1];
	Double_t mult=par[0];
	Double_t k=par[1];
	Double_t d=par[2];
	Double_t rhomeas;
//	qmeas=mult*exp(-(k*rhoraw))*(T0()+p1*T1(smooth)+p2*T2(smooth)+p3*T3(smooth)+p4*T4(smooth));
//  qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)/pow(sq(a*rhoraw)+1,c)));
//  qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)*pow(sq(a*rhoraw)+1,c)));
//  rhomeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)*pow(sq(a*rhoraw)+1,c)));
    rhomeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)));
	return rhomeas;
};

Double_t ChargeAfterSatFix(Double_t* x, Double_t* par)
{
	Double_t sindip=x[0];
	Double_t rhoin=x[1];
    /*
	Double_t p0=par[3];
	Double_t p1=par[4];
	Double_t p2=par[5];
	Double_t p3=par[6];
    */
	Double_t mult=par[0];
	Double_t k=par[1];
	Double_t d=par[2];
//	Double_t a=par[3];
//	Double_t c=par[4];
//    Double_t s=par[5];
    
   
    /*
	Double_t p2=par[3];
	Double_t p3=par[4];
	Double_t p4=par[5];
	Double_t p5=par[6];
    */
//	Double_t smooth=sqrt(sq(sindip)+sq(d));;
	Double_t rhoout;
//	qmeas=mult*exp(-(k*rhoraw))*(T0()+p1*T1(smooth)+p2*T2(smooth)+p3*T3(smooth)+p4*T4(smooth));
 // 
//   qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+(sq(d)/(sq(a*rhoraw)+sq(c)))));
 //  qmeas=mult*exp(-(k*rhoraw)/sqrt(sq(sindip)+sq(d)/pow(sq(a*rhoraw)+1,c)));
//	rhoout=rhoin/mult*(1.+k/mult*rhoin/sqrt(sq(sindip)+sq(d)*pow(sq(a*rhoin/mult)+1,c))+0.5*sq(k/mult*rhoin)/(sq(sindip)+sq(d)*pow(sq(a*rhoin/mult)+1,c)));
	rhoout=rhoin/mult*(1.+k/mult*rhoin/sqrt(sq(sindip)+sq(d))+0.5*sq(k/mult*rhoin)/(sq(sindip)+sq(d)));
//	rhoout-=rhoin/mult*cub(k/mult*rhoin)/pow(sq(sindip)+sq(d),1.5)/3.;
	return rhoout;
};
Double_t CalibrationCharge1(Double_t* x, Double_t* p)
{
    Double_t F=0;
    Double_t y=x[0];
    F=p[0]*(log(1+y)/(1+sq(y))+p[1]*pow(y,p[2]));
    return F;
}
Double_t CalibrationCharge2(Double_t* x, Double_t* p)
{
    Double_t F=0;
    Double_t y=x[0];
    F=p[0]+p[1]*(y-p[2])+p[3]*(cub(y-p[4])-cub(p[2]-p[4]));
    return F;
}
Double_t CalibrationCharge3(Double_t* x, Double_t* p)
{
    Double_t F=0;
    Double_t y=x[0];
    F=p[0]+p[1]*(y-p[2])+p[3]*(1./pow(p[4]-y,p[5])-1./pow(p[4]-p[2],p[5]));
    return F;
}
Double_t CalibrationChargeGlobal(Double_t* x, Double_t* p)
{
    Double_t F=0;
    Double_t y=x[0];
    Double_t r1=p[3];
    Double_t r2=p[7];
    if(y<r1)
    {
        F+=p[0]*(log(1+y)/(1+sq(y))+p[1]*pow(y,p[2]));
    }
    else
    {
        F+=p[0]*(log(1+r1)/(1+sq(r1))+p[1]*pow(r1,p[2]));
        if(y<r2)
        {
            F+=p[4]*(y-r1)+p[5]*(cub(y-p[6])-cub(r1-p[6]));
        }
        else
        {
            F+=p[4]*(r2-r1)+p[5]*(cub(r2-p[6])-cub(r1-p[6]));
            F+=p[8]*(y-r2)+p[9]*(1./pow(p[10]-y,p[11])-1./pow(p[10]-r2,p[11]));
        }
    }
    return F;
}
Double_t squareF(Double_t* x, Double_t* p)
{
    Double_t y=x[0];
    Double_t F=0;
    F=p[0]+p[1]*sq(y-p[2]);
    return F;
}
Double_t XtC(Double_t* x, Double_t* par)
{
    Double_t F=0;
    Double_t y=x[0];
    Double_t r1=par[2];
    Double_t r2=par[3];
    Double_t r3=par[6];    
    Double_t r4=par[10];    
    if(y<r1)      F+=par[0]*log(1+y);       
    else
     {
          F+=par[0]*log(1+r1);
          if(y<r2)  F+=par[1]*pow(y-r1,par[4]);        
         else
         {
            F+=par[1]*pow(r2-r1,par[4]);
            if(y<r3) F+=par[5]*pow(y-r2,par[7]);          
           else
           {
	F+=par[5]*pow(r3-r2,par[7]);
	 if(y<r4)       F+=par[8]*pow(y-r3,par[9]);	 
	 else           F+=par[8]*pow(r4-r3,par[9]); 
           }
       }
    }
    return F;
}
