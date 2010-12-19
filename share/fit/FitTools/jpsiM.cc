#include <iostream.h>
#include "FitTools/jpsiM.h"
#include "FitTools/MathLibrary.h"
#include "FitTools/PidDef.h"
#include <TRandom.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
Double_t M_h_f_(Double_t* a,Double_t *b,Double_t *x)
{
    Double_t msn=0;
    Double_t mult=1;
    Double_t nd;
    for(Int_t n=0;n<=nsequent;n++)
    {
        nd=(Double_t)(n);
        if(n>=1) mult*=((*a+nd-1.0)*(*x)/((*b+nd-1.0)*nd));
        msn+=mult;
    }
    return msn;
};
Double_t M_h_f_d(Double_t* a,Double_t *b,Double_t *x)
{
    Double_t msn=0;
    Double_t mult;
    Double_t nd;
    mult=(*a)/(*b);
    for(Int_t n=0;n<=nsequent;n++)
    {
        nd=(Double_t)(n);
        if(n>=1) mult*=((*a+nd)*(*x)/((*b+nd)*nd));
        msn+=mult;
    }
    return msn;
};
Double_t fsn_(Double_t* t,Double_t* z,Int_t n,Int_t s)
{
    Double_t fsn=1,mf=1;
    if(s>0)
	for(Int_t nc=1;nc<=(n-1);nc+=2)
	{
            mf*=((*t)-nc)*((*t)-(nc+1))/((*z)*(*z)*(nc+1));
            fsn+=mf;
        }
    if(s<0)
	for(Int_t nc=0;nc<=(n-2);nc+=2)
	{
            mf*=(-1)*((*t)+nc)*((*t)+(nc+1))/((*z)*(*z)*(nc+2));
            fsn+=mf;
        }
    return  fsn;
};
Double_t fsn_d(Double_t* t,Double_t* z,Int_t n,Int_t s)
{
    Double_t fsn=0,mf=1;
    if(s>0)
	for(Int_t nc=1;nc<=(n-1);nc+=2)
	{
            mf*=((*t)-nc)*((*t)-(nc+1))/((*z)*(*z)*(nc+1));
            fsn-=mf*(nc+1)/(*z);
        }
    if(s<0)
	for(Int_t nc=0;nc<=(n-2);nc+=2)
	{
            mf*=(-1)*((*t)+nc)*((*t)+(nc+1))/((*z)*(*z)*(nc+2));
            fsn-=mf*(nc+2)/(*z);
        }
    return  fsn;
};
Double_t Fzt_Aprox(Double_t* z,Double_t *t)
{
    Double_t fzt=0,fztr=0;
    fzt=TMath::Gamma(1.+(*t))*exp(-0.5*sq(*z))/sqrt(2*TMath::Pi());
    fzt*=pow((0.31+sq(*z)-0.73*(*z)/sqrt(1+sq((*z)/(1+1.37*(*t))))),-0.5*(*t));
    //    cout<<"GR_:"<<fzt<<"z:"<<zz<<"beta:"<<*t<<endl;
    if((*z)>0){
        fztr=(*t)*pow((*z),(*t))*pow((*z),2.18)/(1.+pow((*z),3.18));
        fztr*=(1+0.5*((*t)-1.)*((*t)-2.)/(sq((*z)-46./(sq(*z)+21.))+2.44+1.5*(*t)));
        fzt+=fztr;
    }
    return fzt;
};
Double_t Fzt_(Double_t* z,Double_t *t)
{
    Double_t fzt;
    Double_t tt;
    Double_t a;
    Double_t b;
    Double_t x;
    tt=(*t)*0.5;
    a=0.5-tt;
    b=0.5;
    x=-0.5*(*z)*(*z);
    fzt=pow(2,tt)*TMath::Gamma(1+tt)*M_h_f_(&a,&b,&x)*sqr_2pi;
    a=1.0-tt;
    b=1.5;
    fzt+=(tt*TMath::Gamma(1+(*t))*(*z)*M_h_f_(&a,&b,&x)/(pow(2,tt)*TMath::Gamma(1+tt)));
    	if((*z)>7.5)
    	{
          fzt=fsn_(t,z,8,1)*(*t)*pow((*z),(*t)-1);
    //    fzt=Fzt_Aprox(z,t);
    	}
    //  fzt=(*t)*pow((*z),(*t)-1)*(1+(1-(*t))*(2-(*t))*0.5/sq(*z));
    	if((*z)<-4.5)
    	{
//  	        fzt=Fzt_Aprox(z,t);
        	fzt=TMath::Gamma(1+(*t))*sqr_2pi*exp(-0.5*sq(*z));
        	fzt*=pow(fabs(*z),-(*t))*fsn_(t,z,8,-1);
    	}
    return fzt;
};
Double_t Fzt_1(Double_t z,Double_t *t)
{
//    Double_t F=Fzt_(&z,t);
    return Fzt_(&z,t);
}
Double_t Fzt_d(Double_t* z,Double_t *t){
    Double_t fzt;
    Double_t fzt_c;
    Double_t tt;
    Double_t a;
    Double_t b;
    Double_t x;
    tt=(*t)*0.5;
    a=0.5-tt;
    b=0.5;
    x=-0.5*(*z)*(*z);
    fzt=-pow(2,tt)*TMath::Gamma(1+tt)*M_h_f_d(&a,&b,&x)*sqr_2pi*(*z);
    a=1.0-tt;
    b=1.5;
    fzt+=(tt*TMath::Gamma(1+(*t))*M_h_f_(&a,&b,&x)/(pow(2,tt)*TMath::Gamma(1+tt)));
    fzt-=(tt*TMath::Gamma(1+(*t))*(*z)*(*z)*M_h_f_d(&a,&b,&x)/(pow(2,tt)*TMath::Gamma(1+tt)));
    if((*z)>7.5)
    {
        fzt=fsn_d(t,z,8,1)*(*t)*pow((*z),(*t)-1);
        fzt+=fsn_(t,z,8,1)*(*t)*((*t)-1.)*pow((*z),(*t)-2.);
    }
    if((*z)<-4.5)
    {
        fzt=TMath::Gamma(1+(*t))*sqr_2pi*exp(-0.5*sq(*z));
        fzt*=pow(fabs(*z),-(*t))*fsn_d(t,z,8,-1);
        fzt_c=TMath::Gamma(1+(*t))*sqr_2pi*exp(-0.5*sq(*z));
        fzt_c*=pow(fabs(*z),(-(*t)-1.))*(*t)*fsn_(t,z,8,-1);
        fzt+=fzt_c;
        fzt_c-=((*z)*TMath::Gamma(1+(*t))*sqr_2pi*exp(-0.5*sq(*z)));
        fzt_c*=pow(fabs(*z),-(*t))*fsn_(t,z,8,-1);
        fzt+=fzt_c;
    }
    return fzt;
};


Double_t xsecnbsim(Double_t* Eb,Double_t* par)
{
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM=0,SiW=0,Mpsi=0,Gee,pi,beta,DeltaE,Delta2,SS,z;
    Double_t xs;

    if(Pid==0)
    {
        Mpsi=MJPsi;
        hM=par[2]+0.5*MJPsi;
	SiW=par[3]*sq(y/hM);
	Gee=GeeJPsi;

    }
    if(Pid!=0)
    {
        Mpsi=MPsiPrime;
    	hM=par[2]+0.5*MPsiPrime;
	SiW=par[3]*sq(y/hM);
	Gee=GeePsiPrime;
    }

//  Mpsi=3686.00;
    //  Mpsi=par[2]*2.;
    pi=3.14159265;
    beta=4*alfa/pi*(log(Mpsi/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    Double_t parPsiPrime[2]={0,1.33};
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;

//    +par[1]*GrIntPsiPrime(y,0.000001,12,parPsiPrime);
//    DeltaE=alfa/pi*(pow(pi,2)-0.5)+0.75*beta;
//  SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
//    z=2*(y-hM)/SiW;
 // xs=par[0]+par[1]*Fzt_(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW*SS;
//    xs=par[0]*sq(hM/y)+par[1]*Fzt_(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;

//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;


};


Double_t xs_bhabha(Double_t Eb,Double_t iscan,Double_t* par)
{

    Double_t y=Eb;
    Double_t hM=0;
    if(Pid==0) hM=1548.44;
    if(Pid==1) hM=1843.;

    Double_t xsbb=0;
    if(iscan==0)
    {
	xsbb=par[0]*sq(hM/y);
    }
    if(iscan==1)
    {
	xsbb=par[1]*sq(hM/y);
    }
    if(iscan==2)
    {
	xsbb=par[2]*sq(hM/y);
    }
    if(iscan==3)
    {
	xsbb=par[3]*sq(hM/y);
    }
    return xsbb;
}
Double_t xsbhabha_(Double_t Eb,Double_t* par)
{
    Double_t y=Eb;

    Double_t hM=0;
    if(Pid==0) hM=1548.44;
   // if(Pid==0) hM=1546.00;
    if(Pid==1) hM=1843.00;
    Double_t xsbb;
    xsbb=par[0]*sq(hM/y);
    return xsbb;
}
Double_t xsbhabha(Double_t Eb,Double_t iscan,Double_t* par)
{

    Double_t y=Eb;
    Double_t xsbb=0;
    Double_t hM=0;
    if(Pid==0) hM=1548.44;
    if(Pid==1) hM=1843.;
    if(iscan==0)
    {
	xsbb=par[13]*sq(hM/y);
    }
    if(iscan==1)
    {
	xsbb=par[13]*sq(hM/y);
    }
    if(iscan==2)
    {
	xsbb=par[13]*sq(hM/y);
    }
    if(iscan==3)
    {
	xsbb=par[13]*sq(hM/y);
    }
    return xsbb;
}


Double_t xsecnbScansSeparate(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi,Gee=0,beta,DeltaE,Delta2,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.12e-3;
	Mpsi=3686.;
    }

    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
//    beta1=beta+1;
   // DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;
};
Double_t xsecnbYpsilon(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=0.5*Eb;
    Double_t hM=0,SiW=0,Gee=0,beta,DeltaE,Delta2,z,SS;
    Double_t xs;
    hM=par[2]+MUpsilon1S*0.5;
    SiW=par[3]*sq(y/hM);
    Gee=GeeUpsilon1S;

    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
   // xs=fabs(par[0])+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;
};
Double_t GrIntYpsilon(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb*0.5;
    Double_t hM=0,SiW=0,Gee=0,beta,DeltaE,Delta2,SS;
    Double_t Gtot=0,Gh=0;
    Double_t ParF[5];
    Double_t xs,Finter;
    hM=par[2]+0.5*MUpsilon1S;
    SiW=par[3]*sq(y/hM);
    Gee=GeeUpsilon1S;
    Gtot=GtotUpsilon1S;
    Double_t M=2*hM;
    beta=4*alfa/TMath::Pi()*(log(M/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
  //  SS=12*TMath::Pi()/sq(M)*Gee*Gh/Gtot/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
  //  SS=12*TMath::Pi()/sq(M)*Gee*0.0011/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    SS=12*TMath::Pi()/sq(M)*0.001216/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
//    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    ParF[0]=beta;
    ParF[1]=M;
    ParF[2]=y*2;
    ParF[3]=SiW;
    ParF[4]=Gtot;
    Double_t rIntegral=0;
//#if 0
//    rIntegral=HANDLE_DGAUSS(FuncR,2*y-9*SiW,2*y+9*SiW,ParF,epsilon);
 //   rIntegral=HANDLE_DGAUSS(FuncRYpsilon,2*y-RangeE,2*y+RangeE,ParF,epsilon_);
//#else
   rIntegral=HANDLE_DGAUSS(FuncRUpsilon,2*y-10,2*y+10,ParF,1.e-7)+
      HANDLE_DGAUSS(FuncRUpsilon,2*y-20.,2*y-10,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRUpsilon,2*y-50.,2*y-20.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncRUpsilon,2*y+10,2*y+20.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRUpsilon,2*y+20.,2*y+50.,ParF,1.e-3);
//#endif
    xs=rIntegral/sqrt(2*TMath::Pi())/SiW*SS*par[1]+fabs(par[0])*sq(hM/y);
    return xs;
};
Double_t FuncRUpsilon(Double_t W, Double_t* parf)
{
      Double_t Func;
      Double_t beta=parf[0];
      Double_t M=parf[1];
      Double_t Wb=parf[2];
      Double_t Sw=parf[3];
      Double_t Gtot=parf[4];
      Double_t Fg=exp(-0.5*sq((W-Wb)/Sw));
      Double_t phi=TMath::ATan2(Gtot*0.5,M-W)*(1.-beta);
      if(phi<0) phi=phi+TMath::Pi();
      Func=pow((0.5*M/sqrt(sq(W-M)+sq(Gtot*0.5))),(1.-beta))*sin(phi)*Fg;
      return Func;
};

Double_t GrInt(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi,Gee=0,beta,DeltaE,Delta2,SS,MD;
    Double_t Gtot=0,Bll=0,R=2.5,Gh=0;
    Double_t ParF[6];
    Double_t xs,Finter;
    MD=0.5*(MDc+MD0);
    if(Pid==0)
    {
        Mpsi=MJPsi;
        hM=par[2]+0.5*Mpsi;
   //! 	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.14e-3;
//!	Mpsi=3096.88;
//	Gtot=0.087;
//	Bll=Gee/Gtot;
	Bll=0.0593;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }

    if(Pid==1)
    {
        Mpsi=MPsiPrime;
    	hM=par[2]+0.5*Mpsi;
      //  hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.12e-3;
	Mpsi=3686.;
//	Gtot=0.277;
//	Bll=Gee/Gtot;
	Bll=0.0090;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }
    if(Pid==2)
    {
        Mpsi=MPsiDoublePrime;
    	hM=par[2]+0.5*Mpsi;
	SiW=par[3]*sq(y/hM);
	Gee=0.26e-3;
//	Bll=Gee/Gtot;
	Bll=0.0000112;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }

    Double_t M=2*hM;
    beta=4*alfa/TMath::Pi()*(log(M/me)-0.5);
//  Finter=-Bll*R*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    Finter=-Bll*R/(1.-2.*Bll)*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=12*TMath::Pi()/sq(M)*Gee*Gh/Gtot/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    ParF[0]=beta;
    ParF[1]=M;
    ParF[2]=y*2;
    ParF[3]=SiW;
    ParF[4]=Gtot;
    ParF[5]=Finter;
    Double_t rIntegral=0;
#if 0
//    rIntegral=HANDLE_DGAUSS(FuncR,2*y-9*SiW,2*y+9*SiW,ParF,epsilon);
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-DGAUSS_RANGE,2*y+DGAUSS_RANGE,ParF,epsilon);
#else
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-3.5,2*y+3.5,ParF,1.e-6)+
      HANDLE_DGAUSS(FuncR,2*y-7.,2*y-3.5,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y-15.,2*y-7.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncR,2*y+3.5,2*y+7.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y+7.,2*y+15.,ParF,1.e-3);
#endif
    xs=rIntegral/sqrt(2*TMath::Pi())/SiW*SS*par[1]+fabs(par[0])*sq(hM/y);
    return xs;
};
Double_t GrIntNew(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi,Gee=0,beta,DeltaE,Delta2,SS,MD;
    Double_t Gtot=0,Bll=0,R=2.5,Gh=0;
    Double_t ParF[6];
    Double_t xs,Finter;

    MD=0.5*(MDc+MD0);
   
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
//	Gtot=0.087;
//	Bll=Gee/Gtot;
	Bll=0.059;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }

    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.12e-3;
	Mpsi=3686.;
//	Gtot=0.277;
//	Bll=Gee/Gtot;
	Bll=0.0090;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }
     if(Pid==2)
    {
    	hM=par[2]+1884.95;
	SiW=par[3]*sq(y/hM);
	Gee=0.26e-3;
        Mpsi=3769.9;
//	Bll=Gee/Gtot;
	Bll=0.0000112;
	Gtot=Gee/Bll;
	Gh=Gtot*(1.-2*Bll);
    }
    Double_t M=2*hM;
    beta=4*alfa/TMath::Pi()*(log(M/me)-0.5);
//    Finter=-Bll*R*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    Finter=-Bll*R/(1.-2.*Bll)*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=12*TMath::Pi()/sq(M)*Gee*Gh/Gtot/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    ParF[0]=beta;
    ParF[1]=M;
    ParF[2]=y*2;
    ParF[3]=SiW;
    ParF[4]=Gtot;
    ParF[5]=Finter;
    Double_t rIntegral=0;
#if 0
//    rIntegral=HANDLE_DGAUSS(FuncR,2*y-9*SiW,2*y+9*SiW,ParF,epsilon);
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-RangeE,2*y+RangeE,ParF,epsilon_);
#else
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-3.5,2*y+3.5,ParF,1.e-7)+
      HANDLE_DGAUSS(FuncR,2*y-7.,2*y-3.5,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y-20.,2*y-7.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncR,2*y+3.5,2*y+7.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y+7.,2*y+20.,ParF,1.e-3);
#endif
    xs=rIntegral/sqrt(2*TMath::Pi())/SiW*SS*par[1]+fabs(par[0])*sq(hM/y);
    return xs;
};
Double_t GrIntPsiPrime(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi,Gee=0,beta,DeltaE,Delta2,SS,MD;
    Double_t Gtot=0,Bll=0,R=2.5,Gh=0;
    Double_t ParF[6];
    Double_t xs,Finter;
    Mpsi=MPsiPrime;
    hM=par[0]+0.5*Mpsi;
    SiW=par[1]*sq(y/hM);
    Gee=2.12e-3;
    //	Gtot=0.277;
    //	Bll=Gee/Gtot;
    Bll=0.0090;
    Gtot=Gee/Bll;
    Gh=Gtot*(1.-2*Bll);
    Double_t M=2*hM;
    beta=4*alfa/TMath::Pi()*(log(M/me)-0.5);
//    Finter=-Bll*R*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    Finter=-Bll*R/(1.-2.*Bll)*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=12*TMath::Pi()/sq(M)*Gee*Gh/Gtot/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    ParF[0]=beta;
    ParF[1]=M;
    ParF[2]=y*2;
    ParF[3]=SiW;
    ParF[4]=Gtot;
    ParF[5]=Finter;
    Double_t rIntegral=0;
#if 0
//    rIntegral=HANDLE_DGAUSS(FuncR,2*y-9*SiW,2*y+9*SiW,ParF,epsilon);
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-RangeE,2*y+RangeE,ParF,epsilon_);
#else
    rIntegral=HANDLE_DGAUSS(FuncR,2*y-3.5,2*y+3.5,ParF,1.e-7)+
      HANDLE_DGAUSS(FuncR,2*y-7.,2*y-3.5,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y-20.,2*y-7.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncR,2*y+3.5,2*y+7.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncR,2*y+7.,2*y+20.,ParF,1.e-3);
#endif
    xs=rIntegral/sqrt(2*TMath::Pi())/SiW*SS;
    return xs;
};

Double_t RPsiDoublePrime(Double_t Eb,Double_t epsilon_,Double_t RangeE,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
//  Cross section in nb for given Ebeam, halfMass and SigmaW
//  par(0) - background cross section
//  par(1) - efficiency for given Gee
//  par(2) - mass shift
//  par(3) - Gtot

//  par(4) - Gee of psi' (relative)
//  par(5) - R0 in fm
//  par(6) - D nonresonant
//  par(7) - psi' mass shift
//  par(8) - sigmaW at psi'
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi,Gee=0,beta,DeltaE,Delta2,SS,S0,MD,betaD,sigmaE;
    Double_t Gtot=0,Bll=0,R=2.5,Gh=0;
    Double_t ParF[6];
    Double_t Finter,xs=0,xRes=0,xCont=0,xD=0;
    Double_t rIntegral=0;
    Double_t parPsiPrime[2]={par[7],par[8]};
    MD=0.5*(MDc+MD0);
    hM=par[2]+MPsiDoublePrime*0.5;
    SiW=par[8]*sq(2*y/MPsiPrime);
 //!   Gee=0.26e-3;
    Gee=GeePsiDoublePrime;
    //	Bll=Gee/Gtot;
    Bll=0.0000112;
    Gtot=par[3];

       //! Gee/Bll;
    Gh=Gtot*(1.-2*Bll);
    Double_t M=2*hM;
    beta=4*alfa/TMath::Pi()*(log(M/me)-0.5);
//  Finter=-Bll*R*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    Finter=-Bll*R/(1.-2.*Bll)*(1.+11./12.*beta)/(1.+0.75*beta)*2.*alfa/3.*sqrt(R*sq(Gtot)/Gee/Gh);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    // SS=12*TMath::Pi()/sq(M)*Gee*Gh/Gtot/M*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    sigmaE=par[8]/sqrt(2);
    S0=3*sqrt(TMath::Pi())*Gee*Gtot*sq(3.862*me/M)*1.e11/sigmaE;
    //    S0=3*sqrt(TMath::Pi())*0.276*0.001*23.6*sq(3.86*0.511/3770)*1.e11/sigmaE;
    ParF[0]=beta;
    ParF[1]=M;
    ParF[2]=y;
    ParF[3]=SiW;
    ParF[4]=Gtot;
    ParF[5]=par[5];
  //  cout<<"Gtot:"<<Gtot<<endl;
    if(Eb>(hM-25))
    {
//!#if 0
//    rIntegral=HANDLE_DGAUSS(FuncR,2*y-9*SiW,2*y+9*SiW,ParF,epsilon);
//    rIntegral=HANDLE_DGAUSS(FuncRDoublePrime,2*y-sigmaE*7,2*y+sigmaE*7,ParF,epsilon_);
    //  rIntegral=HANDLE_DGAUSS(FuncRDoublePrime,-sigmaE*7,sigmaE*7,ParF,0.00001);
      //#else
        rIntegral=HANDLE_DGAUSS(FuncRDoublePrime,-sigmaE*3.5,sigmaE*3.5,ParF,0.001)+
        HANDLE_DGAUSS(FuncRDoublePrime,-sigmaE*7,-sigmaE*3.5,ParF,0.005)+
        HANDLE_DGAUSS(FuncRDoublePrime,3.5*sigmaE,7*sigmaE,ParF,0.005);

//        HANDLE_DGAUSS(FuncRDoublePrime,-sigmaE*3.5,-sigmaE*7,ParF,0.001)+
  //          ;


/*       rIntegral=HANDLE_DGAUSS(FuncRDoublePrime,2*y-3.5,2*y+3.5,ParF,1.e-7)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y-7.,2*y-3.5,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y-20.,2*y-7.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y+3.5,2*y+7.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y+7.,2*y+20.,ParF,1.e-3);*//*
      rIntegral=HANDLE_DGAUSS(FuncRDoublePrime,2*y-3,2*y+35,ParF,1.e-7)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y-70.,2*y-35,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y-200.,2*y-70.,ParF,1.e-3)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y+35,2*y+70.,ParF,1.e-5)+
      HANDLE_DGAUSS(FuncRDoublePrime,2*y+70.,2*y+200.,ParF,1.e-3);*/
//!#endif
    }
    //  xs=rIntegral/sqrt(2*TMath::Pi())/SiW*SS*par[1]+fabs(par[0])*sq(hM/y);
    xRes=rIntegral*S0;
    xCont=fabs(par[0])*sq(MPsiPrime*0.5/y)+GrIntPsiPrime(y,epsilon_,12,parPsiPrime)*par[4];
    if(y>MD)
    {
	   betaD=sqrt(1-sq(MD/y));
           xD=par[6]*cub(betaD)*sq(MD/y);
    }
    else
    {
           xD=0;
    }
    xs=(xD+xRes+xCont)*par[1];
    return xs;
};
Double_t FuncR(Double_t W, Double_t* parf)
{
      Double_t Func;
      Double_t beta=parf[0];
      Double_t M=parf[1];
      Double_t Wb=parf[2];
      Double_t Sw=parf[3];
      Double_t Gtot=parf[4];
      Double_t Finter=parf[5];
      Double_t Fg=exp(-0.5*sq((W-Wb)/Sw));
      Double_t phi=TMath::ATan2(Gtot*0.5,M-W)*(1.-beta);
      if(phi<0) phi=phi+TMath::Pi();
      Func=pow((0.5*M/sqrt(sq(W-M)+sq(Gtot*0.5))),(1.-beta))*(sin(phi)+Finter*cos(phi))*Fg;
      return Func;
};

Double_t FuncRDoublePrime(Double_t dE, Double_t* parf)
{
      Double_t Func;
      Double_t beta=parf[0];
      Double_t M=parf[1];
      Double_t Ebeam=parf[2];
      Double_t Sw=parf[3];
      Double_t Gtot=parf[4];
      Double_t pDc=0,pD0=0,pD0R=0,pDcR=0;
      Double_t R0=parf[5]*1e-2/(3.862*me);
      Double_t delta=0,G=0,z=0,Gt=0;
      Double_t temp=sq(Ebeam+dE)-sq(MDc);
      Double_t pardg[1];
      pardg[0]=beta;
      if(temp<0) temp=0;
      pDc=sqrt(temp);
      temp=sq(Ebeam+dE)-sq(MD0);
      if(temp<0) temp=0;
      pD0=sqrt(temp);
      pD0R=sqrt(sq(MPsiDoublePrime)*0.25-sq(MD0));
      pDcR=sqrt(sq(MPsiDoublePrime)*0.25-sq(MDc));
      Gt=Gtot*(pow(pDc,3)/(1+sq(R0*pDc))+pow(pD0,3)/(1+sq(R0*pD0)))/
          (pow(pDc,3)/(1+sq(R0*pDcR))+pow(pD0R,3)/(1+sq(R0*pD0R)));
      delta=2*alfa/TMath::Pi()*(sq(TMath::Pi())/6-17./36.)+13./12.*beta;
      z=-2*dE/Sw;
      Double_t Fg=exp(-0.5*sq(z));
      G=Fzt_(&z,&beta)+delta*Fg;
      Func=G/(sq(2*(Ebeam+dE)-M)+sq(Gt)*0.25)*pow(2*Sw/M,beta)*sqrt(2*TMath::Pi())*Gt/Gtot;
      return Func;
};
Double_t RPsi(Double_t* Eb,Double_t* par)
{
    Double_t eps=0.000001;
    Double_t y=Eb[0];
    Double_t RangeE=100;
    Double_t xsec=RPsiDoublePrime(y,eps,RangeE,par);
    return  xsec;
}

Double_t HANDLE_DGAUSS(Double_t F(Double_t W,Double_t* parf),Double_t A,Double_t B,Double_t* par, Double_t eps)
{
    	const Double_t Z1 = 1;
  	const Double_t HF = Z1/2;
  	const Double_t CST= 5*Z1/1000;
  	Int_t ichanges=0;
  	Double_t x[12] = { 0.96028985649753623,  0.79666647741362674,
        		   0.52553240991632899,  0.18343464249564980,
                     	   0.98940093499164993,  0.94457502307323258,
                    	   0.86563120238783174,  0.75540440835500303,
                     	   0.61787624440264375,  0.45801677765722739,
                     	   0.28160355077925891,  0.09501250983763744};

  	Double_t w[12] = { 0.10122853629037626,  0.22238103445337447,
         		   0.31370664587788729,  0.36268378337836198,
                     	   0.02715245941175409,  0.06225352393864789,
                     	   0.09515851168249278,  0.12462897125553387,
                     	   0.14959598881657673,  0.16915651939500254,
			   0.18260341504492359,  0.18945061045506850};

      Double_t H=0,BB=0,AA=0,C1=0,C2=0,S8=0,S16=0,U=0,CONST=0;
      if(B==A){ goto end;}
      CONST=CST/fabs(B-A);
      BB=A;
  case1:     AA=BB; BB=B;
  case2:     C1=HF*(BB+AA); C2=HF*(BB-AA);
      S8=0;
      for(int i=0;i<4;i++)
      {
	      U=C2*x[i];
	      S8+=(w[i]*(F(C1+U,par)+F(C1-U,par)));
      }
      S16=0;
      for(int j=4;j<12;j++)
      {
	      U=C2*x[j];
 	      S16+=(w[j]*(F(C1+U,par)+F(C1-U,par)));
      }
      S16*=C2;
      if(fabs(S16-C2*S8) <= eps*(1+fabs(S16)))
      {
       H+=S16;
       if(BB != B) goto case1;
      }
      else
      {
       	BB=C1;
	if(1+CONST*fabs(C2)!=1) goto case2;
	if(ichanges<3)
	{
		cout<<"D103: TOO HIGH ACCURACY REQUIRED I'm trying change accuracy!"<<endl;
		cout<<"A = " << A << " B = " << B << " EPS = " << eps <<endl;
		eps*=1.25;
		ichanges++;
		goto case2;
	}
	else
	{
	    cout<<"D103: TOO HIGH ACCURACY REQUIRED ! I can't to change  accuracy !"<<endl;
	    H=0;
	}
       	goto end;
      }
end: 	return H;
};



long double HANDLE_DGAUSS2(long double F(long double W,long double* parf),long double A,long double B,long double* par, long double eps)
{
    	const long double Z1  = 1.;
  	const long double HF  = Z1/2.;
  	const long double CST = 5*Z1*0.001;
  	Int_t ichanges=0;
  	long double  x[12] = { 0.96028985649753623,  0.79666647741362674,
        		   0.52553240991632899,  0.18343464249564980,
                     	   0.98940093499164993,  0.94457502307323258,
                    	   0.86563120238783174,  0.75540440835500303,
                     	   0.61787624440264375,  0.45801677765722739,
                     	   0.28160355077925891,  0.09501250983763744};

  	long double  w[12] = { 0.10122853629037626,  0.22238103445337447,
         		   0.31370664587788729,  0.36268378337836198,
                     	   0.02715245941175409,  0.06225352393864789,
                     	   0.09515851168249278,  0.12462897125553387,
                     	   0.14959598881657673,  0.16915651939500254,
			   0.18260341504492359,  0.18945061045506850};

      long double  H=0,BB=0,AA=0,C1=0,C2=0,S8=0,S16=0,U=0,CONST=0;
      if(B==A){ goto end;}
      CONST=CST/fabs(B-A);
      BB=A;
  case1:     AA=BB; BB=B;
  case2:     C1=HF*(BB+AA); C2=HF*(BB-AA);
      S8=0;
      for(int i=0;i<4;i++)
      {
	      U=C2*x[i];
	      S8+=(w[i]*(F(C1+U,par)+F(C1-U,par)));
      }
      S16=0;
      for(int j=4;j<12;j++)
      {
	      U=C2*x[j];
 	      S16+=(w[j]*(F(C1+U,par)+F(C1-U,par)));
      }
      S16*=C2;
      if(fabs(S16-C2*S8) <= eps*(1+fabs(S16)))
      {
       H+=S16;
       if(BB != B) goto case1;
      }
      else
      {
       	BB=C1;
	if(1+CONST*fabs(C2)!=1) goto case2;
	if(ichanges<3)
	{
		cout<<"D103: TOO HIGH ACCURACY REQUIRED I'm trying change accuracy!"<<endl;
		cout<<"A = " << A << " B = " << B << " EPS = " << eps <<endl;
		eps*=1.25;
		ichanges++;
		goto case2;
	}
	else
	{
	    cout<<"D103: TOO HIGH ACCURACY REQUIRED ! I can't to change  accuracy !"<<endl;
	    H=0;
	}
       	goto end;
      }
end: 	return H;
};


Double_t xsecnbScansSeparatePsiPrime(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM,SiW,Mpsi,Gee,beta,DeltaE,Delta2,z,SS;
    Double_t xs;
    hM=par[2]+1843.00;
    SiW=par[3]*sq(y/hM);
    Gee=2.14e-3;
    Mpsi=3686.;
    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
  //  beta1=beta+1;
   // DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;
};

Double_t xsecnbScansSeparateBeta(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,beta1,DeltaE,Delta2,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }

    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
    beta1=beta+1;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*((1-par[4]*(2*y-Mpsi))*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW
    	+par[4]*beta/beta1*Fzt_(&z,&beta1)*pow(SiW/y,beta));
    return xs;
};
Double_t xsecnbScansSeparateErf(Double_t Eb,Double_t* par)
{
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,beta1,DeltaE,Delta2,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }

    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
    beta1=beta+1;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*((1-par[4]*(2*y-Mpsi))*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW
    	+par[4]*beta*0.5*(TMath::Erf(z/sqrt(2))+1));
    return xs;
};

Double_t xsecnbFit(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,DeltaE,Delta2,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }
    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
//    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
//    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
//    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
//    z=2*(y-hM)/SiW;
//    xs=fabs(par[0]*sq(hM/y))+par[1]*Fzt_(&z,&beta)*SS*pow(SiW/y,beta)/SiW;
 //   xs=par[0]+par[1]*SS*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    //  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;
};
Double_t xsecnbFitPsiPrime(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM,SiW,Mpsi,Gee,beta,DeltaE,Delta2,z,SS;
    Double_t xs;
    hM=par[2]+1843.00;
    SiW=par[3]*sq(y/hM);
    Gee=2.14e-3;
    Mpsi=3686.;
    beta=4*alfa/TMath::Pi()*(log(hM*2/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)/3-0.5)+0.75*beta;
    Delta2=1/24*sq(beta)*(2/3*log(hM*2/me)+2*sq(TMath::Pi())-37/4);
    SS=6*sq(TMath::Pi()/(hM*2))*Gee*sq(3.86*me)*1.e11*(1.+DeltaE-Delta2);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;
};
Double_t xsecnbFitNew(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,DeltaE,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }

    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+(1+par[4]*(y-Mpsi/2))*par[1]*Fzt_(&z,&beta)*SS*pow(SiW/y,beta)/SiW;
 //   xs=par[0]+par[1]*SS*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnbFitNew_(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];


    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,beta1,DeltaE,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }

    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
    beta1=beta+1;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*((1-par[4]*(2*y-Mpsi))*Fzt_(&z,&beta)*SS*pow(SiW/y,beta)/SiW
    +par[4]*beta/beta1*Fzt_(&z,&beta1)*SS*pow(SiW/y,beta));
//	+par[4]*beta*0.5*(TMath::Erf(z/sqrt(2))+1)*SS);
 //   xs=par[0]+par[1]*SS*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnb(Double_t* Eb,Double_t* par)
{
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,beta1,DeltaE,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM);
	Gee=2.14e-3;
	Mpsi=3686.;
    }
    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
    beta1=beta+1;
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0]*sq(hM/y))+par[1]*SS*Fzt_(&z,&beta)*pow(SiW/y,beta)/SiW;
    return xs;

};

Double_t xsecnbEI(Double_t Eb,Double_t I,Double_t* par){
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Mpsi=0,Gee=0,beta,DeltaE,z,SS;
    Double_t xs;
    if(Pid==0)
    {
    	hM=par[2]+1548.44;
	SiW=par[3]*(1+par[4]*(I-2.134))*sq(y/hM);
	Gee=5.26e-3;
	Mpsi=3096.88;
    }
    if(Pid==1)
    {
    	hM=par[2]+1843.00;
	SiW=par[3]*sq(y/hM)*(1+par[4]*(I-3.2));
	Gee=2.14e-3;
	Mpsi=3686.;
    }
    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0])+(1+par[5]*(y-Mpsi/2))*par[1]*Fzt_(&z,&beta)*SS*pow(2*SiW/Mpsi,beta)/SiW;
  //  xs=par[0]+par[5]*(I-2.08)+par[1]*Fzt_(&z,&beta)*SS*pow(2*SiW/Mpsi,beta)/SiW;
    //   xs=par[0]+par[1]*SS*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnbEIFix(Double_t* x,Double_t* par){
//     Routine arguments
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-29} cm^2]
    Double_t y=x[0];
    Double_t hM,SiW,Mpsi,Gee,beta,DeltaE,z,SS;
    Double_t xs;
    hM=par[2]+1548.44;
 //   hM=par[2]+1843;

    SiW=par[3]*(1+par[4]*(x[1]-2.1))*sq(y/hM);
 //   SiW=par[3];
    Mpsi=3096.88;
    Gee=5.26e-3;
    beta=4*alfa/TMath::Pi()*(log(Mpsi/me)-0.5);
    DeltaE=alfa/TMath::Pi()*(pow(TMath::Pi(),2)-0.5)+0.75*beta;
    SS=6*sq(TMath::Pi()/Mpsi)*Gee*sq(3.86*me)*1.e11*(1.+DeltaE);
    z=2*(y-hM)/SiW;
    xs=fabs(par[0])+par[1]*Fzt_(&z,&beta)*SS*pow(2*SiW/Mpsi,beta)/SiW;
 //   xs=par[0]+par[1]*SS*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnb_gaus(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full  xs [10^{-29} cm^2]
    Double_t y=Eb[0];
    Double_t hM,SiW,Mpsi,Gee,pi,beta,DeltaE,z;
    Double_t xs;
    hM=par[2];
    SiW=par[3];
    Mpsi=3096.88;
//    Mpsi=3686.00;
  //  Mpsi=par[2]*2.;
    Gee=5.26e-3;
    pi=3.14159265;
    beta=4*alfa/pi*(log(Mpsi/me)-0.5);
    DeltaE=alfa/pi*(pow(pi,2)-0.5)+0.75*beta;
    z=2*(y-hM)/SiW;
    xs=par[0]+par[1]*TMath::Gaus(y,hM,SiW)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnb_d_p3(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = sigma full
    Double_t y=Eb[0];
    Double_t hM,SiW,Mpsi,Gee,pi,beta,DeltaE,z;
    Double_t xs;
    hM=par[2];
    SiW=par[3];
   Mpsi=3096.88;
//    Mpsi=3686.00;
  //  Mpsi=par[2]*2.;
    Gee=5.26e-3;
    pi=3.14159265;
    beta=4*alfa/pi*(log(Mpsi/me)-0.5);
    DeltaE=alfa/pi*(pow(pi,2)-0.5)+0.75*beta;
    z=2*(y-hM)/SiW;
    xs=par[1]*Fzt_(&z,&beta)*2*(beta-1)*pow(2*SiW/Mpsi,beta)/sq(SiW)/Mpsi;
    xs+=par[1]*Fzt_d(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnb_0(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = full width  .
    Double_t y=Eb[0];
    Double_t hM,SiW,Mpsi,Gee,pi,beta,DeltaE,z;
    Double_t xs;
    hM=par[2]+1843;
    SiW=par[3];
   // Mpsi=3097.0;
    Mpsi=9460.0;
  //  Mpsi=par[2]*2.;
    Gee=5.26e-3;
    pi=3.14159265;
    beta=0;
    DeltaE=alfa/pi*(pow(pi,2)-0.5);
    z=2*(y-hM)/SiW;
    xs=par[0]+par[1]*Fzt_(&z,&beta)/SiW;
//  xs=par[0]+par[1]*Fzt_Aprox(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t xsecnb_d(Double_t* Eb,Double_t* par){
    //     Routine arguments
    //     par(0) = constant, par(1) = total cross section,
    //     par(2) = mass/2,   par(3) = width.
    Double_t y=Eb[0];
    Double_t hM,SiW,Mpsi,Gee,pi,beta,DeltaE,z;
    Double_t xs;
    hM=par[2] ;
    SiW=par[3] ;
    Mpsi=3097.0;
 // Mpsi=par[2]*2.;
    Gee=5.26e-3;
    pi=3.14159265;
    beta=4*alfa/pi*(log(Mpsi/me)-0.5);
    DeltaE=alfa/pi*(pow(pi,2)-0.5)+0.75*beta;
    z=2*(y-hM)/SiW;
    xs=par[1]*Fzt_d(&z,&beta)*pow(2*SiW/Mpsi,beta)/SiW;
    return xs;
};
Double_t SEvsI(Double_t icur,Double_t *par){
    Double_t se;
    se=par[0]+par[1]*(icur-par[2]);
    return se;
};
struct ESt{
    Int_t np;
    Double_t** Energy;
    Double_t** St;
};
class scan_R
{
    private:
//        Int_t ns;
        ESt* es;
       // Int_t* Nst;
       // Double_t tp;
       // Double_t lump;

    public:
        ~scan_R(){delete [] es;};
        scan_R(Int_t nnp);
        scan_R(Int_t nnp,Double_t *E,Double_t* S,Double_t* Ee,Double_t* Se);
        scan_R(Int_t nnp,Double_t *E,Double_t* S,Double_t Ee,Double_t Se);
        scan_R(Int_t nnp,Int_t* dimp,const scan_R&scan_Rold);
        scan_R(Int_t nnp,Double_t* EnLabel,scan_R& scan_Rold);
        void print_out();
        void GetDims(Double_t* E,Double_t* S,Double_t* Ee,Double_t* Se);
//        void GetAverDims(Int_t npd,Int_t* dimp,Double_t* E,Double_t* S,Double_t* Ee,Double_t* Se);
//        void Draw();
        void Draw(TCanvas* ScanG,Int_t cl);
//       void Fit(Double_t F_SvsE(Double_t* Eb,Double_t* par),TCanvas* ScanG,const char* NameFit,Int_t qf);
        void Fit(TCanvas* ScanG, const char* NameFit,Int_t qf);
        void Fit(const char* NameFit);
//        void ParSimCalc(Double_t Imax,Double_t Imin,Double_t tlive,Double_t l_0);
        Int_t GetNp(){return es->np;};
        Double_t  GetMinEnergy();
        Double_t  GetMaxEnergy();
};
Double_t scan_R::GetMinEnergy()
{
    Double_t EMin;
    EMin=es->Energy[0][0];
    for(int i=1;i<es->np;i++){
        if((es->Energy[0][i])<EMin) EMin=(es->Energy[0][i]);
    }
    return EMin;

};
Double_t scan_R::GetMaxEnergy(){
    Double_t EMax;
    EMax=es->Energy[0][0];
    for(int i=1; i<es->np;i++){
        if((es->Energy[0][i])>EMax) EMax=es->Energy[0][i];
    }
    return EMax;
};

scan_R::scan_R(Int_t nnp){
    es=new ESt ;
    es->np=nnp;
    es->Energy=new Double_t* [2];
    es->St=new Double_t* [2];
//    es->EorF=new Int_t [nnp];
    for(Int_t i=0;i<2;i++){
        es->Energy[i]=new Double_t [nnp];
        es->St[i]=new Double_t [nnp];
    }
};
scan_R::scan_R(Int_t nnp,Double_t *E,Double_t* S,Double_t* Ee,Double_t* Se){
        cout<<"num points  in scan  :"<<nnp<<endl;

    Double_t  EnMin;
    Double_t* EnergyOld;
    Double_t* EnergyNew;
    Int_t* ind;
    es=new ESt ;

    es->np=nnp;
    es->Energy=new Double_t* [2];
    es->St=new Double_t* [2];
    EnergyOld=new Double_t [nnp];
    EnergyNew=new Double_t [nnp];
    ind=new Int_t [nnp];
    cout<<"scan init ok"<<endl;
    for(Int_t i=0;i<2;i++)
    {
        es->Energy[i]=new Double_t [nnp];
        es->St[i]=new Double_t [nnp];
    }
    for(Int_t i=0;i<nnp;i++)
    {
        ind[i]=i;
        EnergyNew[i]=E[i];
    }
    for(int i=0;i<nnp;i++)
    {
        EnMin=EnergyNew[i];
        for(int j=i+1;j<nnp;j++){
            if(EnergyNew[j]<EnMin){
                swapI(ind[i],ind[j]);
                EnMin=EnergyNew[j];
                swapD(EnergyNew[i],EnergyNew[j]);

            }
        }
    }
    for(Int_t i=0;i<nnp;i++){
        es->Energy[0][i]=E[ind[i]];
        es->St[0][i]=S[ind[i]];
        es->Energy[1][i]=Ee[ind[i]];
        es->St[1][i]=Se[ind[i]];
    }

};
scan_R::scan_R(Int_t nnp,Double_t *E,Double_t* S,Double_t Ee,Double_t Se){
       cout<<"num points  in scan_  :"<<nnp<<endl;

    TRandom rndE,rndS;
    es=new ESt ;
    es->np=nnp;
    es->Energy=new Double_t* [2];
    es->St=new Double_t* [2];
    for(Int_t i=0;i<2;i++){
        es->Energy[i]=new Double_t [2];
        es->St[i]=new Double_t [2];
    }
    for(Int_t i=0;i<nnp;i++){
        es->Energy[0][i]=E[i];
        es->St[0][i]=S[i];
        es->Energy[1][i]=rndE.Gaus(E[i],Ee);
        es->St[1][i]=rndS.Gaus(S[i],Se);
    }
};

scan_R::scan_R(Int_t nnp,Int_t* dimp,const scan_R& scan_Rold){
    es=new ESt;
    es->np=nnp;
    es->Energy=new Double_t* [2];
    es->St=new Double_t* [2];
    for(Int_t j=0;j<2;j++){
        es->Energy[j]=new Double_t [nnp];
        es->St[j]=new Double_t [nnp];
    }
    aver_value(nnp,dimp,scan_Rold.es->Energy[0],scan_Rold.es->Energy[1],es->Energy[0],es->Energy[1]);
    aver_value(nnp,dimp,scan_Rold.es->St[0],scan_Rold.es->St[1],es->St[0],es->St[1]);
};

scan_R::scan_R(Int_t nnp,Double_t* EnLabel,scan_R& scan_Rold){
    Int_t* dimp;
    Int_t  l=0;
    Double_t tmp;
    es=new ESt;
    es->np=nnp;
    es->Energy=new Double_t* [2];
    es->St=new Double_t* [2];
    dimp=new Int_t [nnp];
    for(Int_t j=0;j<2;j++){
        es->Energy[j]=new Double_t [nnp];
        es->St[j]=new Double_t [nnp];
    }
    for(Int_t k=0;k<(nnp-1);k++){
        tmp=(EnLabel[k]+EnLabel[k+1])/2.;
        while(scan_Rold.es->Energy[0][l]<=tmp){
            dimp[k]++;
            l++;
        }
    }
    dimp[nnp-1]=scan_Rold.GetNp()-l;
    aver_value(nnp,dimp,scan_Rold.es->Energy[0],scan_Rold.es->Energy[1],es->Energy[0],es->Energy[1]);
    aver_value(nnp,dimp,scan_Rold.es->St[0],scan_Rold.es->St[1],es->St[0],es->St[1]);
};

void scan_R::print_out(){
    cout<<"\n np="<<es->np<<endl;
    for(Int_t i=0;i<es->np;i++){
        cout<<"\n  Energy["<<i<<"]="<<(es->Energy[0][i])<<endl;
        cout<<"\n  EnergyE["<<i<<"]="<<es->Energy[1][i]<<endl;
        cout<<"\n  St["<<i<<"]="<<es->St[0][i]<<endl;
        cout<<"\n  StE["<<i<<"]="<<es->St[1][i]<<endl;
    }
};

void scan_R::GetDims(Double_t* E,Double_t* S,Double_t* Ee,Double_t* Se)
{
    for(Int_t i=0;i<es->np;i++)
    {
        E[i]=es->Energy[0][i];
        Ee[i]=es->Energy[1][i];
        S[i]=es->St[0][i];
        Se[i]=es->St[1][i];
    }
};
void scan_R::Draw(TCanvas* ScanG,Int_t cl)
{
    TGraphErrors*  SvsE=new TGraphErrors(es->np,es->Energy[0],es->St[0],es->Energy[1],es->St[1]);
    SvsE->SetMarkerStyle(20);
    SvsE->SetMarkerColor(cl);
    ScanG->cd();
    SvsE->Draw("AP");
    ScanG->Update();
};
void scan_R::Fit(TCanvas* ScanG,const char* NameFit,Int_t qf)
{
    TGraphErrors*  SvsE=new TGraphErrors(es->np,es->Energy[0],es->St[0],es->Energy[1],es->St[1]);
// out<<"\n ok1"<<endl;

//SvsE->Fit(NameFit,"VRIM","O");
//cout<<"\n ok2"<<endl;
    //   SvsE->Fit(NameFit,"LIMREV");
//     SvsE->Fit(NameFit,"V","LIRE");

    if(qf==1)
    {
        SvsE->SetMarkerStyle(20);
        SvsE->SetMarkerColor(2);
	ScanG->cd();
	SvsE->Draw("AP");
	ScanG->Update();
	SvsE->Fit(NameFit,"LIMREQ");
	ScanG->Update();
    }
    if(qf==2)
    {
        SvsE->SetMarkerStyle(20);
        SvsE->SetMarkerColor(2);
//        ScanG->cd();
	SvsE->Draw("P");
//	ScanG->Update();
	SvsE->Fit(NameFit,"LIMREQ");
        ScanG->Update();
    }
//     SvsE->~TGraphErrors();
};
void scan_R::Fit(const char* NameFit)
{
    TGraphErrors*  SvsE=new TGraphErrors(es->np,es->Energy[0],es->St[0],es->Energy[1],es->St[1]);
    SvsE->Fit(NameFit,"NLIMREQ");

    SvsE->~TGraphErrors();
    /*
    if(qf){
        SvsE->SetMarkerStyle(20);
        SvsE->SetMarkerColor(2);
        ScanG->cd();
        ScanG->Clear();
        SvsE->Draw("AP");
        ScanG->Update();
        }
  */
};
/*
void ParSimCalc(Double_t Imax,Double_t Imin,Double_t tlive,Double_t l_0,){
    Double_t tp;
    Double_t lump;
    tp=tlive*ln(Imax/Imin);
    lump=l_0*0.5*tlive*(1-exp(-2*tp/tlive));
};


void GetAverDims(Int_t npd,Int_t* dimp,Double_t* E,Double_t* S,Double_t* Ee,Double_t* Se){
    aver_value(npd,dimp,Energy[0],Energy[1],E,Ee);
    aver_value(npd,dimp,St[0],St[1],S,Se);
};
*/

