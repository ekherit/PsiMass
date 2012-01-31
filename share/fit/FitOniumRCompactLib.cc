//--------------------------------------------------------------------------
// File and Version Information:
// FitOniumRCompactLib.cc
//  Description:
//  Library for fitting onium resonances (short version)
// Environment:
//      Software developed for the KEDR Detector at the VEPP-4M.
//      Branch for BES analysis
// Author List:
//      Korneliy Todyshev               Originator
// Copyright Information:
//      Copyright (C) 2000               KEDR
//
//------------------------------------------------------------------------

#include <iostream>
#include <math.h>   
#include <gsl/gsl_integration.h>
#include <complex>
#include <iomanip>
#include "FitTools/MathLibrary.h"
#include "FitOniumRCompactLib.hh"
using namespace std;
using std::setprecision;
double rangescale=1.5;
static const Double_t xgauss16[8] = {
    9.50125098376374401877e-02,    2.81603550779258913231e-01,
    4.58016777657227386350e-01,    6.17876244402643748452e-01,
    7.55404408355003033891e-01,    8.65631202387831743866e-01,
    9.44575023073232576090e-01,    9.89400934991649932601e-01
};
static const Double_t wgauss16[8] = {
    1.89450610455068496287e-01,    1.82603415044923588872e-01,
    1.69156519395002538183e-01,    1.49595988816576732080e-01,
    1.24628971255533872056e-01,    9.51585116824927848073e-02,
    6.22535239386478928628e-02,    2.71524594117540948514e-02
};
static const Double_t xgauss32[16] = {
    4.83076656877383162364e-02,    1.44471961582796493484e-01,
    2.39287362252137074544e-01,    3.31868602282127649782e-01,
    4.21351276130635345353e-01,    5.06899908932229390044e-01,
    5.87715757240762329066e-01,    6.63044266930215200960e-01,
    7.32182118740289680412e-01,    7.94483795967942406965e-01,
    8.49367613732569970160e-01,    8.96321155766052123971e-01,
    9.34906075937739689159e-01,    9.64762255587506430761e-01,
    9.85611511545268335400e-01,    9.97263861849481563534e-01
};

static const Double_t wgauss32[16] = {
    9.65400885147278005666e-02,    9.56387200792748594185e-02,
    9.38443990808045656367e-02,    9.11738786957638847129e-02,
    8.76520930044038111450e-02,    8.33119242269467552223e-02,
    7.81938957870703064685e-02,    7.23457941088485062287e-02,
    6.58222227763618468406e-02,    5.86840934785355471448e-02,
    5.09980592623761761959e-02,    4.28358980222266806557e-02,
    3.42738629130214331033e-02,    2.53920653092620594561e-02,
    1.62743947309056706058e-02,    7.01861000947009660028e-03
};

Double_t CrossSBhabhaPP(Double_t Eb,Double_t* par)
{
    Double_t y=Eb;
    Double_t hM=0;
    hM=1843.0;
    Double_t xsbb;
    xsbb=par[0]*sq(hM/y);
    return xsbb;
};

Double_t CrSOniumRAzimov(Int_t Id,Double_t Eb,Double_t* par)
{   
    Double_t y=Eb;
    Double_t hM=0,SiW=0,Gee=0,beta,DeltaE,Delta2,SS,FSS,MD;
    Double_t Gtot=0,Bll=0,Gh=0;
    Double_t  ParF[idNPar];
    Double_t xs,z,MR=0;
    Double_t Bh=0.0;
    if(Id==_IdJPsi){
      MR=_MJPsi;  
      Gee=_GeeJPsi;
      Bll=_BllJPsi;      
      SiW=par[idRSw]*sq(2.*y/_MJPsi);//1 
    }
    if(Id==_IdPsiPrime){
      MR=_MPsiPrime;  
      Gee=_GeePsiPrime;  
      SiW=par[idRSw]*sq(2.*y/_MPsiPrime); //1 
 
    }
    hM=par[idRM]+0.5*MR; //0     
    if(Id==_IdJPsi){
      Gtot=Gee/Bll;   
      Gh=Gtot*(1.-2.*Bll);
      Bh=1.-2.*Bll;
      ParF[idBhadr]=Bh;
    }

    if(Id==_IdPsiPrime){
      //      Gh=Gtot*(1.-3.*Bll);
      //      Gtot=0.286;//!!!!!
      Gee=_GeePsiPrime;
      Gtot=0.317;//!!!!!
      Bh=0.9785;// more precisely -> 0.13% //    
      ParF[idBhadr]=Bh;
    }     
    Double_t M=2.*hM;
    beta=4.*_alpha/TMath::Pi()*(log(2.*y/_me)-0.5);
    Double_t SSA=0;
    if(Id==_IdJPsi){
         SSA= 12.*TMath::Pi()*Gee*Bh;	 
    }
    else
      {
        SSA=Gee*Bh*_ConvConstant;	 
      }
   

    ParF[idbeta]=beta;
    ParF[idM]=M;
    ParF[idW]=y*2.;
    ParF[idSw]=SiW;
    ParF[idGt]=Gtot;
    ParF[ideff]=par[idReff];
    ParF[idefftau]=par[idRTauEff];
    ParF[idGee]=Gee; 
    ParF[idFreeInt]=par[idRFreeInt];
    ParF[idLambda]=par[idRLambda];

    Double_t RS=rangescale*SiW;
    Double_t rIntegral=0.0;
    if(Id==_IdPsiPrime){
      rIntegral= HANDLEDGAUSS(K_FuncRInterfPsiP,ParF[idW]-3.5*RS,ParF[idW]+3.5*RS,ParF,1.e-12)+
      HANDLEDGAUSS(K_FuncRInterfPsiP,ParF[idW]-7.0*RS,ParF[idW]-3.5*RS,ParF,1.e-12)+
      HANDLEDGAUSS(K_FuncRInterfPsiP,ParF[idW]-20.*RS,ParF[idW]-7.0*RS,ParF,1.e-12)+
      HANDLEDGAUSS(K_FuncRInterfPsiP,ParF[idW]+3.5*RS,ParF[idW]+7.0*RS,ParF,1.e-12)+
      HANDLEDGAUSS(K_FuncRInterfPsiP,ParF[idW]+7.0*RS,ParF[idW]+20.*RS,ParF,1.e-12);
      }
    else if(Id==_IdJPsi){
       rIntegral=         
         HANDLEDGAUSS(K_FuncRInterfJPsi,ParF[idW]-3.5*RS,ParF[idW]+3.5*RS,ParF,1e-8)+
         HANDLEDGAUSS(K_FuncRInterfJPsi,ParF[idW]-7.0*RS,ParF[idW]-3.5*RS,ParF,1.e-7)+
         HANDLEDGAUSS(K_FuncRInterfJPsi,ParF[idW]-20.*RS,ParF[idW]-7.0*RS,ParF,1.e-4)+
         HANDLEDGAUSS(K_FuncRInterfJPsi,ParF[idW]+3.5*RS,ParF[idW]+7.0*RS,ParF,1.e-7)+
         HANDLEDGAUSS(K_FuncRInterfJPsi,ParF[idW]+7.0*RS,ParF[idW]+20.*RS,ParF,1.e-4);
    }
    Double_t xRes=rIntegral*SSA/(SiW*sqrt(2.*TMath::Pi()));   
    return xRes;      
};

Double_t K_FuncRInterfPsiP(Double_t W, Double_t* parf)
{
 
  Double_t Bll=0.00772; //psi'
  Double_t Ratio=0;      
  Double_t  Rrat=2.14;
  
  Double_t  Bhadr=parf[idBhadr];    
  Double_t  beta=parf[idbeta];
  Double_t  M=parf[idM];
  Double_t  Wb=parf[idW];
  Double_t  Sw=parf[idSw];
  Double_t  Gtot=parf[idGt];   
  Double_t  GeeBhadr=parf[idGee];  // Gee*Bhadr    for fixed  efficiency
  Double_t  effh=parf[ideff];
  Double_t  efftt=parf[idefftau];
  Double_t  vy=sqrt(1.-sq(2*_MTau/W));  
  Double_t  Rtt=(3.-sq(vy))*0.5*vy*(1.-0.00015245952901473550*(0.5*W-1839.01));//simple approximation

     if(parf[idFreeInt]==1)
      { 
        Ratio=fabs(parf[idLambda]);
      }
     else
       {
//       Ratio=(sqrt(Rrat*Bll)+efftt/effh*sqrt(Rtt*_BttPP/Rrat))/sqrt(0.9785+efftt/effh*_BttPP);
//fixed lambda
         Ratio=(sqrt(Rrat)+efftt/effh*sqrt(Rtt*_BttPP/Rrat/Bll))/sqrt(1.+efftt/effh*_BttPP/0.9785);
       }    
      Double_t  Fg=1;    
#ifdef _ParChromOn_
      Double_t mchrom=1.+parf[idChrom]*(W-Wb);
      if(mchrom<0) mchrom=0;
      Fg=exp(-0.5*sq((W-Wb)/Sw))*mchrom;      
#else
      Fg=exp(-0.5*sq((W-Wb)/Sw));
#endif
      Double_t  betaW= 4.*_alpha/TMath::Pi()*(log(W/_me)-0.5);   
      Double_t  Delta1=_alpha/TMath::Pi()*(pow(TMath::Pi(),2.)/3.-0.5)-1./24.*sq(betaW)*(2./3.*log(W/_me)+2.*sq(TMath::Pi())-37./4.);
      Double_t  DeltaE=Delta1+0.75*betaW;          
      complex<double> f=0;      
      complex<double> fp=0;      
      Double_t  PiHre=-0.6e-2;//  ! Berends, Kommen for W from 2 to 4 GeV
      complex<double> Pi10C(1.-(PiLre(W,_me)+PiLre(W,_mmu)+PiLre(W,_MTau)+PiHre),(PiLim(W,_me)+PiLim(W,_mmu)+PiLim(W,_MTau)));         
      complex<double> M2(sq(M/W)-1.,-Gtot*M/W/W);
      complex<double> beta_11(betaW-1.,0.);
      f=pow(M2,beta_11);
      fp=pow(M2,beta_11)/Pi10C;             

      Double_t rI=0.;
      Double_t s=W*W;
      Double_t M2_S=sq(M/W);
      Double_t delta=(Gtot*M)/s;        

 
      rI= 12.*TMath::Pi()*Fg/s/M*effh*(
                   imag(f)*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)
                   -betaW*0.5*((atan(M/Gtot)-atan((M-W*W/M)/Gtot))*(1.+sq(M/W)))
                   +0.25*betaW*Gtot*M/W/W*(log((sq(M2_S)+sq(delta))/(sq(M2_S-1.)+sq(delta))))                                      
//             -2./3.*_alpha*sqrt(Rrat*Gtot/GeeBhadr)*Ratio*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)*M/W*real(fp)

            //                                           free lambda 
                     //         -2./3.*_alpha*sqrt(Rrat*Gtot/GeeBhadr)*Ratio*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)*M/W*real(fp)                    
                     //fixed
                   -2./3.*_alpha*sqrt(Rrat)/0.9785*Ratio*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)*M/W*real(fp)					   );

      

  return rI;
};


Double_t K_FuncRInterfJPsi(Double_t W, Double_t* parf)
{
 
      Double_t  Bll=_BllJPsi; 
      Double_t  Ratio=0;      
      Double_t  Rrat=2.14;
      Double_t  Bhadr=parf[idBhadr];    
      Double_t  beta=parf[idbeta];
      Double_t  M=parf[idM];
      Double_t  Wb=parf[idW];
      Double_t  Sw=parf[idSw];
      Double_t  Gtot=parf[idGt];   
      Double_t  GeeBhadr=parf[idGee];  // Gee*Bhadr    for fixed  efficiency 
      Double_t  effh=parf[ideff];


     if(parf[idFreeInt]==1)
      { 
        Ratio=fabs(parf[idLambda]);
      }
     else
       {
//         Ratio=sqrt(Rrat*Bll/(1.-_BllPPSum)); 
//fixed lambda
         Ratio=sqrt(Rrat);
       }    
      Double_t  Fg=1.0;
     

#ifdef _ParChromOn_
      Double_t mchrom=1.+parf[idChrom]*(W-Wb);
      if(mchrom<0) mchrom=0;
      Fg=exp(-0.5*sq((W-Wb)/Sw))*mchrom;      
#else
      Fg=exp(-0.5*sq((W-Wb)/Sw));
#endif
      Double_t  betaW= 4.*_alpha/TMath::Pi()*(log(W/_me)-0.5);   
      Double_t  Delta1=_alpha/TMath::Pi()*(pow(TMath::Pi(),2.)/3.-0.5)-1./24.*sq(betaW)*(2./3.*log(W/_me)+2.*sq(TMath::Pi())-37./4.);
      Double_t  DeltaE=Delta1+0.75*betaW;          
      complex<double> f=0;      
      complex<double> fp=0;      
      Double_t  PiHre=-0.6e-2;//  ! Berends, Kommen for W from 2 to 4 GeV
      complex<double> Pi10C(1.-(PiLre(W,_me)+PiLre(W,_mmu)+PiHre),(PiLim(W,_me)+PiLim(W,_mmu)));   
      
#ifdef RELATIVE	 
      complex<double> M2(sq(M/W)-1.,-Gtot*M/W/W);
      complex<double> beta_11(betaW-1.,0.);
      f=pow(M2,beta_11);
      fp=pow(M2,beta_11)/Pi10C;
      //fp=pow(M2,beta_11);         
#else 
         complex<double> PW_M_G(M/W+1.,0.5*Gtot/W);
         complex<double> MW_M_G(M/W-1.,-0.5*Gtot/W);
         complex<double> M_G(M/W,0.5*Gtot/W);      
         complex<double> beta_1(betaW-1.,0.);      
         f=pow(PW_M_G,beta_1)*pow(MW_M_G,beta_1)*M_G;   
#endif 	 
         Double_t rI=0.;
         Double_t s=W*W;
         Double_t M2_S=sq(M/W);
         Double_t delta=(Gtot*M)/s;        
        
         
          rI= Fg/s/M*effh*(
                         
				   
                         imag(f)*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)
                         -betaW*0.5*((atan(M/Gtot)-atan((M-W*W/M)/Gtot))*(1.+sq(M/W)))
                         +0.25*betaW*Gtot*M/W/W*(log((sq(M2_S)+sq(delta))/(sq(M2_S-1.)+sq(delta))))                                        				   
                    
                       //                                           free lambda 
                         ///        -2./3.*_alpha*sqrt(Rrat*Gtot/GeeBhadr)*Ratio*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)*M/W*real(fp)                    
                         //fixed
                         -2./3.*_alpha*sqrt(Rrat)/0.877*Ratio*TMath::Pi()*betaW/sin(TMath::Pi()*betaW)*(1.+DeltaE)*M/W*real(fp)
                                                                                                                                                                                  
                         );
          
         

         
          
return rI*_ConvConstant;
};







Double_t FCrSPPrimeAzimov(Double_t* Eb,Double_t* par)
{
    Double_t y=Eb[0]/(1.0*ScaleEGr);
    Double_t xsec=CrSOniumR(1, _IdPsiPrime ,y,par);
    //    if(y>1842.9&&y<1843.1) cout<<"Eb:"<< y<< " s:"<<xsec<<endl;   
    return  xsec;
};

Double_t FCrSJpsiAzimov(Double_t* Eb,Double_t* par)
{
    Double_t y=Eb[0]/(1.0*ScaleEGr);
    Double_t xsec=CrSOniumR(1,_IdJPsi,y,par);
    //    if(y>1842.9&&y<1843.1) cout<<"Eb:"<< y<< " s:"<<xsec<<endl;   
    return  xsec;
};

Double_t HANDLEDGAUSS(Double_t F(Double_t W,Double_t* parf),Double_t A,Double_t B,Double_t* par, Double_t epsc)
{
  Int_t ichanges=0;  
  Double_t H=0,BB=0,AA=0,C1=0,C2=0,S16=0,S32=0,S64=0,U=0,CONST=0;
  Double_t eps=epsc;
  if(epsc<0){
    cout<<"HDG error epsc:"<<epsc<<endl;
  }
  if(B<=A){goto end;}
  CONST=0.005/fabs(B-A);
  BB=A;
case1:     AA=BB; BB=B;
case2:     C1=0.5*(BB+AA); C2=0.5*(BB-AA);   
           S16=wgauss16[0]*(F(C1+C2*xgauss16[0],par)+F(C1-C2*xgauss16[0],par))+
             wgauss16[1]*(F(C1+C2*xgauss16[1],par)+F(C1-C2*xgauss16[1],par))+
             wgauss16[2]*(F(C1+C2*xgauss16[2],par)+F(C1-C2*xgauss16[2],par))+
             wgauss16[3]*(F(C1+C2*xgauss16[3],par)+F(C1-C2*xgauss16[3],par))+
             wgauss16[4]*(F(C1+C2*xgauss16[4],par)+F(C1-C2*xgauss16[4],par))+
             wgauss16[5]*(F(C1+C2*xgauss16[5],par)+F(C1-C2*xgauss16[5],par))+
             wgauss16[6]*(F(C1+C2*xgauss16[6],par)+F(C1-C2*xgauss16[6],par))+
             wgauss16[7]*(F(C1+C2*xgauss16[7],par)+F(C1-C2*xgauss16[7],par));  

           S32=wgauss32[0]*(F(C1+C2*xgauss32[0],par)+F(C1-C2*xgauss32[0],par))+
             wgauss32[1]*(F(C1+C2*xgauss32[1],par)+F(C1-C2*xgauss32[1],par))+
             wgauss32[2]*(F(C1+C2*xgauss32[2],par)+F(C1-C2*xgauss32[2],par))+
             wgauss32[3]*(F(C1+C2*xgauss32[3],par)+F(C1-C2*xgauss32[3],par))+
             wgauss32[4]*(F(C1+C2*xgauss32[4],par)+F(C1-C2*xgauss32[4],par))+
             wgauss32[5]*(F(C1+C2*xgauss32[5],par)+F(C1-C2*xgauss32[5],par))+
             wgauss32[6]*(F(C1+C2*xgauss32[6],par)+F(C1-C2*xgauss32[6],par))+
             wgauss32[7]*(F(C1+C2*xgauss32[7],par)+F(C1-C2*xgauss32[7],par))+
             wgauss32[8]*(F(C1+C2*xgauss32[8],par)+F(C1-C2*xgauss32[8],par))+
             wgauss32[9]*(F(C1+C2*xgauss32[9],par)+F(C1-C2*xgauss32[9],par))+
             wgauss32[10]*(F(C1+C2*xgauss32[10],par)+F(C1-C2*xgauss32[10],par))+
             wgauss32[11]*(F(C1+C2*xgauss32[11],par)+F(C1-C2*xgauss32[11],par))+
             wgauss32[12]*(F(C1+C2*xgauss32[12],par)+F(C1-C2*xgauss32[12],par))+
             wgauss32[13]*(F(C1+C2*xgauss32[13],par)+F(C1-C2*xgauss32[13],par))+
             wgauss32[14]*(F(C1+C2*xgauss32[14],par)+F(C1-C2*xgauss32[14],par))+
             wgauss32[15]*(F(C1+C2*xgauss32[15],par)+F(C1-C2*xgauss32[15],par));  

           S32*=C2;    
           if(fabs(S32-C2*S16) <= eps*(1.+fabs(S32)))
           {
             H+=S32;
             if(BB != B) goto case1;
           }     
           else
           {      
             BB=C1;
             if(1.+CONST*fabs(C2)!=1.) goto case2;
             H=S16;
             //if(ichanges<3)
             //{
             //  cout<<"D103: TOO HIGH ACCURACY REQUIRED I'm trying change accuracy!"<<endl;
             //  cout<<"A = " << A << " B = " << B << " EPS = " << eps <<endl;
             //  eps*=10.;
             //  ichanges++;
             //  goto case2;
             //}
             //else
             //{
             //  cout<<"D103: TOO HIGH ACCURACY REQUIRED ! I can't to change  accuracy !"<<endl;
             //  H=0;
             //}
             goto end;
           }

end: 	return H;
};



//Tools for cross section
Double_t FCrossSection(Double_t* Eb,Double_t* par)
{
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-33} cm^2]
  double y=Eb[0]/ScaleEGr;
  double cross=CrSOniumR(_MethodAzimov,_IdJPsi,y,par);
   // double cross=CrSOniumR(_MethodDIntegral,_IdJPsi,y,par);
  //   double cross=CrSOniumR(_MethodAzimov,_IdPsiPrime,y,par);
  // double cross=CrSOniumR(_MethodDIntegral,_IdPsiPrime,y,par);
   //double cross=TMath::Pi()*_alpha*_alpha/3./y/y*_ConvConstant*2.2;  
   cout<<"Eb[0]:"<<Eb[0]<<" cross:"<<cross<<endl;
  // double cross=3.*86.8/4./y/y*1000000.;  
  return cross;
}

Double_t myJPsiCrossSection(Double_t* Eb,Double_t* par)
{
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-33} cm^2]
  double y=Eb[0];
  double cross=CrSOniumR(_MethodAzimov,_IdJPsi,y,par);
  return cross;
}

Double_t myPsiPrimeCrossSection(Double_t* Eb,Double_t* par)
{
//     par(0) = constant, par(1) = total cross section,
//     par(2) = mass/2,   par(3) = sigma full  //xs [10^{-33} cm^2]
  double y=Eb[0];
  double cross=CrSOniumR(_MethodAzimov,_IdPsiPrime,y,par);
  return cross;
}

Double_t CrSOniumR(Int_t Method,Int_t Id,Double_t Eb,Double_t* par)
{
  Double_t y=Eb;
  Double_t hM=0,SiW=0,Gee=0,beta,DeltaE,Delta2,SS,MD;
  Double_t Gtot=0,Bll=0,R=2.5,Gh=0;
  Double_t xs,MR=0;     
  Double_t rIntegral=0;
  Double_t fixenergy=1843.;
  if(Id==_IdJPsi)
  {
    fixenergy=1548.;
  }
  if(Id==_IdPsiPrime)
  {
    fixenergy=1843.;
  }
  rIntegral=CrSOniumRAzimov(Id, y, par);          
  xs= fabs(par[idRbg]*sq(fixenergy/y))+par[idReff]*rIntegral;             
  return xs;
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

int GetNumRows(const char *FileName,int npar)
{
  int counter=-1;
  Double_t Spool;
  ifstream test(FileName,ios::in);//|ios::floatfield);
  if(!test){
    cout<<"there is no file:"<<FileName<<endl;
    counter=0;
    test.close();
  }
  else
  {
    // counter=0;     
    while(test.get()!=EOF) // find the end off file  runs.par
    {
      counter++;             // number of lines
      for(int j=0;j<npar;j++)
      {
        test>>Spool;
      }
    }
    test.close();
  } 
  return counter;
}

void FillArrayFromFile(const char* FileName,Double_t** Array,int npar,int nps)
{
  ifstream readingfile(FileName,ios::in);
  if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
  else
  {
    for(int i=0;i<nps;i++)
    {
      Array[i]=new Double_t [npar];
      for(int j=0;j<npar;j++)
      {
        readingfile>>Array[i][j];
      }
    }
    cout<<endl;
  }
  readingfile.close();
}

void FillArrayFromFile(const char* FileName,Double_t** Array,int npar,int nparNew,int nps)
{
    ifstream readingfile(FileName,ios::in);
    
     if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
     else
     {
     for(int i=0;i<nps;i++)
     {
     	Array[i]=new Double_t [nparNew];
     	for(int j=0;j<npar;j++)
	{
	    readingfile>>Array[i][j];
	}
     }
     cout<<endl;
     }
     readingfile.close();
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
Double_t  PiLre(Double_t W,Double_t M){
  Double_t x=sq(2*M/W);
  Double_t sqr=0;
  Double_t PiLre=0;
  if(x<=1) {
    sqr=sqrt(1.-x);
    PiLre=-5./9.-x/3.+(2.+x)/6.*sqr*log(fabs((1.+sqr)/(1.-sqr)));
  }
  else{
    sqr=sqrt(x-1.);
    PiLre=-5./9.-x/3.+(2.+x)/6.*sqr*atan(1./sqr);
  }
  return PiLre*_alpha/M_PI;   
}


Double_t PiLim(Double_t W,Double_t M)
{
  Double_t  x=sq(2.*M/W);
  Double_t  phi=0;
  if(x<=1) {
    phi=(2.+x)/6.*sqrt(1.-x);
  }
  else
    {
      phi=0;
    }    
  return -(1./3.+phi)*_alpha;
}
