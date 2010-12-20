#define _HESSE_ 1 // calculate errors with HESSE
#define _MINOs_ 1 // use minos errors
//#include"KdDCSim/curve.hh"

#include<iostream>
#include<stdlib.h>
#include<iomanip>
#include<math.h>
#include <complex.h>
#include<string>
#include<assert.h>
#include<stdio.h>
#include<fstream>
#include<unistd.h>
#include <error.h>
#include <argp.h>
#include <TROOT.h>
#include <TH1.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TAxis.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TClonesArray.h>
#include <TPad.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFormula.h>
#include <TString.h>
#include <TVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TPostScript.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TProfile.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRandom.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
R__EXTERN TSystem *gSystem;
extern void InitGui();
VoidFuncPtr_t initfuncs[] = { InitGui, 0 };
TROOT root("ROOTPROGRAM","ROOTPROGRAM", initfuncs);
#include"FitOniumRCompactLib.hh"
#include<TCanvas.h>
#include"TPostScript.h"
#include<TApplication.h>
using namespace std;
//Double_t epsilon=0.000001;
//#include "FitTools/jpsiM.h"
#define  NumMaxP 120
Double_t EMaxRange=1550.;
Double_t EMinRange=1530.;
Int_t    JpsiFitOnly=0;
Int_t    NpPP=0;
Int_t     UseLumBB=0;
Double_t CrossBhabha=265;

Double_t CrossSInScan[NumMaxP];
Double_t CrossSErrInScan[NumMaxP];
Double_t CrossSBBInScan[NumMaxP];
Double_t CrossSBBErrInScan[NumMaxP];
Double_t LumInScan[NumMaxP];
Double_t NmhInScan[NumMaxP];
Double_t LumLgammaInScan[NumMaxP];
Double_t NbbInScan[NumMaxP];
Double_t EInScan[NumMaxP];
Double_t WInScan[NumMaxP];
Double_t EErrInScan[NumMaxP];
Double_t WErrInScan[NumMaxP];
Int_t    NumEpoints=0;
Double_t MinChi2=1e+7;
Int_t   FreeEnergyFit=0;
//id fcnResChi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
/* Program documentation. */

static char doc[] ="Поиск узких резонансов";

/* A description of the arguments we accept. */
static char args_doc[] = "";
#define OPT_MStep             1
#define OPT_ERangeL           2
#define OPT_ERangeR           3
#define OPT_NeeFlag           4
#define OPT_NmhFlag           5
#define OPT_CrBhabha          6
#define OPT_dEmin             7
#define OPT_Emin              8
#define OPT_Emax              9
#define OPT_scan             10

#define OPT_sel             11
#define OPT_lumbb           12




static struct argp_option options[] = {
    {"verbose",'v',0,0,"produce verbose information",100},
    {"quickly",'q',0,0," fast mode",100},
    {"FreeEnergy",'E',0,0,"FreeEnergy",100},
    {"JpsiFit",'J',0,0," fit Jpsi",100},
    {"Chi2",'C',0,0," chi2",100},
    {"MStep",OPT_MStep,"MeV",0,"step",1},
    {"ERangeL",OPT_ERangeL,"MeV",0,"ERangeL",2},
    {"ERangeR",OPT_ERangeR,"MeV",0,"ERangeR",3},
    {"NeeFlag",OPT_NeeFlag,"sel",0,"NeeFlag",4}, 
    {"NmhFlag",OPT_NmhFlag,"sel",0,"NmhFlag",5}, 
    {"CrBhabha",OPT_CrBhabha,"nb",0,"CrBhabha",6}, 
    {"dEmin",OPT_dEmin,"MeV",0,"dEmin",7}, 
    {"Emin",OPT_Emin,"MeV",0,"Emin",8}, 
    {"Emax",OPT_Emax,"MeV",0,"Emax",9}, 
    {"scan",OPT_scan,"scan",0,"scan",10}, 
    {"sel",OPT_sel,"sel",0,"sel",11}, 
    {"lumbb",OPT_lumbb,"lumbb",0,"lumbb",12},     
    
   
    {0}
};

struct arguments
{
  int verbose;
  int quickly;
  int FreeEnergy;
  int JpsiFit;
  double MStep;
  double ERangeL;
  double ERangeR;
  int    NeeFlag;
  int    NmhFlag;
  double CrBhabha;
  double dEmin;
  double Emin;
  double Emax;
  int Chi2;

  int sel;
  int lumbb;
  int scan;

};
static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    union ArG
    {
    	struct arguments *arguments;
    	void *void_arg;
    };
    ArG arg_union;
    arg_union.void_arg=state->input;
    switch (key)
    {
    case 'v': arg_union.arguments->verbose += 1; break;
    case 'q': arg_union.arguments->quickly += 1; break;
    case 'J': arg_union.arguments->JpsiFit += 1; break;
    case 'E': arg_union.arguments->FreeEnergy += 1; break;
    case 'C': arg_union.arguments->Chi2 += 1; break;
    case OPT_MStep:arg_union.arguments->MStep=atof(arg);break;
    case OPT_CrBhabha:arg_union.arguments->CrBhabha=atof(arg);break;
    case OPT_ERangeL:arg_union.arguments->ERangeL =atof(arg);break;
    case OPT_ERangeR:arg_union.arguments->ERangeR =atof(arg);break;
    case OPT_NeeFlag:arg_union.arguments->NeeFlag=atoi(arg);break;
    case OPT_NmhFlag:arg_union.arguments->NmhFlag=atoi(arg);break;            
    case OPT_dEmin:arg_union.arguments->dEmin=atof(arg);break;            
    case OPT_Emin:arg_union.arguments->Emin=atof(arg);break;            
    case OPT_Emax:arg_union.arguments->Emax=atof(arg);break;            
    case OPT_scan:arg_union.arguments->scan=atoi(arg);break;            
    case OPT_sel:arg_union.arguments->sel=atoi(arg);break;            
    case OPT_lumbb:arg_union.arguments->lumbb=atoi(arg);break;            
         
  
    }
    return 0;
};
static struct argp argp = { options, parse_opt, args_doc, doc };


int main(int argc, char **argv)
{

  struct arguments arguments;
  arguments.verbose=0;
  arguments.quickly=0;
  arguments.scan=1;
  arguments.sel=0;
  arguments.MStep=  0.1;
  arguments.ERangeL= 6.5;
  arguments.ERangeR= 6.5;
  arguments.NeeFlag= 0;
  arguments.NmhFlag= 0;
  arguments.FreeEnergy= 0;
  arguments.JpsiFit= 1;
  arguments.Emin= 500.;
  arguments.Emax= 1600.;  
  arguments.Chi2= 0;
  arguments.lumbb= 0;
  //arguments.CrBhabha=62;  
  // arguments.CrBhabha=535;  
 
  arguments.CrBhabha=80; 
 
  arguments.dEmin=0.02 ;  
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  TApplication* theApp=0;

  cout<<"JpsiFit:"<<arguments.JpsiFit<<"Ch2:"<<arguments.Chi2<<endl;
  JpsiFitOnly=arguments.JpsiFit;
  UseLumBB=arguments.lumbb;
  CrossBhabha= arguments.CrBhabha;
  FreeEnergyFit= arguments.FreeEnergy;
  if(arguments.quickly==0)
    {
      theApp=new  TApplication("App", &argc, argv);
    };  
  TMinuit*  MinuitRes=0; 
  TF1* FitRes=0;
  TF1* FitResBG=0;
  TGraphErrors* GrRes=0;
  TLine* LineRes=0;
  int dimMHFile=7;
  int MHRun=0;
  int MHLum=1;
  int MHEnergy=2;
  int MHEnergyErr=3;
  int MHNee=5;
  int MHNgg=6; 
  int MHNh=4; 
  double** AllMH=0;
  int      npMHFile=0; 
  if(arguments.scan==1){
    npMHFile=GetNumRows("scan1.txt",dimMHFile);       
  }
  else if(arguments.scan==2)
    {
      npMHFile=GetNumRows("scan2.txt",dimMHFile);          
    }
  else if(arguments.scan==3){
    npMHFile=GetNumRows("scan12.txt",dimMHFile); 
  }
  AllMH=new double* [npMHFile];
  if(arguments.scan==1){
    FillArrayFromFile("scan1.txt",AllMH,dimMHFile,npMHFile);  
  }
  else if(arguments.scan==2)
    {
      FillArrayFromFile("scan2.txt",AllMH,dimMHFile,npMHFile);  

    }
  else if(arguments.scan==3){
    FillArrayFromFile("scan12.txt",AllMH,dimMHFile,npMHFile);  
  }
  cout<<"npMHFile:"<<npMHFile<<endl;
  
  int dimAP= 15;
  int ARun=0;
  int AEnergy =1;
  int AEnergyErr=2;
  int AMHEv  =3;
  int AMHEvH =4;
  int AEE    =5;
  int ALe    =6;
  int ALp    =7;
  int ACr    =8;
  int ACrErr =9;
  int ALRatE  =10;
  int ALRatP  =11;
  int ALRatEErr  =12;
  int ALRatPErr  =13;
  int ASE  =14;  
  int npAP=0;    
  npAP=npMHFile;
  double** AP=new double* [npAP];
  int Aind=0;
   
  for(int i=0;i<npMHFile;i++ )
    {      
      AP[Aind]=new double [dimAP];    
      AP[Aind][ARun]=AllMH[i][MHRun];                    
      AP[Aind][AEnergy]=AllMH[i][MHEnergy]*0.5;        
      AP[Aind][AEnergyErr]=AllMH[i][MHEnergyErr]*0.5;       
      AP[Aind][ALe]=AllMH[i][MHLum];
      AP[Aind][ALp]=AllMH[i][MHLum];
      AP[Aind][AMHEv]=AllMH[i][MHNh];
      AP[Aind][AEE]=AllMH[i][MHNee];            
      if(arguments.NeeFlag==1){
        AP[Aind][AEE]=AllMH[i][MHNgg];
      }     
      Aind++;
    }
    
  
  
  if(Aind!=npAP) cout<<"PROBLEM !!!!!!!!!"<<endl;
  int   NumParForC[3];  
  int   CUpDown[3];
  NumParForC[0] = AEnergy;
  CUpDown[0]  = 1;
  NumParForC[1] = ARun;
  CUpDown[1]  = 1;
 
  for(Int_t i=0;i<npAP;i++) // simplest bubble's method for each condition
    {
      for(Int_t j=0;j<i;j++)
        {
          for(Int_t ic=1;ic>=0;ic--)
            {
      
              if(comparDRows(CUpDown[ic],AP[i],AP[j],NumParForC[ic])<0)
                {
                  swapDRows(AP[i],AP[j],dimAP);

                }
            }
        }
    }
    
  char ouputstring[140];
  ofstream AA("test.txt", ios::out);
  for(int i=0;i<npAP;i++)
    {
      sprintf(ouputstring,"  %7.4f %.0f %.0f %.0f",AP[i][AEnergy],AP[i][ARun],AP[i][AMHEv],AP[i][AEE]); 
      AA<<ouputstring<<"\n";	            
    }
  AA.close();        
  Double_t EMin=AP[0][AEnergy]-1 ;      

  Double_t EMax=AP[npAP-1][AEnergy]+1;
  //  EMin=Emin-6.0;
  //  EMax=1843+5.0;
  //  double   Energy=EMin; 
  //uble   Energy=EMin; 




   
 
 
 
 
     
  int np=npAP ;   
  cout<<"EMin:"<<EMin<<"EMax:"<<EMax<<"np:"<<np<<endl;
   
  Double_t *En_=new Double_t[np];
  Double_t *Eerr_=new Double_t[np];
  Double_t *Nmh_=new Double_t[np];
  Double_t *Nbb_=new Double_t[np];
  Double_t *Le_=new Double_t[np];
  Double_t *Lp_=new Double_t[np];    
  Int_t*        Euse=new Int_t[np];
  np=0;
  
  Double_t SigmaW=1.0;
  for(int i=0;i<npAP;i++)
    {   
      En_[np]=AP[i][AEnergy];
      Eerr_[np]=AP[i][AEnergyErr];
      Nmh_[np]=AP[i][AMHEv];
      Nbb_[np]=AP[i][AEE];  
      Le_[np]=AP[i][ALe];  
      Lp_[np]=AP[i][ALp];  
      np++;      
    }
  int  NEp=0;
  SeparatePointsPartNew(np,Nbb_,&NEp,Euse,En_,arguments.dEmin);      
  cout<<"!!!dEmin:"<<arguments.dEmin<<endl;
  Double_t *En  =new Double_t[NEp];
  Double_t *Eerr=new Double_t[NEp];
  Double_t *Nmh=new  Double_t[NEp];
  Double_t *Nbb=new  Double_t[NEp];
  Double_t *Le=new  Double_t[NEp];
  Double_t *Lp=new  Double_t[NEp];
  
  SumPointsByQuant(np,NEp,Euse,En_,Nbb_,Eerr_,En,Eerr,false);
  SumPointsSimple(np,NEp,Euse,Nmh_,Nmh);
  SumPointsSimple(np,NEp,Euse,Nbb_,Nbb); 
  SumPointsSimple(np,NEp,Euse,Le_,Le); 
  SumPointsSimple(np,NEp,Euse,Lp_,Lp);     
  NumEpoints=NEp;
  int numpar=4;
  if(FreeEnergyFit==1) numpar+=NEp;
  Double_t ECorrBB=0;
  Double_t LG=0,Lee=0;
  for(int is=0;is<NEp;is++)
    {
      EInScan[is]=En[is];
      WInScan[is]=2.*En[is];
      EErrInScan[is]=Eerr[is];
      WErrInScan[is]=Eerr[is]*2.;
      NmhInScan[is]=Nmh[is];   
      NbbInScan[is]=Nbb[is];       
      ECorrBB=1./CrossSBhabhaPP(En[is],&arguments.CrBhabha);
      LumInScan[is]=NbbInScan[is]*ECorrBB;
      //if(arguments.NeeFlag==1) LumInScan[is]*=9.54; //this is for gamma gamma luminocity
      LG+=TMath::Max(Le[is],Lp[is]);
      Lee+= LumInScan[is];
      cout<<"LumG:"<<Le[is]<<" Lp:"<<Lp[is]<<" LumInScan[is]:"<<LumInScan[is]<<endl;
      LumLgammaInScan[is]=TMath::Max(Le[is],Lp[is]);              
        
        
      if(UseLumBB==0){          
	CrossSInScan[is]=        Nmh[is]/LumLgammaInScan[is];
	if(Nmh[is]>4){
	  CrossSErrInScan[is]=sqrt(Nmh[is])/LumLgammaInScan[is];        
	}
	else
	  {
	    CrossSErrInScan[is]=sqrt(Nmh[is]+4)/LumLgammaInScan[is];     
	  }
          
        CrossSBBInScan[is]=Nbb[is]/LumLgammaInScan[is];
        if(Nbb[is]>4){
          CrossSBBErrInScan[is]=sqrt(Nbb[is])/LumLgammaInScan[is];        
        }

        else
          {
            CrossSBBErrInScan[is]=sqrt(Nbb[is]+4)/LumLgammaInScan[is];         
          }
      }
      else 
        {
          CrossSInScan[is]=         Nmh[is]/LumInScan[is];
          if(Nmh[is]>4){
            CrossSErrInScan[is]=sqrt(Nmh[is])/LumInScan[is];        
          }
          else
            {
              CrossSErrInScan[is]=sqrt(Nmh[is]+4)/LumInScan[is];     
            }
       
          CrossSBBInScan[is]=Nbb[is]/LumInScan[is];
          if(Nbb[is]>4){
            CrossSBBErrInScan[is]=sqrt(Nbb[is])/LumInScan[is];        
          }
          else
            {
              CrossSBBErrInScan[is]=sqrt(Nbb[is]+4)/LumInScan[is];         
            }

	}
        
    }
  cout<<"LG:"<<LG<<" Lee:"<<Lee<<endl;
  if(arguments.verbose){

    for(int is=0;is<NEp;is++)
      {
	cout<<"point:"<<is<<" Energy:"<<En[is]<<"NMH:"<<Nmh[is]<<endl;	    
          
      }
  }    
  //    GrRes=new TGraphErrors(NEp,WInScan,CrossSInScan,WErrInScan,CrossSErrInScan);
  //   GrRes->Fit("pol0","Q0");
  ///    double par0=GrRes->GetFunction("pol0")->GetParameter(0);
  //    cout<<"par0:"<<par0<<endl;
  

  
  MinuitRes = new TMinuit(numpar);
  if(arguments.Chi2==0){
    MinuitRes->SetFCN(fcnResMult);  

  }
  else {
    // !!     MinuitRes->SetFCN(fcnResChi2);  
    //!!      cout<<"USING CHI2 !!!!"<<endl;
  }
    
  Double_t arglistRes[numpar*2];
    
  Int_t ierflgRes = 0;
  if(arguments.verbose>0)
    {
      arglistRes[0]=-1;
      MinuitRes->mnexcm("SET PRINT", arglistRes,1,ierflgRes);
      arglistRes[0] =0;
      MinuitRes->mnexcm("SET NOW", arglistRes,0,ierflgRes);
    }
  else
    {
      //   arglistRes[0]=0;
      arglistRes[0]=-1;
      MinuitRes->mnexcm("SET PRINT", arglistRes,1,ierflgRes);
      arglistRes[0] = 0;
      MinuitRes->mnexcm("SET NOW", arglistRes,0,ierflgRes);
    }
  MinuitRes->mnexcm("SET ERR", arglistRes,1,ierflgRes);
   
  Double_t vstartRes[5]= {50,0.8,0,1.3,arguments.CrBhabha};   
      
  Double_t stepRes[5] =  {0.10,0.01,0.01,0.01,0.0};
        
     
     
  MinuitRes->SetMaxIterations(20000);                         
  MinuitRes->DefineParameter(0,"bg",vstartRes[0],stepRes[0],-150,150.0);
  MinuitRes->DefineParameter(1,"eff",vstartRes[1],stepRes[1],0.1,1.0);      
  MinuitRes->DefineParameter(2,"dM/2.",vstartRes[2],stepRes[2],-1.,1);      
  MinuitRes->DefineParameter(3,"SigmaW",vstartRes[3],stepRes[3],0.5,1.78);        
  if(FreeEnergyFit==1){
    for(int j=0;j<NEp;j++){
      char  NameP[10];
      sprintf(NameP,"dE%d",j);         
      MinuitRes->DefineParameter(j+4,NameP,0,0.01,-0.3,0.3);        
    }
  }
   
      
  MinuitRes->mnexcm("MIGRAD", arglistRes,numpar,ierflgRes);
#if _HESSE_ 
  MinuitRes->mnexcm("HESSE", arglistRes,0,ierflgRes);
#endif
#if _MINOs_
  MinuitRes->mnexcm("MINOs 1000 3 3", arglistRes,0,ierflgRes);
#endif
  // Print results
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  Double_t * parRes= new Double_t [numpar] ;
  Double_t*       parErrRes= new Double_t [numpar] ;    
  MinuitRes->mnstat(aminRes,edmRes,errdefRes,nvparRes,nparxRes,icstatRes);
  MinuitRes->mnprin(numpar,aminRes);
  for(Int_t i=0;i<numpar;i++)
    {
      MinuitRes->GetParameter(i,parRes[i],parErrRes[i]);
    }
      
  int nf=MinuitRes->GetNumFreePars();
  if(FreeEnergyFit==1) nf-=NEp;
 
  cout<<"Minuit Mass= "<<3686.111+parRes[2]*2.<<endl;
  cout<<"parRes[2]*2.:"<<parRes[2]*2.<<endl;
  Double_t* parPsiPF    = new Double_t [idRNP];
  parPsiPF[idRbg]=parRes[0];
  parPsiPF[idRM]=parRes[2];//parPsiP[Iscan][ippeff];   
  parPsiPF[idRSw]=parRes[3];   
  parPsiPF[idReff]=parRes[1];
  parPsiPF[idRFreeGee]=0;
  parPsiPF[idRTauEff]=0;
     
  for(int is=0;is<NEp;is++)
    {       
      if(arguments.FreeEnergy==1){ 
	WInScan[is]=2.*(En[is]);//+parRes[is+4]);
	WErrInScan[is]=WErrInScan[is];
      }
      else {
	WInScan[is]=2.*En[is];
	WErrInScan[is]=WErrInScan[is];
      }
      
       
    }
  GrRes=new TGraphErrors(NEp,WInScan,CrossSInScan,WErrInScan,CrossSErrInScan);
  TF1* FitPsiP=new TF1("FitPsiP",FCrSPPrimeAzimov,1836.*ScaleEGr,1855*ScaleEGr,idRNP);  
  FitPsiP->SetParameters( parPsiPF);
  TCanvas* TestCanv=new TCanvas("TestCanv","TestCanv",900,700); 
  TestCanv->cd();
  GrRes->Draw("AP");     
  FitPsiP->Draw("SAME");
  char Info1[100];
  TLatex*  latexM1=new TLatex();
  latexM1->SetTextSize(0.038);
  latexM1->SetTextColor(2);       
  sprintf(Info1,"#\chi^{2}_{#\psi(2S)}=%3.3f/ (%d -%d) =%3.3f",MinChi2,NpPP,nf,MinChi2/(NpPP-nf)); 
  double xx=1837*ScaleEGr;
  latexM1->DrawLatex(xx,100,Info1);
  sprintf(Info1,"#\delta M_{#\psi(2S)}=%3.3f#\pm%3.3f [MeV]",_MPsiPrime+parRes[2]*2.-3686.111,parErrRes[2]*2.);
  latexM1->DrawLatex(xx,80,Info1);
  TestCanv->Update();
  delete [] En_;   
  delete [] Eerr_;
  delete [] Nmh_;
  delete [] Nbb_;                   
  delete [] En;   
  delete [] Eerr;
  delete [] Nmh;
  delete [] Nbb;               
  delete FitRes;
   
 

    
 
  theApp->Run();

  
  return 0;
}

static int FCNcall=0;

int FCallCheck=0;

Double_t parprev[10]={0,0,0, 0,0,0, 0,0,0, 0};
void fcnResMult(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    //calculate chisquare
    	Double_t chisq =0;
	Double_t chisqbb = 0;
	Double_t chisqmh =0 ;
	Double_t sigmaFull;
	Double_t sigmaMH;
	Double_t sigmaBB;
	Double_t nFull;
	Double_t lumFull;
	Double_t Energy;
        Double_t parmh[idRNP];
        Double_t EnergyChi2=0;
     	Double_t parbb[1];        
        parmh[idRbg]=par[0];
        parmh[idReff]=par[1];
        parmh[idRM]=par[2];
        parmh[idRSw]=par[3];   
        parmh[idRFreeGee]=0;       
        parmh[idRTauEff]=0;       
        parbb[0]=CrossBhabha;        
          for (Int_t i=0;i<NumEpoints;i++)
          {
            Energy=EInScan[i];
            if(FreeEnergyFit==1){ 
              Energy+=par[4+i];
              EnergyChi2+=(par[4+i]*par[4+i]/(EErrInScan[i]*EErrInScan[i]));
            }
            sigmaMH=CrSOniumR(_MethodAzimov,_IdPsiPrime,Energy,parmh);  
            sigmaBB=CrossSBhabhaPP(Energy,parbb);                        
            sigmaFull=sigmaMH+sigmaBB;
            nFull=NbbInScan[i]+NmhInScan[i];                                
            lumFull=nFull/sigmaFull;
            if(UseLumBB==0)  {
              lumFull=LumLgammaInScan[i];
              nFull=NmhInScan[i];
              sigmaFull=sigmaMH;
            }
	     
          if( FCNcall ==0 ){
            cout<<"Point "<< i << " E " << Energy
                << " Nmh " << NmhInScan[i]<<"Nbb:"<<NbbInScan[i]<<
             " sigmaBB:"<<sigmaBB<<
             " sigmaMH:"<<sigmaMH<<" rmh:"<<NmhInScan[i]/ sigmaMH<<
              " rbb: "<<NbbInScan[i] / sigmaBB<<" LumFull:"<<lumFull<<
              " LG:"<<LumLgammaInScan[i]<<
              " parmh0:"<<parmh[0]<<"parmh[1]:"<<parmh[1]<<
              "parmh[2]:"<<parmh[2]<<"parmh[3]:"<<parmh[3]
                <<endl;
               NpPP++;
            }
          chisqmh=0;
           if(NmhInScan[i]>0)
             {
               chisqmh= (NmhInScan[i]*log(NmhInScan[i]/(sigmaMH*lumFull))
                         +sigmaMH*lumFull-NmhInScan[i]);            
             }
           else if(NmhInScan[i]==0)  
           {
             chisqmh =  sigmaMH*lumFull;
           }
            chisq+=  (2*chisqmh);
                     
            
            if(UseLumBB!=0){
              if(NbbInScan[i]>0)
                {
                  chisqbb= (NbbInScan[i]*log(NbbInScan[i]/(sigmaBB*lumFull))+
                            sigmaBB*lumFull-NbbInScan[i]);
               
                }
              else if(NbbInScan[i]==0)
              {
                 chisqbb =  sigmaBB*lumFull;
                
              }
              chisq+=  (2*chisqbb);
              
              
            }
              /*            chisqbb= sq(sigmaBB*LumLgammaInScan[i]-NbbInScan[i])/sq(NbbInScan[i]+1);
            chisqbb= sq(sigmaBB-CrossSBBInScan[i])/sq(CrossSBBErrInScan[i]);
           
            chisq+=  chisqbb;        

            }*/               
          }

          
          f = chisq;
          if(FreeEnergyFit==1) f+=EnergyChi2;//;    
          if(MinChi2>f) MinChi2=f;
          FCNcall=1;
     
}
