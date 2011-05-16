#define _HESSE_ 1 // calculate errors with HESSE
#define _MINOs_ 1 // use minos errors

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
#include <TLegend.h>
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
Double_t SigmaWInScan[NumMaxP];
Double_t dSigmaWInScan[NumMaxP];
Double_t WInScan[NumMaxP];
Double_t EErrInScan[NumMaxP];
Double_t WErrInScan[NumMaxP];
Int_t    NumEpoints=0;
Double_t MinChi2=1e+7;
void fcnResMult    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult2    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult_both2    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcnResMult3    (Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

/* Program documentation. */
static char doc[] ="Searching narrow resonance";

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
#define OPT_both           13
#define OPT_SW             14
#define OPT_USE_CHI2             15


bool USE_CBS_SIGMAW = false;
bool USE_CHI2 = false;
Int_t   FreeEnergyFit=0;

enum LuminosityType
{
  BESLUM, NEELUM, NGGLUM
};

LuminosityType LUMINOSITY;


double CHI2_TOTAL;
double CHI2_LUM;
double CHI2_ENERGY;
double CHI2_SIGNAL;


inline double sq(double x) { return x*x; }

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
  {"both",OPT_both,"both",0,"both",13},     
  {"SW",OPT_SW,"SW",0,"SW",14},     
  {"use-chi2",OPT_USE_CHI2,"use-chi2",0,"use-chi2",15},     
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
  int both;
  int scan;
  int SW; //use cbs measure sigmaW
  int use_chi2; //use chi2 algorithm
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
    case OPT_both:arg_union.arguments->both=atoi(arg);break;            
    case OPT_SW:arg_union.arguments->SW=atoi(arg);break;            
    case OPT_USE_CHI2:arg_union.arguments->use_chi2=atoi(arg);break;            
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
  arguments.SW=0;
  arguments.use_chi2=0;
  //arguments.CrBhabha=62;  
  // arguments.CrBhabha=535;  

  arguments.CrBhabha=80; 
  arguments.both=0; 



  arguments.dEmin=0.02 ;  
  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  TApplication* theApp=0;

  cout<<"JpsiFit:"<<arguments.JpsiFit<<"Ch2:"<<arguments.Chi2<<endl;
  JpsiFitOnly=arguments.JpsiFit;
  UseLumBB=arguments.lumbb;
  double Kee;
  double Kgg;

  if(arguments.SW==1) USE_CBS_SIGMAW=true;
  if(arguments.use_chi2==1) USE_CHI2=true;
  ifstream kfile("CrBhabha.txt");
  if(kfile) 
  {
    kfile >> Kee;
    kfile >> Kgg;
  }
  kfile.close();
  //arguments.CrBhabha=Kee;
  arguments.CrBhabha=0;
  //Fast luminosity chosen
  switch(arguments.lumbb)
  {
    case 1:
      arguments.NeeFlag=1;
      break;
    case 2:
      arguments.NeeFlag=2;
      break;
    default:
      break;
  }

  switch(arguments.NeeFlag)
  {
    case 0:
      LUMINOSITY = BESLUM;
      arguments.CrBhabha=Kee;
      break;
    case 1:
      LUMINOSITY = NEELUM;
      arguments.CrBhabha=Kee;
      break;
    case 2:
      LUMINOSITY = NGGLUM;
      arguments.CrBhabha=Kgg;
      break;
    default:
      LUMINOSITY = BESLUM;
      arguments.CrBhabha=Kee;
      break;
  }
  CrossBhabha= arguments.CrBhabha;
  cout << "Cross section of luminosity measurement process: " << CrossBhabha << " nb" << endl;
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
  int dimMHFile=9;
  int MHRun=0;
  int MHLum=1;
  int MHEnergy=2;
  int MHEnergyErr=3;
  int MHSigmaW=4;
  int MHdSigmaW=5;
  int MHNh=6; 
  int MHNee=7;
  int MHNgg=8; 
  double** AllMH=0;
  int      npMHFile=0; 

  /* ***** Reading data from file ****************** */
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

  int dimAP= 17;
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
  int ASigmaW = 15;
  int AdSigmaW = 16;
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
    AP[Aind][ALe]=AllMH[i][MHLum];//Fill luminocity from bes online lum
    AP[Aind][ALp]=AllMH[i][MHLum];//This is the same for both electron and positron
    AP[Aind][AMHEv]=AllMH[i][MHNh];
    AP[Aind][ASigmaW] = AllMH[i][MHSigmaW];
    AP[Aind][AdSigmaW] = AllMH[i][MHdSigmaW];
    switch(arguments.NeeFlag)
    {
      //Nee luminosity
      case 0:
      case 1:
        AP[Aind][AEE]=AllMH[i][MHNee];
        break;
      //Ngg luminosity
      case 2:
        AP[Aind][AEE]=AllMH[i][MHNgg];
        break;
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

  /*  SKIP SORTING OF POINTS IF separate scan*/
  if(arguments.both==0)
  {
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
  int np=npAP ;   
  cout<<"EMin:"<<EMin<<"EMax:"<<EMax<<"np:"<<np<<endl;

  Double_t *En_=new Double_t[np];
  Double_t *Eerr_=new Double_t[np];
  Double_t *SigmaW_=new Double_t[np];
  Double_t *dSigmaW_=new Double_t[np];
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
    Nbb_[np]=AP[i][AEE];  //here will be Nee or Ngg
    Le_[np]=AP[i][ALe];   //here will be BES online lum
    Lp_[np]=AP[i][ALp];  
    SigmaW_[np]=AP[i][ASigmaW];
    dSigmaW_[np]=AP[i][AdSigmaW];
    np++;      
  }
  int  NEp=0;
  SeparatePointsPartNew(np,Nbb_,&NEp,Euse,En_,arguments.dEmin);      
  cout<<"!!!dEmin:"<<arguments.dEmin<<endl;
  Double_t *En  =new Double_t[NEp];
  Double_t *Eerr=new Double_t[NEp];
  Double_t *SW  =new Double_t[NEp];
  Double_t *dSW =new Double_t[NEp];
  Double_t *Nmh=new  Double_t[NEp];
  Double_t *Nbb=new  Double_t[NEp];
  Double_t *Le=new  Double_t[NEp];
  Double_t *Lp=new  Double_t[NEp];

  SumPointsByQuant(np,NEp,Euse,En_,Nbb_,Eerr_,En,Eerr,false);
  SumPointsSimple(np,NEp,Euse,Nmh_,Nmh);
  SumPointsSimple(np,NEp,Euse,Nbb_,Nbb); 
  SumPointsSimple(np,NEp,Euse,Le_,Le); 
  SumPointsSimple(np,NEp,Euse,Lp_,Lp);     
  SumPointsSimple(np,NEp,Euse,SigmaW_,SW);     
  SumPointsSimple(np,NEp,Euse,dSigmaW_,dSW);     
  NumEpoints=NEp;
  int numpar=4;
  if(FreeEnergyFit==1) numpar+=NEp;
  if(arguments.both==1) numpar+=2;
  if(arguments.both==2) numpar+=3;
  Double_t ECorrBB=0;
  Double_t LG=0,Lee=0;
  //fill global arrayes which needed to access from the FCN function
  for(int is=0;is<NEp;is++)
  {
    EInScan[is]=En[is];
    WInScan[is]=2.*En[is];
    EErrInScan[is]=Eerr[is];
    WErrInScan[is]=Eerr[is]*2.;
    SigmaWInScan[is]=SW[is];
    dSigmaWInScan[is]=dSW[is];
    NmhInScan[is]=Nmh[is];   
    NbbInScan[is]=Nbb[is];       
    ECorrBB=1./CrossSBhabhaPP(En[is],&arguments.CrBhabha);
    LumInScan[is]=NbbInScan[is]*ECorrBB;
    LG+=TMath::Max(Le[is],Lp[is]);
    Lee+= LumInScan[is];
    LumLgammaInScan[is]=TMath::Max(Le[is],Lp[is]); //BES lum
    cout<<"BESLUM:"<<LumLgammaInScan[is]<<" LumInScan:"<<LumInScan[is]<<endl;
    double lum=0; //temporary luminosity could be ee, gg and bes lum
    if(UseLumBB==0)
    {
      switch(LUMINOSITY)
      {
        case BESLUM:
          lum=LumLgammaInScan[is];
          break;
        case NEELUM:
        case NGGLUM:
          lum=LumInScan[is];
          break;
      }
    }
    else lum=LumInScan[is];
    //calculate mhadr cross section
    CrossSInScan[is]=Nmh[is]/lum;
    if(Nmh[is]>4) CrossSErrInScan[is]=sqrt(Nmh[is]*(1.+Nmh[is]/NbbInScan[is]))/lum;        
    else CrossSErrInScan[is]=sqrt(Nmh[is]+4)/lum;     
    //calculate bhabha or gg cross section
    CrossSBBInScan[is]=Nbb[is]/lum;
    if(Nbb[is]>4) CrossSBBErrInScan[is]=sqrt(Nbb[is])/lum;        
    else CrossSBBErrInScan[is]=sqrt(Nbb[is]+4)/lum;         
  }
  cout<<"Total luminosity: BES lum:"<<LG<<" Lee or Lgg:"<<Lee<<endl;
  if(arguments.verbose)
  {
    for(int is=0;is<NEp;is++)
    {
      cout<<"point:"<<is<<" Energy:"<<En[is]<<"NMH:"<<Nmh[is]<<endl;	    
    }
  }    

  MinuitRes = new TMinuit(numpar);
  if(arguments.Chi2==0)
  {
    MinuitRes->SetFCN(fcnResMult);  
    if(arguments.both==1 && FreeEnergyFit==1) MinuitRes->SetFCN(fcnResMult2);  
    if(arguments.both==1 && FreeEnergyFit==0) MinuitRes->SetFCN(fcnResMult3);  
    if(arguments.both==2 && FreeEnergyFit==1) MinuitRes->SetFCN(fcnResMult_both2);  
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
  arglistRes[0] = 2;
  MinuitRes->mnexcm("SET STRATEGY", arglistRes,1,ierflgRes);

  Double_t vstartRes[5]= {0.8,0.01,0.0,1.59,arguments.CrBhabha};   

  Double_t stepRes[5] =  {0.1,0.01,0.01,0.01,0.0};



  MinuitRes->SetMaxIterations(100000000);                         
  //MinuitRes->DefineParameter(0,"bg",vstartRes[0],stepRes[0],-150,150.0);
  //MinuitRes->DefineParameter(1,"eff",vstartRes[1],stepRes[1],0.01,1.0);      
  //MinuitRes->DefineParameter(2,"dM/2.",vstartRes[2],stepRes[2],-0.5,0.5);      
  //MinuitRes->DefineParameter(3,"SigmaW",vstartRes[3],stepRes[3],0.5,1.8);
  MinuitRes->DefineParameter(0,"bg",vstartRes[0],stepRes[0],-100,+100);
  MinuitRes->DefineParameter(1,"eff",vstartRes[1],stepRes[1],0.0,1.0);      
  MinuitRes->DefineParameter(2,"dM/2.",vstartRes[2],stepRes[2],-0.5,0.5);      
  MinuitRes->DefineParameter(3,"SigmaW",vstartRes[3],stepRes[3],0.5,2.0);
  //if(USE_CBS_SIGMAW) MinuitRes->FixParameter(3);
  if(FreeEnergyFit==1)
  {
    for(int j=0;j<NEp;j++)
    {
      char  NameP[10];
      sprintf(NameP,"dE%d",j);         
      MinuitRes->DefineParameter(j+4,NameP,0,0.01,-0.5,+0.5);        
    }
  }

  int nep = FreeEnergyFit==0 ? 0 : NEp;
  if(arguments.scan==3)
  {
    cout << "Fit scan1 and scan2"<< endl;
    if(arguments.both==1)
    {
      MinuitRes->DefineParameter(0+4+nep,"bg2",vstartRes[0],stepRes[0],-150,150.0);
      MinuitRes->DefineParameter(1+4+nep,"eff2",vstartRes[1],stepRes[1],0.1,1.0);      
      cout << "Variant with own bg and eff" << endl;
    }
    if(arguments.both==2)
    {
      MinuitRes->DefineParameter(0+4+nep,"bg2",vstartRes[0],stepRes[0],-150,150.0);
      MinuitRes->DefineParameter(1+4+nep,"eff2",vstartRes[1],stepRes[1],0.1,1.0);      
      MinuitRes->DefineParameter(2+4+nep,"SigmaW2",vstartRes[3],stepRes[3],0.5,1.8);      
      cout << "Variant with own bg, eff and sigmaW" << endl;
    }

  }

  MinuitRes->mnexcm("MIGRAD", arglistRes,numpar,ierflgRes);
  MinuitRes->mnimpr();
  MinuitRes->mnexcm("HESSE", arglistRes,0,ierflgRes);
  MinuitRes->mnexcm("MINOs 10000000 3 3", arglistRes,0,ierflgRes);
  // Print results
  Double_t aminRes,edmRes,errdefRes;
  Int_t nvparRes,nparxRes,icstatRes;
  Double_t * parRes= new Double_t [numpar] ;
  Double_t*       parErrRes= new Double_t [numpar] ;    
  for(Int_t i=0;i<numpar;i++)
  {
    MinuitRes->GetParameter(i,parRes[i],parErrRes[i]);
    //cout << i << " " << parRes[i] << endl;
  }
  Int_t npar=numpar;
  Double_t grad=0;
  Double_t fval=0;
  Int_t flag=1;
  MinuitRes->Eval(npar, &grad, fval, parRes, flag);
  cout << "fval = " << fval << endl;

  MinuitRes->mnstat(aminRes,edmRes,errdefRes,nvparRes,nparxRes,icstatRes);
  MinuitRes->mnprin(numpar,aminRes);

  int nf=MinuitRes->GetNumFreePars();
  if(FreeEnergyFit==1) nf-=NEp;

  cout.precision(15);
  cout<<"Minuit Mass= "<<_MPsiPrime+parRes[2]*2.<<endl;
  cout<<"PDG Mass= "<<_MPsiPrime<<endl;
  cout<<"parRes[2]*2.:"<<parRes[2]*2.<< " +- " << parErrRes[2]*2. <<  " MeV. chi2/ndf = " <<MinChi2 << "/(" << NpPP<<"-"<<nf<<") = "  << MinChi2/(NpPP-nf) <<endl;
  cout << "Contribution to chi square:" << endl;
  cout.precision(4);
  cout << "chi2 signal: " << CHI2_SIGNAL <<  " or " << CHI2_SIGNAL/CHI2_TOTAL*100 << "%" << endl;
  cout << "chi2 energy: " << CHI2_ENERGY <<  " or " << CHI2_ENERGY/CHI2_TOTAL*100 << "%" <<  endl;
  cout << "chi2 lum: "    << CHI2_LUM    <<  " or " << CHI2_LUM/CHI2_TOTAL*100 << "%" << endl;
  cout << "Total chi2: " << CHI2_TOTAL << endl;
  Double_t* parPsiPF    = new Double_t [idRNP];
  Double_t* parPsiPF2    = new Double_t [idRNP];
  parPsiPF[idRbg]=parRes[0];
  parPsiPF[idRM]=parRes[2];//parPsiP[Iscan][ippeff];   
  parPsiPF[idRSw]=parRes[3];   
  parPsiPF[idReff]=parRes[1];
  parPsiPF[idRFreeGee]=0;
  parPsiPF[idRTauEff]=0;


  if(arguments.scan==3)
  {
    int nep = arguments.FreeEnergy == 0 ? 0 : NEp;
    if(arguments.both==1)
    {
      parPsiPF2[idRbg]=parRes[0+4+nep];
      parPsiPF2[idReff]=parRes[1+4+nep];
      parPsiPF2[idRM]=parRes[2];
      parPsiPF2[idRSw]=parRes[3];   
      parPsiPF2[idRFreeGee]=0;
      parPsiPF2[idRTauEff]=0;
      for(unsigned i=0; i<idRNP;i++)
      {
        cout << i << " " << parPsiPF2[i] << endl;
      }
    }
    if(arguments.both==2)
    {
      parPsiPF2[idRbg]=parRes[0+4+nep];
      parPsiPF2[idReff]=parRes[1+4+nep];
      parPsiPF2[idRM]=parRes[2];
      parPsiPF2[idRSw]=parRes[2+4+nep];   
      parPsiPF2[idRFreeGee]=0;
      parPsiPF2[idRTauEff]=0;
      for(unsigned i=0; i<idRNP;i++)
      {
        cout << i << " " << parPsiPF2[i] << endl;
      }
    }
  }


  double EnergyChi2=0;
  for(int is=0;is<NEp;is++)
  {       
    //WInScan[is]=2.*En[is];
    WErrInScan[is]=WErrInScan[is];
    //if(arguments.FreeEnergy==1) WInScan[is]+=2.*+parRes[is+4];
    EnergyChi2+= sq(parRes[is+4]/EErrInScan[is]);
  }

  cout << "Energy Chi2 contribution:" << EnergyChi2 << endl;

  GrRes=new TGraphErrors(NEp,WInScan,CrossSInScan,WErrInScan,CrossSErrInScan);
  TF1* FitPsiP=new TF1("FitPsiP",FCrSPPrimeAzimov,1836.*ScaleEGr,1855*ScaleEGr,idRNP);  
  TF1* FitPsiP2=new TF1("FitPsiP",FCrSPPrimeAzimov,1836.*ScaleEGr,1855*ScaleEGr,idRNP);  
  FitPsiP->SetParameters( parPsiPF);
  FitPsiP2->SetParameters( parPsiPF2);
  TCanvas* TestCanv=new TCanvas("BES_PSIP_SCAN","BES PsiP Scan",900,700); 
  TestCanv->SetGridx();
  TestCanv->SetGridy();
  TestCanv->SetFillColor(0);
  TestCanv->SetBorderMode(0);
  //TestCanv->cd();
  GrRes->Draw("AP");
  gPad->SetBorderMode(0);
  GrRes->GetXaxis()->SetTitle("W [MeV]");
  GrRes->GetYaxis()->SetTitle("#sigma [nb]");
  GrRes->SetTitle("#psi' scan");
  FitPsiP->Draw("SAME");
  if(arguments.both==1 || arguments.both==2)
  {
    FitPsiP2->SetLineColor(kBlue);
    FitPsiP2->Draw("SAME");
    TLegend * l= new TLegend(0.8,0.9,1.0,1.0);
    l->AddEntry(FitPsiP,"First scan", "l");
    l->AddEntry(FitPsiP2,"Second scan", "l");
    l->Draw();
  }
  char Info1[100];
  TLatex*  latexM1=new TLatex();
  latexM1->SetTextSize(0.038);
  latexM1->SetTextColor(2);       
  sprintf(Info1,"#chi^{2}_{#psi(2S)}=%3.3f/ (%d -%d) =%3.3f",MinChi2,NpPP,nf,MinChi2/(NpPP-nf)); 
  double xx=1838*ScaleEGr;
  //double yy = parRes[0]*5;
  double yy = FitPsiP->GetMaximum();
  latexM1->DrawLatex(xx,yy,Info1);
  sprintf(Info1,"#Delta M_{#psi(2S)}=%3.3f#pm%3.3f [MeV]",parRes[2]*2.,parErrRes[2]*2.);
  latexM1->DrawLatex(xx,yy*0.8,Info1);
  TLatex * latexSw = new TLatex();
  latexSw->SetTextSize(0.038);
  latexSw->SetTextColor(2);
  sprintf(Info1,"#sigma_{W}=%1.3f #pm %1.3f [MeV]",parRes[3],parErrRes[3]); 
  latexSw->DrawLatex(xx,yy*0.6,Info1);

  TLatex * latexProb = new TLatex();
  latexProb->SetTextSize(0.038);
  latexProb->SetTextColor(2);
  sprintf(Info1,"P(#chi^{2})=%1.4f",TMath::Prob(MinChi2,NpPP-nf)); 
  latexProb->DrawLatex(xx,yy*0.4,Info1);

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
  Double_t SigmaWChi2=0;
  Double_t LumChi2=0;
  Double_t SignalChi2=0;
  Double_t parbb[1];        
  parmh[idRbg]=par[0];
  parmh[idReff]=par[1];
  parmh[idRM]=par[2];
  parmh[idRSw]=par[3];   
  parmh[idRFreeGee]=0;       
  parmh[idRTauEff]=0;       
  parbb[0]=CrossBhabha;        
  if(iflag==1)
  {
    cout << setw(4) << "Parnumber" <<  setw(20) << "value" << endl;
    for(int i=0;i<npar;i++)
    {
      cout << setw(4)<< i << setw(20) << par[i] << endl;
    }
    cout << setw(4)<< "point#" << setw(15)<< "energy, MeV" <<  setw(15) << "energy chi2" << setw(15) <<  "signal chi2" << endl;
  }
  for (Int_t i=0;i<NumEpoints;i++)
  {
    Energy=EInScan[i];
    double echi2=0;
    if(FreeEnergyFit==1)
    { 
      Energy+=par[4+i];
      echi2=(par[4+i]*par[4+i]/(EErrInScan[i]*EErrInScan[i]));
      EnergyChi2+=echi2;
    }
    if(USE_CBS_SIGMAW)
    {
      //parmh[idRSw]=SigmaWInScan[i];
      SigmaWChi2+= sq((parmh[idRSw] - SigmaWInScan[i])/dSigmaWInScan[i]);
    }
    sigmaMH=CrSOniumR(_MethodAzimov,_IdPsiPrime,Energy,parmh);  
    sigmaBB=CrossSBhabhaPP(Energy,parbb);                        
    sigmaFull=sigmaMH+sigmaBB;
    nFull=NbbInScan[i]+NmhInScan[i]; 
    lumFull=nFull/sigmaFull;
    if(UseLumBB==0)
    {
      switch(LUMINOSITY)
      {
        case BESLUM:
          lumFull=LumLgammaInScan[i];
          break;
        case NEELUM:
          lumFull=LumInScan[i];
          break;
        default:
          lumFull=LumInScan[i];
          break;
      }
      nFull=NmhInScan[i];
      sigmaFull=sigmaMH;
    }

    if( FCNcall ==0 )
    {
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
    //calculate signal contribution
    if(NmhInScan[i]>0)
    {
      double N=sigmaMH*lumFull; //expected number of events
      if(USE_CHI2)
      {
        chisqmh = sq(N - NmhInScan[i])/(NmhInScan[i]*(1. + NmhInScan[i]/(NbbInScan[i]/10.))); //just chi square
      }      
      else chisqmh = 2*(NmhInScan[i]*log(NmhInScan[i]/N) +N - NmhInScan[i]); //likelihood for low statistics
    }
    else if(NmhInScan[i]==0) chisqmh =  sigmaMH*lumFull;
    SignalChi2+=chisqmh;
    chisq+=chisqmh; 
    if(iflag==1)
    {
      cout << setw(4)<< i << setw(15)<< EInScan[i] << setw(15) << echi2 << setw(15) <<  chisqmh << endl;
    }


    if(UseLumBB!=0)
    {
      if(NbbInScan[i]>0)
      {
        double N =sigmaBB*lumFull;
        if(USE_CHI2) chisqbb = sq(N - NbbInScan[i])/NbbInScan[i]; //just chi square
        else chisqbb = 2*(NbbInScan[i]*log(NbbInScan[i]/N)+ N - NbbInScan[i]); //likelihood for low statistics
      }
      else if(NbbInScan[i]==0)
      {
        chisqbb =  sigmaBB*lumFull;
      }
      chisq+=chisqbb;
      LumChi2 += chisqbb;
    }
  }
  f = chisq;
  if(FreeEnergyFit==1) f+=EnergyChi2;//;    
  if(USE_CBS_SIGMAW) f+=SigmaWChi2;    
  if(MinChi2>f)
  {
    MinChi2=f;
    CHI2_ENERGY = EnergyChi2;
    CHI2_TOTAL=f;
    CHI2_LUM = LumChi2;
    CHI2_SIGNAL = SignalChi2;
  }
  FCNcall=1;
}

void fcnResMult2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
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
  int Nscan1=6; //number of point in first scan
  for (Int_t i=0;i<Nscan1;i++)
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
    if(UseLumBB==0)
    {
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
  }
  parmh[idRbg]=par[4+NumEpoints];
  parmh[idReff]=par[5+NumEpoints];

  for (Int_t i=Nscan1;i<NumEpoints;i++)
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
    if(UseLumBB==0)
    {
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
      //if(USE_GAUS)
      //{
      //  chisqmh = 
      //}
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
  }
  f = chisq;
  if(FreeEnergyFit==1) f+=EnergyChi2;//;    
  if(MinChi2>f) MinChi2=f;
  FCNcall=1;
}

void fcnResMult3(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
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
  for (Int_t i=0;i<6;i++)
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
  }
  parmh[idRbg]=par[4];
  parmh[idReff]=par[5];

  for (Int_t i=6;i<NumEpoints;i++)
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
  }
  f = chisq;
  if(FreeEnergyFit==1) f+=EnergyChi2;//;    
  if(MinChi2>f) MinChi2=f;
  FCNcall=1;
}

void fcnResMult_both2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
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
  int Nscan1=6; //number of point in first scan
  for (Int_t i=0;i<Nscan1;i++)
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
  }
  parmh[idRbg]=par[0+4+NumEpoints];
  parmh[idReff]=par[1+4+NumEpoints];
  parmh[idRSw]=par[2+4+NumEpoints];

  for (Int_t i=Nscan1;i<NumEpoints;i++)
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
  }
  f = chisq;
  if(FreeEnergyFit==1) f+=EnergyChi2;//;    
  if(MinChi2>f) MinChi2=f;
  FCNcall=1;
}

