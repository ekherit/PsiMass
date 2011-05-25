/*
 * =====================================================================================
 *
 *       Filename:  first-scan.cpp
 *
 *    Description:  Make result for first scan of psi prime at BES-3
 *
 *        Version:  1.0
 *        Created:  06.12.2010 18:47:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TCut.h>

#include <iomanip>
#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "averager.h"
using namespace std;



const double BEPC_ALPHA=0.022; //BEPC crossing angle
const double Me=0.510998918; //Electron mass, MeV

#include "first-scan.h"
TTree * get_tree(const char * file)
{
	TFile * mc_file = new TFile(file);
  if(mc_file->IsOpen())
  {
    TTree * mc = (TTree*)mc_file->Get("mhadr");
    //mc->AddFriend("dedx");
    return mc;
  }
  else
    return 0;
}

inline double sq(double x) {return  x*x; }
struct RunInfo_t
{
  unsigned run;
  unsigned NBee;
  unsigned NEee;
  unsigned Nuu;
  unsigned Ngg;
  unsigned Nh;
  double SNR;
  double lum;
  double xsec;
  RunInfo_t(void)
  {
    run=0;
    NBee=0;
    NEee=0;
    Nuu=0;
    Ngg=0;
    Nh=0;
    SNR=0;
    lum=0;
    xsec=0;
  }
  RunInfo_t(unsigned r, unsigned nbee, unsigned beee, unsigned nuu, unsigned ngg, unsigned nh, double snr, double l, double xs)
  {
    run=r;
    NBee=nbee;
    NEee=beee;
    Nuu=nuu;
    Ngg=ngg;
    Nh=nh;
    SNR=snr;
    lum=l;
    xsec=xs;
  }
};

struct ScanPoint_t
{
  list <unsigned> runs; //list of runs
  list <RunInfo_t> ri;
  unsigned scan; //scan number
  double lum;
  double W; //cm energy
  double dW; //cm energy error
  double Ee,dEe; //electron energy and error
  double Ep,dEp; //positron energy and error
  double Sw, dSw; //SigmaW, and error
  double Se, dSe; //Sigma for electron
  double Sp, dSp; //Sigma for positron
  unsigned long Nh; //number of multihadron
  unsigned long Nee; //number of e+e- (bhabha)
  unsigned long Ngg;  //number of gamma gamma
  ibn::averager<double> Nchtr, Nntr; //number of charged tracks and it rms
  unsigned pn;
  ScanPoint_t(void)
  {
    //BES part
    lum=0;
    scan=0;
    //Selection part
    Nh=0;
    Nee=0;
    Ngg=0;
    //CBS part
    W=0;
    dW=0;
    pn=0;
    Ee=0;
    dEe=0;
    Ep=0;
    dEp=0;
    Sw=0;
    dSw=0;
    Se=0;
    dSe=0;
    Sp=0;
    dSp=0;
  }

  ScanPoint_t & operator=(const ScanPoint_t & sp)
  {
    runs=sp.runs; //list of runs
    ri=sp.ri;
    lum=sp.lum;
    pn = sp.pn;
    scan=sp.scan;

    //cm energy
    W=sp.W;
    dW=sp.dW;
    Sw=sp.Sw;
    dSw=sp.dSw; //SigmaW, and error

    Ee=sp.Ee;
    dEe=sp.dEe;
    Se=sp.Se;
    dSe=sp.dSe; //Sigma for electron

    Ep=sp.Ep;
    dEp=sp.dEp;
    Sp=sp.Sp;
    dSp=sp.dSp; //Sigma for positron

    Nh=sp.Nh; //number of multihadron
    Nee=sp.Nee; //number of e+e- (bhabha)
    Ngg=sp.Ngg;  //number of gamma gamma
    return *this;
  }
};


ScanPoint_t * get_scan_point_by_run(vector <ScanPoint_t> & pv, unsigned run)
{
  for(unsigned pn=0; pn<pv.size(); pn++)
  {
    for(list<unsigned>::iterator i=pv[pn].runs.begin();i!=pv[pn].runs.end(); ++i)
      if(*i==run) return &pv[pn];
  }
  return 0;
}

void track_number(void)
{
  TChain * chain = new TChain("head","head");
  for(unsigned run=20334; run!=20368;++run)
  {
    char buf[1024];
    sprintf(buf,"psip-%d.root",run);
    cout << run << " test" << endl;
    chain->AddFile(buf);
    //TFile file(buf);
    //if(file.IsOpen())
    //{
    //  TTree * head = (TTree*)file.Get("head");
    //}
  }
  chain->Draw("nchtr:nchtr_rms/sqrt(nsel):run:0","","goff");
  TGraphErrors * g = new TGraphErrors(chain->GetSelectedRows(), chain->GetV3(), chain->GetV1(), chain->GetV4(), chain->GetV2());
  g->Draw("a*");
  chain->Draw("nntr:nntr_rms/sqrt(nsel):run:0","","goff");
  TGraphErrors * gn = new TGraphErrors(chain->GetSelectedRows(), chain->GetV3(), chain->GetV1(), chain->GetV4(), chain->GetV2());
  gn->Draw("a*");
}

void set_alias(TTree * t)
{
  //t->SetAlias("ecut", "mdc.E>=0.05");
  //t->SetAlias("ngt","Sum$(mdc.E>=0.05)");
  //t->SetAlias("ngt_Eemc","Sum$((mdc.E>=0.05)*E)");

  //t->SetAlias("atheta","theta[0]+theta[1]  - 3.1415926535");
  //t->SetAlias("aphi",  "abs(phi[0]-phi[1]) - 3.1415926535");
}


void make_runinfo(list<RunInfo_t> & runinfo)
{
  runinfo.push_back(RunInfo_t());
  runinfo.push_back(RunInfo_t(20334,	8352	,15402	,314	,852	,2219	  ,0.00795939	  ,68.8556	,26.137));
  runinfo.push_back(RunInfo_t(20335,	4348	,8052	  ,147	,457	,1121	  ,0.000277673	,35.7011	,0));
  runinfo.push_back(RunInfo_t(20336,	34133	,61667	,1029	,3697	,8049	  ,0.00847031	  ,277.217	,28.9759));
  runinfo.push_back(RunInfo_t(20337,	850	  ,1488	  ,28	  ,105	,216	  ,0.00666954	  ,6.65356	,20.6225));
  runinfo.push_back(RunInfo_t(20339,	22829	,41971	,781	,2381	,5773	  ,0.00654966	  ,187.524	,43.6648));
  runinfo.push_back(RunInfo_t(20340,	26523	,48033	,736	,2681	,22819	,0.0196156	  ,214.428	,122.52));
  runinfo.push_back(RunInfo_t(20341,	14282	,25865	,410	,1574	,12487	,0.013441	    ,115.654	,116.369));
  runinfo.push_back(RunInfo_t(20342,	15217	,26503	,451	,1487	,56684	,0.11754	    ,115.343	,759.517));
  runinfo.push_back(RunInfo_t(20343,	12164	,21310	,351	,1242	,46841	,0.0963056	  ,92.8011	,720.264));
  runinfo.push_back(RunInfo_t(20344,	20069	,34940	,710	,1997	,75500	,0.0706677	  ,152.242	,820.795));
  runinfo.push_back(RunInfo_t(20346,	24059	,42934	,1112	,2502	,46403	,0.0498852	  ,191.064	,337.807));
  runinfo.push_back(RunInfo_t(20347,	18444	,33420	,860	,1930	,35233	,0.032694	    ,148.45	  ,339.611));
  runinfo.push_back(RunInfo_t(20348,	17771	,33069	,650	,1877	,7028	  ,0.0089732	  ,150.38	  ,38.4824));
  runinfo.push_back(RunInfo_t(20349,	7719	,13829	,277	,820	,3124	  ,0.00758064	  ,62.8171	,46.0948));
  runinfo.push_back(RunInfo_t(20350,	14918	,27594	,599	,1597	,6143	  ,0.00656115	  ,125.361	,61.7188));
  runinfo.push_back(RunInfo_t(20351,	8638	,15893	,312	,915	,2252	  ,0.00437843	  ,71.025	  ,18.5507));
  runinfo.push_back(RunInfo_t(20353,	17976	,32131	,689	,1811	,4737	  ,0.00429332	  ,143.833	,34.1345));
  runinfo.push_back(RunInfo_t(20354,	10800	,19357	,384	,1150	,3202	  ,0.00195352	  ,86.676	  ,8.83348));
  runinfo.push_back(RunInfo_t(20355,	20772	,37246	,631	,2130	,7991	  ,0.00651277	  ,172.162	,42.5681));
  runinfo.push_back(RunInfo_t(20356,	1694	,3086	  ,72	  ,194	,608	  ,0.00605984	  ,13.8865	,35.9984));
  runinfo.push_back(RunInfo_t(20357,	8273	,14988	,270	,884	,3287	  ,0.00543686	  ,67.2076	,57.1118));
  runinfo.push_back(RunInfo_t(20358,	23085	,41467	,627	,2422	,38660	,0.0350919	  ,184.104	,301.457));
  runinfo.push_back(RunInfo_t(20359,	15081	,26522	,368	,1569	,25672	,0.0246531	  ,117.939	,301.608));
  runinfo.push_back(RunInfo_t(20360,	21167	,36828	,693	,2134	,78398	,0.0711638	  ,160.531	,689.09));
  runinfo.push_back(RunInfo_t(20361,	14509	,25591	,514	,1391	,53789	,0.0545199	  ,111.575	,704.419));
  runinfo.push_back(RunInfo_t(20362,	23313	,40911	,836	,2401	,76084	,0.0616514	  ,179.373	,655.239));
  runinfo.push_back(RunInfo_t(20363,	13483	,23712	,614	,1373	,44075	,0.0501659	  ,104.171	,654.185));
  runinfo.push_back(RunInfo_t(20364,	20787	,37980	,779	,2179	,17062	,0.0191568	  ,170.787	,124.873));
  runinfo.push_back(RunInfo_t(20365,	18771	,34205	,757	,2013	,15877	,0.0165697	  ,154.353	,136.972));
  runinfo.push_back(RunInfo_t(20366,	24053	,43467	,907	,2463	,11139	,0.0105485	  ,196.233	,49.6093));
  runinfo.push_back(RunInfo_t(20367,	19150	,35054	,689	,2037	,9176	  ,0.0092355	  ,158.675	,74.4855));
}


/*  calculate c.m.s energy */
double cm_energy(double Ee, double Ep)
{
  double pe = sqrt(Ee*Ee-Me*Me);
  double pp = sqrt(Ep*Ep-Me*Me);
  return sqrt((Ee+Ep)*(Ee+Ep) - pe*pe - pp*pp + 2*pe*pp*cos(BEPC_ALPHA));
  //return 2*sqrt(Ee*Ep)*cos(alpha/2.);
}

double dW_dE1(double E1, double E2)
{
  double W=cm_energy(E1,E2);
  double p1 = sqrt(E1*E1-Me*Me);
  double p2 = sqrt(E2*E2-Me*Me);
  return (E2+p2/p1*E1*cos(BEPC_ALPHA))/W;
}

void print_data(const vector<ScanPoint_t> &pv)
{
  cout << setw(4) << "run" << setw(5) <<"point" 
    << setw(10) << "Ee" << setw(10) << "dEe" << setw(10) << "SigmaWe" << setw(10) << "dSigmaWe" 
    << setw(10) << "Ep" << setw(10) << "dEp" << setw(10) << "SigmaWp" << setw(10) << "dSigmaWp" 
    << setw(10) << "Wcm" << setw(10) << "dWcm"
    << setw(10) << "lum" << setw(16) << "BES run list";
  cout << endl;
  for(unsigned i=0;i<pv.size();i++)
  {
    cout<<setw(4) << i<< setw(5)<<pv[i].pn
      << setw(10) << pv[i].Ee << setw(10)<< pv[i].dEe << setw(10) << pv[i].Se << setw(10)<< pv[i].dSe
      << setw(10) << pv[i].Ep << setw(10)<< pv[i].dEp << setw(10) << pv[i].Sp << setw(10)<< pv[i].dSp
      << setw(10) << pv[i].W  << setw(10) << pv[i].dW
      << setw(10) << pv[i].lum;
    for(list<unsigned>::const_iterator I = pv[i].runs.begin(); I!=pv[i].runs.end(); ++I)
    {
      cout << setw(6) << *I;
    }
    cout << endl;
  }

}

void read_energy(const char * filename,vector <ScanPoint_t> &pv)
{
  ifstream Efile(filename);
  Efile.ignore(1024,'\n');
  cout << "Reading file with cbs measurement: " << filename << endl;
  cout << setw(4) << "run" << setw(5) <<"point" 
    << setw(10) << "Ee" << setw(10) << "dEe" << setw(10) << "SigmaWe" << setw(10) << "dSigmaWe" 
    << setw(10) << "Ep" << setw(10) << "dEp" << setw(10) << "SigmaWp" << setw(10) << "dSigmaWp" 
    << setw(10) << "Wcm" << setw(10) << "dWcm"
    << setw(10) << "lum" << setw(16) << "BES run list";
  cout << endl;
  for(unsigned i=0;i<pv.size();i++)
  {
    int n;
    int pn;
    double Ee, dEe, Se, dSe, Ep, dEp, Sp, dSp;
    Efile >> n >> pn >> Ee>> dEe >> Se >> dSe >> Ep >> dEp >> Sp >> dSp;
    //Spread should be in MeV
    Se*=1e-3;
    Sp*=1e-3;
    dSe*=1e-3;
    dSp*=1e-3;
    cout<<setw(4) << n<< setw(5)<<pn<<setw(10)<<Ee<<setw(10)<<dEe<<setw(10)<<Se
      <<setw(10)<<dSe<<setw(10)<<Ep<<setw(10)<<dEp<<setw(10)<<Sp<<setw(10)<<dSp;
    double Wcm = cm_energy(Ee, Ep);
    double dWcm = sqrt(sq(dW_dE1(Ee,Ep)*dEe) + sq(dW_dE1(Ep,Ee)*dEp));

    cout << setw(10) << Wcm << setw(10) << dWcm;
    cout << setw(10) << pv[i].lum;
    for(list<unsigned>::iterator I = pv[i].runs.begin(); I!=pv[i].runs.end(); ++I)
    {
      cout << setw(6) << *I;
    }

    //fill c.m.energy
    pv[i].W=Wcm;
    pv[i].dW = sqrt(dWcm*dWcm);
    //fill c.m. spread
    pv[i].Sw = sqrt(Sp*Sp+Se*Se);
    pv[i].dSw = sqrt(sq(Sp*dSp)+sq(Se*dSe))/pv[i].Sw;
    //fill electron and positron spread.
    pv[i].Se = Se;
    pv[i].dSe = dSe;
    pv[i].Sp = Sp;
    pv[i].dSp = dSp;
    pv[i].Ee = Ee;
    pv[i].Ep = Ep;
    pv[i].dEe = dEe;
    pv[i].dEp = dEp;
    cout << endl;
  }
}

void read_scan_info(const char * filename, vector<ScanPoint_t> &pv)
{
  ifstream file(filename);
  if(!file)
  {
    cerr << "Unable to open file with scan information: " << filename << endl;
    exit(1);
  }
  pv.resize(0);
  string s;
  cout << "Reading scan informaion: " << endl;
  cout << setw(3) << "#" << setw(3) << "pnt id" << setw(15) << "W,MeV" << setw(15) << "dW,MeV" << setw(15) << "lum" << setw(16) << "BES runs list"<< endl;
  while(getline(file,s))
  {
    istringstream is(s);
    char c;
    is.get(c); 
    if(c=='#') continue;
    is.putback(c);
    int i; //number of point
    ScanPoint_t sp;
    is >> i >> sp.scan >>  sp.pn >> sp.W >> sp.dW >> sp.lum;
    unsigned run;
    while(is>>run)
    {
      sp.runs.push_back(run);
    }
    pv.push_back(sp);
    cout << setw(3) << i << setw(3) << sp.pn << setw(15) << sp.W << setw(15) << sp.dW << setw(15) << sp.lum;
    for(list<unsigned>::iterator I = sp.runs.begin(); I!=sp.runs.end(); ++I)
    {
      cout << setw(8) << *I;
    }
    cout << endl;
  }
}



void make_scan_points(vector <ScanPoint_t> &pv)
{
  vector<ScanPoint_t> pvtmp;
  read_scan_info("share/scan_info.txt",pvtmp);
  //read_energy("share/cbs-energy2.txt",pvtmp);
  //read_energy("share/cbs-energy3.txt",pvtmp);
  //read_energy("share/cbs-energy4.txt",pvtmp);
  //read_energy("share/cbs-energy5.txt",pvtmp);
  //read_energy("share/cbs-energy6.txt",pvtmp);
  read_energy("share/cbs-energy7.txt",pvtmp);
  int var=2; //combine
  int idx=0;
  switch(var)
  {
    case 0:
      pv.resize(pvtmp.size());
      for(unsigned i=0;i<pvtmp.size();++i,++idx) pv[idx] = pvtmp[i];
      break;
    case 1:
      pv.resize(pvtmp.size()-1);
      for(unsigned i=1;i<pvtmp.size();++i,++idx) pv[idx] = pvtmp[i];
      break;
    case 2:
      cout << "Combining first two run" << endl;
      /*  Combine first points number 1 */
      pv.resize(pvtmp.size()-1);
      pv[0]=pvtmp[0];
      //merge list of runs
      pv[0].runs.merge(pvtmp[1].runs);
      pv[0].lum+=pvtmp[1].lum; //integrated luminosity
      ibn::averager <double> Ea,Ee,Ep;
      Ea.add(pvtmp[0].W, 1./sq(pvtmp[0].dW)*pvtmp[0].lum);
      Ea.add(pvtmp[1].W, 1./sq(pvtmp[1].dW)*pvtmp[1].lum);
      Ee.add(pvtmp[0].Ee, 1./sq(pvtmp[0].dEe)*pvtmp[0].lum);
      Ee.add(pvtmp[1].Ee, 1./sq(pvtmp[1].dEe)*pvtmp[1].lum);
      Ep.add(pvtmp[0].Ep, 1./sq(pvtmp[0].dEp)*pvtmp[0].lum);
      Ep.add(pvtmp[1].Ep, 1./sq(pvtmp[1].dEp)*pvtmp[1].lum);
      pv[0].W = Ea.average();
      pv[0].dW=Ea.sigma_average();
      pv[0].Ee = Ee.average();
      pv[0].dEe=Ee.sigma_average();
      pv[0].Ep = Ep.average();
      pv[0].dEp=Ep.sigma_average();
      unsigned i2=1;
      for(unsigned i=2;i<pvtmp.size();++i,++i2) pv[i2] = pvtmp[i];
      print_data(pv);
      break;
  };
}


void set_selection(int selection_version, TCut & mh_cut, TCut & ee_cut , TCut & gg_cut)
{
  TCut mh_base_cut; //base cut for signal
  TCut mh_strict_cut; //strict cut for signal
  TCut ee_base_cut; //base cut for bhabha
  TCut ee_ext_cut;
  TCut ee_theta_cut;
  TCut mh_theta_cut;
  TCut gg_base_cut;
  TCut gg_theta_cut;
  TCut hp_cut;
  TCut rv_cut;
  TCut ip_cut[3]; //interaction point cut
  TCut mdcEcut;
  TCut ggEcut;
  switch(selection_version)
  {
    case 3:
      mh_base_cut = "nemc>2  && S>0.05 && Eemc<2.5 && Emdc<4";
      mh_strict_cut = "nemc>3  && S>0.05 && Eemc<2.5 && Emdc<4";
      ee_base_cut = "nemc==2 && S<0.05 && Emdc<5 && Eemc>2.5";
      ee_ext_cut = "(nemc==2 || nemc==3) && S<0.05 && Emdc<5 && Eemc>2.5";
      ee_theta_cut  = "Sum$(sin(theta)<0.55)==mdc.ntrack";
      mh_theta_cut  = "Sum$(sin(theta)>0.45)==mdc.ntrack";
      gg_base_cut = "Etotal > 3.3 && Etotal < 4  && sqrt((Sum$(x)-2)**2 + Sum$(y)**2)<4 && abs(Sum$(z))<9";
      gg_theta_cut  = "Sum$(sin(theta)>0.45)==2";
      hp_cut = "Sum$(abs(hpz)<3)==2 && Sum$(hpr<0.25)==2";
      rv_cut = "Sum$(rvxy[hpidx]<0.5)==2 && Sum$(abs(rvz[hpidx])<5)==2";
      /* this is strict cut */
      //mh_cut = mh_strict_cut  && rv_cut && mh_theta_cut && "pt100";
      /* bha bha ext cut */
      //ee_cut = ee_ext_cut  && rv_cut && ee_theta_cut;
      /* very strict cut */
      //mh_cut = mh_strict_cut && rv_cut && mh_theta_cut && "emc.ntrack==0";
      mdcEcut = "Sum$(E>0.05)==nemc";
      ggEcut = "Sum$(E>0.05)==2";

      //mh_cut = mh_base_cut && rv_cut;
      //ee_cut = ee_base_cut && rv_cut &&  ee_theta_cut;
      //gg_cut = gg_base_cut  && gg_theta_cut && ggEcut;

      /* new test cut for E cut */
      //mh_cut = "Sum$(E>0.02)>=3 && S>=0.06 && Eemc<2.5 && Emdc<5"  && rv_cut;
      //ee_cut = "Sum$(E>0.02)==2 && S<=0.05 && Emdc<5 && Eemc>2.5"  && rv_cut && ee_theta_cut;
      gg_cut = gg_base_cut  && gg_theta_cut && "Sum$(E>0.02)==2";

      mh_cut = "ngt > 2 &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut && "nemc==ntrack";
      //mh_cut = "nemc > 2 &&  S>=0.06 && Eemc<2.5 && Emdc<5" && rv_cut && mdcEcut;
      //mh_cut = ("ngt > 2 &&  S>=0.06" || ( "ngt==2 && S>=0.06" && rv_cut && mh_theta_cut)) && "ngt_Eemc<2.5 && Emdc<5";
      ee_cut = "ngt == 2 &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut;
      //ee_cut = "nemc == 2 &&  S<=0.05 && Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut &&mdcEcut;

      //strict cut
      //mh_cut = "ngt >= 4  &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut && mh_theta_cut && "pt100";
      //ee_cut = "ngt == 2  &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut && "emc.ntrack==0";
      //gg_cut = gg_base_cut  && gg_theta_cut;
      break;
    case 4:
      {
        mh_base_cut = "mdc.ntrack>2 && S>0.06 && Eemc<2.5 && Emdc<4 && pt50";
        mh_strict_cut = "mdc.ntrack>3 && S>0.06 && Eemc<2.5 && Emdc<4 && pt100";
        mh_theta_cut  = "Sum$(sin(theta)>0.45)==mdc.ntrack";
        ee_base_cut = "mdc.ntrack==2 && S<0.05 && Eemc>2.5 && Emdc<5";
        ee_theta_cut = "Sum$(sin(theta)<0.45)==mdc.ntrack";
        gg_base_cut = "Etotal > 3.3 && Etotal < 4  && sqrt((Sum$(x)-2)**2 + Sum$(y)**2)<4 && abs(Sum$(z))<9";
        gg_theta_cut = "Sum$(theta)>0.45";
        TCut gg_cos_cut = "acos(cos)>3.1";
        //TCut tof_cut = "Sum$(tof>0)==tof.tof.ntrack && Sum$(tof<10)==tof.tof.ntrack";
        //TCut tof_cut = "Sum$(t0>540 && t0<660)==Length$(t0)";
        TCut tof_cut = "Sum$(tof>1&&tof<6)==Length$(tof) && Sum$(t0>540 && t0<640)==Length$(t0)";

        for(int i=0;i<3;i++)
        {
          char buf[1024];
          double IP_R=0.5;//cm
          double IP_Z=5;//cm
          sprintf(buf,"rvxy[%d]<%f&&abs(rvz[%d])<%f", i, IP_R, i, IP_Z);
          ip_cut[i]=TCut(buf);
        }
        //mh_cut = mh_base_cut  && ip_cut[0] && ip_cut[1] && ip_cut[2]&& mh_theta_cut;
        mh_cut = mh_base_cut  && "Sum$(mdc.rvxy<0.5)==mdc.ntrack && Sum$(abs(mdc.rvz)<5)==mdc.ntrack" && mh_theta_cut;
        ee_cut = ee_base_cut  && ip_cut[0] && ip_cut[1] && ee_theta_cut && "q[0]*q[1]<0";
        gg_cut = gg_base_cut  && gg_theta_cut &&gg_cos_cut;
      }
      break;

    case 6:
      {
        TCut good_track = "Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack && Sum$(mdc.rvxy<1)==mdc.ntrack";
        mh_base_cut  =   "mdc.ntrack>3" && good_track;
        mh_cut = mh_base_cut;

        TCut ee_acol =   "abs(atheta)<0.03 && -0.06 < aphi && aphi<0.01";
        TCut ee_endcup = "abs(cos(mdc.theta[0]))>0.86 && abs(cos(mdc.theta[1]))>0.86";
        TCut ee_E("mdc.E[0]/Eb>0.8 && mdc.E[1]/Eb>0.8 && Emdc<5");
        ee_base_cut = "mdc.ntrack>1 && mdc.ntrack<4" && good_track && ee_endcup && ee_acol && ee_E;
        ee_cut = ee_base_cut  && "q[0]*q[1]<0";

        TCut gg_acol = "abs(atheta) < 0.05 &&  aphi>-0.08 && aphi<0.04";
        TCut gg_barel = "Sum$(cos(theta)<0.8)";
        TCut gg_E("E[0]/Eb && E[1]/Eb>0.8");
        gg_base_cut = gg_acol && gg_barel && gg_E;
        gg_cut = gg_base_cut;
      }
      break;
    case 7:
      {
        TCut good_track = "Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack && Sum$(mdc.rvxy<1)==mdc.ntrack";
        TCut mh_p = "Sum$(p<2.5)==mdc.ntrack";
        TCut mh_2good = "mdc.rvxy[0]<1 && mdc.rvxy[1]<1&& Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack";
        TCut mh_2good2 = "mdc.rvxy[0]<1 && mdc.rvxy[1]<1&& abs(cos(mdc.theta[0]))<0.93 &&abs(cos(mdc.theta[1]))<0.93";
        mh_base_cut  =   "mdc.ntrack>3" && good_track && "S>0.06"&& mh_p;
        //mh_base_cut  =   "mdc.ntrack>2" && good_track; //week cut
        //mh_base_cut  =   "mdc.ntrack>4" && good_track && "S>0.06" && "Sum$(cos(theta)<1.8&&pt>0.05)==mdc.ntrack"; //strong cut
        //mh_base_cut = "mdc.ntrack>2" && "S>0.06" && mh_2good2;
        mh_cut = mh_base_cut;

        TCut ee_acol =   "abs(atheta)<0.03 && -0.06 < aphi && aphi<0.01";
        TCut ee_endcup = "abs(cos(mdc.theta[0]))>0.86 && abs(cos(mdc.theta[1]))>0.86";
        TCut ee_E("mdc.E[0]/Eb>0.8 && mdc.E[1]/Eb>0.8 && mdc.E[0]/Eb<1.2 && mdc.E[1]/Eb<1.2");
        TCut ee_p =  "mdc.p[0]<2.5 && mdc.p[1]<2.5 && mdc.p[0]/Eb>0.9 && mdc.p[1]/Eb>0.9";
        ee_base_cut = "mdc.ntrack>1 && mdc.ntrack<4" && good_track && ee_endcup && ee_acol && ee_E && ee_p;
        ee_cut = ee_base_cut  && "q[0]*q[1]<0";
        
        cout << "Bhabha cut selection: " << endl;
        cout << ee_cut << endl;

        TCut gg_n = "ngct==0 && ngt>1";
        TCut gg_acol = "abs(atheta) < 0.05 &&  aphi>-0.06 && aphi<0.02";
        TCut gg_barel = "abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8";
        TCut gg_E = "E[0]/Eb>0.8 && E[1]/Eb>0.8 && E[0]/Eb<1.2 && E[1]/Eb<1.2";
        gg_base_cut = gg_acol && gg_barel && gg_E && gg_n;
        gg_cut = gg_base_cut;
      }
      break;
  }
}


void make_result(void)
{

  //Cut used in selection;
  TCut  mh_cut;
  TCut  ee_cut;
  TCut  gg_cut;


  int SELECTION_VERSION=7;
  set_selection(SELECTION_VERSION, mh_cut, ee_cut, gg_cut);
  
  list<RunInfo_t> runinfo;
  make_runinfo(runinfo);
  vector <ScanPoint_t> pv;
  make_scan_points(pv);
  //reset luminosity in order to fill it from runinfo table.
  for(unsigned point=0;point<pv.size();++point)
  {
    pv[point].lum=0;
    for(list<unsigned>::iterator i=pv[point].runs.begin();i!=pv[point].runs.end(); ++i)
    {
      for(list<RunInfo_t>::iterator ri=runinfo.begin(); ri!=runinfo.end(); ++ri)
      {
        if(*i==ri->run)
        {
          pv[point].lum+=ri->lum;
          pv[point].ri.push_back(*ri);
        }
      }
    }
  }
  cout << setw(6) << "pnt#" << setw(10) << "Nmh" << setw(10) << "Nee" << setw(10) << "Ngg" << setw(10) << "Nee/Ngg"<< endl;
  for(unsigned pn=0; pn<pv.size(); pn++)
  {
    TChain * mdc = new TChain("mdc","mdc");
    TChain * emc = new TChain("emc","emc");
    TChain * gg  = new TChain("gg","gg");
    //Tune the beem energy
    char Ebeem_str[1024];
    sprintf(Ebeem_str,"(%f*1)",pv[pn].W/2.*1e-3);
    mdc->SetAlias("Eb",Ebeem_str);
    gg->SetAlias("Eb",Ebeem_str);
    //for each run for this scan point
    for(list<unsigned>::iterator i=pv[pn].runs.begin();i!=pv[pn].runs.end(); ++i)
    {
      unsigned run=*i;
      char file_name[1024];
      sprintf(file_name,"psip-%d.root",run);
      mdc->AddFile(file_name);
      emc->AddFile(file_name);
      gg->AddFile(file_name);
    }
    mdc->Draw("ntrack",mh_cut,"goff");
    unsigned Nsignal = mdc->GetSelectedRows();
    mdc->Draw("ntrack",ee_cut,"goff");
    unsigned Nbhabha = mdc->GetSelectedRows();
    gg->Draw("ntrack",gg_cut,"goff");
    unsigned Ngg = gg->GetSelectedRows();
    pv[pn].Nh=Nsignal;
    pv[pn].Nee=Nbhabha;
    pv[pn].Ngg=Ngg;
    cout << setw(6) << pn+1 << setw(10) << Nsignal << setw(10)<< Nbhabha << setw(10)<< Ngg 
      <<  setw(10)<< setprecision(4) << double(Nbhabha)/double(Ngg) << endl;
    delete mdc;
    delete emc;
    delete gg;
  }

  int fw = 14; //format width
  cout.precision(8);
  cout << setw(3) << "#point" << setw(fw) << "lum, nb^-1" 
    << setw(fw) << "W, MeV" << setw(fw)  << "dW, MeV"
    << setw(fw) << "Sw, MeV" << setw(fw)  << "dSw, MeV"
    << setw(fw) << "mhadr" << setw(fw) << "Nee" << setw(fw) << "Ngg" << 
    setw(fw)  << "Nee/lum" << setw(fw) << "Nee/Ngg"<< endl;
  TGraphErrors * nchtr_g = new TGraphErrors;
  TGraphErrors * nntr_g = new TGraphErrors;
  TGraphErrors * bb_lum_g = new TGraphErrors;
  TGraphErrors * bb_gg_g = new TGraphErrors;
  ofstream scan1("scan1.txt");
  ofstream scan2("scan2.txt");
  ofstream scan12("scan12.txt");
  for(unsigned i=0; i<pv.size(); i++)
  {
    ostringstream os;
    os << setw(2) << pv[i].pn << setw(fw) << pv[i].lum  << 
      setw(fw)<< pv[i].W   << setw(fw) << pv[i].dW << 
      setw(fw)<< pv[i].Sw  << setw(fw) << pv[i].dSw << 
      setw(fw)<< pv[i].Nh  << 
      setw(fw)<< pv[i].Nee <<
      setw(fw)<< pv[i].Ngg;
    cout << os.str() << setw(fw) << double(pv[i].Nee)/pv[i].lum << setw(fw) << double(pv[i].Nee)/double(pv[i].Ngg) << endl;
    scan12 << os.str() << endl;
    switch(pv[i].scan)
    {
      case 0:
        break;
      case 1:
        scan1 << os.str()<<endl;
        break;
      case 2:
        scan2 << os.str()<<endl;
        break;
      default:
        break;
    }
    nchtr_g->SetPoint(i,i, pv[i].Nchtr.average());
    nchtr_g->SetPointError(i, 0,pv[i].Nchtr.sigma());

    nntr_g->SetPoint(i,i, pv[i].Nntr.average());
    nntr_g->SetPointError(i, 0,pv[i].Nntr.sigma());
    bb_lum_g->SetPoint(i,pv[i].W, pv[i].Nee/pv[i].lum);
    bb_lum_g->SetPointError(i,pv[i].dW, sqrt(double(pv[i].Nee))/pv[i].lum);
    bb_lum_g->SetTitle("N_{ee} / L");
    bb_lum_g->GetXaxis()->SetTitle("point");
    bb_gg_g->SetPoint(i,pv[i].W, double(pv[i].Nee)/pv[i].Ngg);
    bb_gg_g->SetPointError(i,0, double(pv[i].Nee)/pv[i].Ngg*sqrt(1./pv[i].Nee+1./pv[i].Ngg));
    bb_gg_g->SetTitle("N_{ee} / N_{#gamma#gamma}");
    bb_gg_g->GetXaxis()->SetTitle("W [MeV]");
  }
  //TCanvas * ch_c= new TCanvas("tracks_number","Number of tracks");
  //ch_c->Divide(1,2);
  //ch_c->cd(1);
  //nchtr_g->Draw("a*");
  //nchtr_g->GetXaxis()->SetTitle("point");
  //nchtr_g->GetYaxis()->SetTitle("charged tracks");
  //ch_c->cd(2);
  //nntr_g->Draw("a*");
  //nntr_g->GetXaxis()->SetTitle("point");
  //nntr_g->GetYaxis()->SetTitle("neutral tracks");

  //TCanvas * tr_c = new TCanvas("tracks_number","Number of tracks");
  //tr_c->Divide(1,2);
  //tr_c->cd(1);
  //nchtr2_g->SetMarkerStyle(21);
  //nchtr2_g->Draw("ap");
  //tr_c->cd(2);
  //nntr2_g->SetMarkerStyle(22);
  //nntr2_g->Draw("ap");

  TCanvas * rat_c = new TCanvas;
  rat_c->Divide(1,2);
  rat_c->cd(1);
  gStyle->SetOptFit();
  bb_lum_g->Draw("a*");
  bb_lum_g->SetMarkerStyle(21);
  bb_lum_g->SetMarkerSize(1.5);
  bb_lum_g->SetLineWidth(2);
  bb_lum_g->Fit("pol0");
  TF1 * flum = bb_lum_g->GetFunction("pol0");
  flum->Print();
  double Kee =   flum->GetParameter(0);
  cout << "bb_lum: " << flum->GetParameter(0) << "+-"<< flum->GetParError(0) 
    << " ch2/ndf=" << flum->GetChisquare()/(flum->GetNumberFitPoints()-flum->GetNumberFreeParameters()) << endl;
  rat_c->cd(2);
  TCanvas * ee_gg_c = new TCanvas;
  ee_gg_c->SetFillColor(kWhite);
  bb_gg_g->Draw("ap");
  bb_gg_g->SetMarkerStyle(21);
  bb_gg_g->SetMarkerSize(1.5);
  bb_gg_g->SetLineWidth(2);
  bb_gg_g->Fit("pol0");
  TF1 * fgg = bb_gg_g->GetFunction("pol0");
  fgg->Print();
  double Kbbgg =   fgg->GetParameter(0);
  //double Kgg = Kbbgg/Kee;
  double Kgg = Kee/Kbbgg;
  cout << "Kgg=" << Kgg << endl;
  cout << "bb_gg: " << fgg->GetParameter(0) << "+-"<< fgg->GetParError(0) 
    << " ch2/ndf=" << fgg->GetChisquare()/(fgg->GetNumberFitPoints()-fgg->GetNumberFreeParameters()) << endl;
  //Draw result for different lums
  TMultiGraph * mg =new TMultiGraph;
  TGraphErrors * RLg[2]; 
  TGraphErrors * REEg[2];
  TGraphErrors * RGGg[2];
  TLegend * l = new TLegend(0.6,0.8,1.0,1.0);
  for(int i=0;i<2;i++)
  {
    RLg[i] = new TGraphErrors;
    REEg[i] = new TGraphErrors;
    RGGg[i] = new TGraphErrors;
    RLg[i]->SetMarkerStyle(20+i);
    REEg[i]->SetMarkerStyle(20+i);
    RGGg[i]->SetMarkerStyle(20+i);
    RLg[i]->SetMarkerColor(1);
    REEg[i]->SetMarkerColor(2);
    RGGg[i]->SetMarkerColor(4);
    mg->Add(RLg[i],"p");
    mg->Add(REEg[i],"p");
    mg->Add(RGGg[i],"p");
    char buf[1024];
    sprintf(buf, "scan %d,  online lum",i+1);
    l->AddEntry(RLg[i],buf,"p");
    sprintf(buf, "scan %d,  bhabha lum",i+1);
    l->AddEntry(REEg[i],buf,"p");
    sprintf(buf, "scan %d,  #gamma#gamma lum",i+1);
    l->AddEntry(RGGg[i],buf,"p");
  }
  unsigned Nscan1=6;
  for(unsigned i=0; i<Nscan1; i++)
  {
    RLg[0]->SetPoint(i, pv[i].W, double(pv[i].Nh)/pv[i].lum);
    REEg[0]->SetPoint(i,pv[i].W, double(pv[i].Nh)/pv[i].Nee*Kee);
    RGGg[0]->SetPoint(i,pv[i].W, double(pv[i].Nh)/pv[i].Ngg*Kgg);
    RLg[0]->SetPointError(i, pv[i].dW, sqrt(double(pv[i].Nh))/pv[i].lum);
    REEg[0]->SetPointError(i,pv[i].dW,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
    RGGg[0]->SetPointError(i,pv[i].dW,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
  }
  for(unsigned i=Nscan1; i<pv.size(); i++)
  {
    RLg[1]->SetPoint (i-Nscan1,pv[i].W, double(pv[i].Nh)/pv[i].lum);
    REEg[1]->SetPoint(i-Nscan1,pv[i].W, double(pv[i].Nh)/pv[i].Nee*Kee);
    RGGg[1]->SetPoint(i-Nscan1,pv[i].W, double(pv[i].Nh)/pv[i].Ngg*Kgg);
     RLg[1]->SetPointError(i-Nscan1,pv[i].dW, sqrt(double(pv[i].Nh))/pv[i].lum);
    REEg[1]->SetPointError(i-Nscan1,pv[i].dW,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
    RGGg[1]->SetPointError(i-Nscan1,pv[i].dW,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
  }

  TCanvas * result_c = new TCanvas;
  result_c->SetGridx();
  result_c->SetGridy();
  mg->Draw("a");
  l->Draw();

  ofstream kfile("CrBhabha.txt");
  kfile << Kee << endl;
  kfile << Kgg << endl;
}

