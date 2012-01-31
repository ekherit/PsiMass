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
#include "utils.h"
#include "selection.h"
#include "interference.h"
#include "RunInfo.h"
#include "energy.h"
using namespace std;




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
	double lum_ee_cor; //correction to Nee luminosity
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

void set_alias(TTree * )
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

bool ScanPoint_energy_comp(const ScanPoint_t & sp1, const ScanPoint_t &sp2)
{
  return sp1.W < sp2.W;
}

void combine(double dW, list<ScanPoint_t> & spl, list<ScanPoint_t> &cmb)
{
  ibn::averager<double> Wa; //averager for the cm energy
  ibn::averager<double> Sa; //averager for the energy spread
  double lastW=0,lastdW=0;
  ScanPoint_t sp;
  unsigned pn=0;
  for(list<ScanPoint_t>::iterator i=spl.begin();i!=spl.end();++i)
  {
    double W = i->W;
    if(i!=spl.begin() && abs(W-Wa.average())>dW)
    {
      sp.W = Wa.waverage();
      sp.dW=  Wa.wsigma_average();
      sp.Sw = Sa.waverage();
      sp.dSw = Sa.wsigma_average();
      sp.pn=pn++;
      cmb.push_back(sp);
      Wa.reset();
      ScanPoint_t tmpsp;
      sp=tmpsp;
    }
    if(i->W!=lastW && i->dW != lastdW) 
    {
      Wa.add(W,1./sq(i->dW));
      Sa.add(i->Sw,1./sq(i->dSw));
    }
    lastW=i->W;
    lastdW=i->dW;
    sp.runs.insert(sp.runs.end(),i->runs.begin(),i->runs.end());
    sp.ri.insert(sp.ri.end(),i->ri.begin(),i->ri.end());
    sp.lum+=i->lum;
    sp.Nh+=i->Nh;
    sp.Nee+=i->Nee;
    sp.Ngg+=i->Ngg;
  }
  sp.W = Wa.waverage();
  sp.dW= Wa.wsigma_average();
  sp.Sw = Sa.waverage();
  sp.dSw = Sa.wsigma_average();
  sp.pn=pn++;
  cmb.push_back(sp);
  Wa.reset();
  ScanPoint_t tmpsp;
  sp=tmpsp;
}

void make_scan_points2(const char * runinfo_filename, vector <ScanPoint_t> &pv)
{
  cout << "Reading run info" << endl;
  std::vector<RunInfo_t> RI;
  read_run_info(runinfo_filename, RI);
  cout << "Done" << endl;
  pv.resize(RI.size());
  for(unsigned i=0; i<pv.size(); ++i)
  {
    pv[i].pn=i;
    pv[i].scan=1; //scan number the same
    pv[i].lum = RI[i].lum;
    pv[i].runs.push_back(RI[i].run); //only one run
    pv[i].Ee=RI[i].e.E;
    pv[i].dEe=RI[i].e.dE;
    pv[i].Ep=RI[i].p.E;
    pv[i].dEp=RI[i].p.dE;
    pv[i].Se=RI[i].e.S;
    pv[i].dSe=RI[i].e.dS;
    pv[i].Sp=RI[i].p.S;
    pv[i].dSp=RI[i].p.dS;
    pv[i].W = RI[i].W;
    pv[i].dW = RI[i].dW;
    //pv[i].W = cm_energy(pv[i].Ee, pv[i].Ep);
    //pv[i].dW = sqrt(sq(dW_dE1(pv[i].Ee,pv[i].Ep)*pv[i].dEe) + sq(dW_dE1(pv[i].Ep,pv[i].Ee)*pv[i].dEp));
    pv[i].Sw = sqrt(pv[i].Sp*pv[i].Sp+pv[i].Se*pv[i].Se);
    pv[i].dSw = sqrt(sq(pv[i].Sp*pv[i].dSp)+sq(pv[i].Se*pv[i].dSe))/pv[i].Sw;
    cout << i << " run="<< RI[i].run << " ";
    cout << "W="<<pv[i].W << "+-" << pv[i].dW;
    cout << "lum=" << pv[i].lum;
    cout << endl;
  }
  cout << "end of make scan points " << endl;
  list<ScanPoint_t> spl;
  for(unsigned i=0;i<pv.size();++i)
  {
    spl.push_back(pv[i]);
  }
  spl.sort(ScanPoint_energy_comp);
  //for(list<ScanPoint_t>::iterator i=spl.begin();i!=spl.end();++i)
  //{
  //  cout << "sorted: " << i->runs.front() << " " << i->W << "+-" << i->dW << endl;
  //}
  list <ScanPoint_t> cmb;
  combine(0.2,spl,cmb);
  unsigned pn=0;
  for(list<ScanPoint_t>::iterator i=cmb.begin();i!=cmb.end();++i)
  {
    cout << setw(3)<< ++pn<< " ";
    unsigned prevrun=0;
    i->runs.sort();
    for(list<unsigned>::iterator r=i->runs.begin();r!=i->runs.end();++r)
    {
      if(*r==(prevrun+1)) 
      {
        if(++r!=i->runs.end())
        {
          cout << "-";
          r--;
        }
        else 
        {
          r--;
          cout << *r;
        }
      }
      else 
      {
        if(prevrun==0)
          cout <<*r;
        else 
          cout << prevrun<<","<<*r;
      }
      prevrun=*r;
    }
    cout << " "  << setprecision(7)<< i->W <<  "+-" << setprecision(4)<<i->dW << " " << i->lum<<endl;
  }
  pv.resize(cmb.size());
  pn=0;
  for(list<ScanPoint_t>::iterator i=cmb.begin();i!=cmb.end();++i)
  {
    pv[pn] = *i;
    pn++;
  }
}


void draw_energy_vs_time(const char * runinfo_filename)
{
  vector <RunInfo_t> pv;
  read_run_info(runinfo_filename,pv);
  TGraphErrors * besrung=new TGraphErrors(pv.size());
  TGraphErrors * emsruneg=new TGraphErrors(pv.size());
  TGraphErrors * emsrunpg=new TGraphErrors(pv.size());
  for(unsigned i=0;i<pv.size();++i)
  {
    double t = pv[i].begin_time/2.+pv[i].end_time/2.;
    double dt = pv[i].end_time/2. - pv[i].begin_time/2.;
    cout << t << " " << dt << endl;
    besrung->SetPoint(i,t,pv[i].W/2.+0.1);
    besrung->SetPointError(i,dt,pv[i].dW/2.);
    t = pv[i].e.begin_time/2.+pv[i].e.end_time/2.;
    dt = pv[i].e.end_time/2. - pv[i].e.begin_time/2.;
    cout << t << " " << dt << " " << pv[i].e.begin_time << " " << pv[i].e.end_time << endl;
    emsruneg->SetPoint(i,t,pv[i].e.E);
    emsruneg->SetPointError(i,dt,pv[i].e.dE);
    t = pv[i].p.begin_time/2.+pv[i].p.end_time/2.;
    dt = pv[i].p.end_time/2. - pv[i].p.begin_time/2.;
    cout << t << " " << dt << endl;
    emsrunpg->SetPoint(i,t,pv[i].p.E);
    emsrunpg->SetPointError(i,dt,pv[i].p.dE);
  }
  besrung->SetMarkerStyle(20);
  emsruneg->SetMarkerStyle(21);
  emsrunpg->SetMarkerStyle(22);
  besrung->Draw("ap");
  besrung->GetXaxis()->SetTimeDisplay(1);
  besrung->GetXaxis()->SetTimeFormat("%H:%M%F1970-01-01");
  emsruneg->Draw("p");
  emsruneg->SetMarkerColor(kBlue);
  emsrunpg->SetMarkerColor(kRed);
  emsrunpg->Draw("p");
}

void make_result(const char * runinfo_filename)
{
  //Cut used in selection;
  TCut  mh_cut;
  TCut  ee_cut;
  TCut  gg_cut;


  int SELECTION_VERSION=7;
  set_selection(SELECTION_VERSION, mh_cut, ee_cut, gg_cut);
  ofstream cuts_file("cuts.txt");
  cuts_file << "Multihadronic cut: " << mh_cut << endl;
  cuts_file << "Bhabha cut: " << ee_cut << endl;
  cuts_file << "Gamma gamma cut: " << gg_cut << endl;
  
  list<RunInfo_t> runinfo;
  make_runinfo(runinfo);
  vector <ScanPoint_t> pv;
  make_scan_points2(runinfo_filename, pv);
  //reset luminosity in order to fill it from runinfo table.
  //for(unsigned point=0;point<pv.size();++point)
  //{
  //  pv[point].lum=0;
  //  for(list<unsigned>::iterator i=pv[point].runs.begin();i!=pv[point].runs.end(); ++i)
  //  {
  //    for(list<RunInfo_t>::iterator ri=runinfo.begin(); ri!=runinfo.end(); ++ri)
  //    {
  //      if(*i==ri->run)
  //      {
  //        pv[point].lum+=ri->lum;
  //        pv[point].ri.push_back(*ri);
  //      }
  //    }
  //  }
  //}
  cout << setw(6) << "pnt#" << setw(8)<< "lum"<<setw(10) << "Nmh" << setw(10) << "Nee" << setw(10) << "Ngg" << setw(10) << "Nee/Ngg"<< endl;
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
      sprintf(file_name,"data/%d.root",run);
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
		pv[pn].lum_ee_cor = BBIntCor(pv[pn].W);
    pv[pn].Nee=Nbhabha;
    pv[pn].Ngg=Ngg;
    cout << setw(6) << pn+1<< setw(8) << pv[pn].lum << setw(10) << Nsignal << setw(10)<< Nbhabha << setw(10)<< Ngg 
      <<  setw(10)<< setprecision(4) << double(Nbhabha)/double(Ngg) 
      <<  setw(8) << (pv[pn].lum_ee_cor-1)*100.; ;
    //cout << setw(10) << pv[pn].Sw<<"+-"<<pv[pn].dSw;
    cout << setw(10) << pv[pn].W<<"+-"<<pv[pn].dW;
    cout << endl;
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
      setw(fw)<< pv[i].Ngg <<
			setw(fw)<< pv[i].lum_ee_cor;
    cout << os.str() << setw(fw) << double(pv[i].Nee)/pv[i].lum << setw(fw) << double(pv[i].Nee)/double(pv[i].Ngg);
    cout << setw(fw) << pv[i].Sw<<"+-"<<pv[i].dSw;
    cout << endl;
    scan12 << os.str() << endl;
    switch(pv[i].scan)
    {
      default:
        scan1 << os.str()<<endl;
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
  for(int i=0;i<1;i++)
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
  unsigned Nscan1=100;
  for(unsigned i=0; i<pv.size(); i++)
  {
    RLg[0]->SetPoint(i, pv[i].W/2, double(pv[i].Nh)/pv[i].lum);
    REEg[0]->SetPoint(i,pv[i].W/2, double(pv[i].Nh)/pv[i].Nee*Kee);
    RGGg[0]->SetPoint(i,pv[i].W/2, double(pv[i].Nh)/pv[i].Ngg*Kgg);
    RLg[0]->SetPointError(i, pv[i].dW/2, sqrt(double(pv[i].Nh))/pv[i].lum);
    REEg[0]->SetPointError(i,pv[i].dW/2,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
    RGGg[0]->SetPointError(i,pv[i].dW/2,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
  }
  //for(unsigned i=Nscan1; i<pv.size(); i++)
  //{
  //  RLg[1]->SetPoint (i-Nscan1,pv[i].W/2, double(pv[i].Nh)/pv[i].lum);
  //  REEg[1]->SetPoint(i-Nscan1,pv[i].W/2, double(pv[i].Nh)/pv[i].Nee*Kee);
  //  RGGg[1]->SetPoint(i-Nscan1,pv[i].W/2, double(pv[i].Nh)/pv[i].Ngg*Kgg);
  //   RLg[1]->SetPointError(i-Nscan1,pv[i].dW/2, sqrt(double(pv[i].Nh))/pv[i].lum);
  //  REEg[1]->SetPointError(i-Nscan1,pv[i].dW/2,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
  //  RGGg[1]->SetPointError(i-Nscan1,pv[i].dW/2,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
  //}

  TCanvas * result_c = new TCanvas;
  result_c->SetGridx();
  result_c->SetGridy();
  mg->Draw("a");
  l->Draw();

  ofstream kfile("CrBhabha.txt");
  kfile << Kee << endl;
  kfile << Kgg << endl;
}

