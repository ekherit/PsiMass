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
struct ScanPoint_t
{
  list <unsigned> runs; //list of runs
  double lum;
  double E;
  double Eerror;
  unsigned long Nh; //number of multihadron
  unsigned long Nee; //number of e+e- (bhabha)
  unsigned long Ngg;  //number of gamma gamma
  ibn::averager<double> Nchtr, Nntr; //number of charged tracks and it rms
  unsigned pn;
  ScanPoint_t(void)
  {
    Nh=0;
    Nee=0;
    Ngg=0;
    lum=0;
    E=0;
    Eerror=0;
    pn=0;
  }
};

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
  t->SetAlias("ecut", "mdc.E>=0.02");
  t->SetAlias("ngt","Sum$(mdc.E>=0.02)");
  t->SetAlias("ngt_Eemc","Sum$((mdc.E>=0.02)*E)");
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


void make_scan_points(vector <ScanPoint_t> &pv)
{
  pv.resize(13);
  ScanPoint_t * sp=&pv[0];
  sp=&pv[0]; sp->pn=1; sp->E=3677.872; sp->lum=104.561; sp->Eerror=0.241; sp->runs.push_back(20334); sp->runs.push_back(20335); 
  sp=&pv[1]; sp->pn=1; sp->E=3678.055; sp->lum=464.741;sp->Eerror=0.146;sp->runs.push_back(20339);sp->runs.push_back(20336);
  sp=&pv[2]; sp->pn=3; sp->E=3682.836; sp->lum=330.082;sp->Eerror=0.141;sp->runs.push_back(20340);sp->runs.push_back(20341);
  sp=&pv[3]; sp->pn=5; sp->E=3686.298; sp->lum=360.3861;sp->Eerror=0.124;sp->runs.push_back(20344);sp->runs.push_back(20342); sp->runs.push_back(20343); 
  sp=&pv[4]; sp->pn=7; sp->E=3688.277; sp->lum=339.514;sp->Eerror=0.192;sp->runs.push_back(20346);sp->runs.push_back(20347); 
  sp=&pv[5]; sp->pn=9; sp->E=3696.883; sp->lum=339.000;sp->Eerror=0.176;sp->runs.push_back(20350);sp->runs.push_back(20348); sp->runs.push_back(20349); 
  sp=&pv[6]; sp->pn=1; sp->E=3676.428; sp->lum=230.509;sp->Eerror=0.128;sp->runs.push_back(20353);sp->runs.push_back(20354); 
  sp=&pv[7]; sp->pn=2; sp->E=3681.665; sp->lum=239.3696;sp->Eerror=0.131;sp->runs.push_back(20357);sp->runs.push_back(20355);
  sp=&pv[8]; sp->pn=4; sp->E=3683.509; sp->lum=302.043;sp->Eerror=0.122;sp->runs.push_back(20358);sp->runs.push_back(20359); 
  sp=&pv[9]; sp->pn=5; sp->E=3686.004; sp->lum=272.106;sp->Eerror=0.126;sp->runs.push_back(20361);sp->runs.push_back(20360); 
  sp=&pv[10]; sp->pn=6; sp->E=3687.093; sp->lum=283.544; sp->Eerror=0.155; sp->runs.push_back(20362);sp->runs.push_back(20363);
  sp=&pv[11]; sp->pn=8; sp->E=3690.152; sp->lum=325.140;sp->Eerror=0.117; sp->runs.push_back(20365);sp->runs.push_back(20364);
  sp=&pv[12]; sp->pn=9; sp->E=3693.086; sp->lum=354.908;sp->Eerror=0.143; sp->runs.push_back(20366);sp->runs.push_back(20367);
}

void make_result(void)
{
  TCut mh_base_cut = "nemc>2  && S>0.05 && Eemc<2.5 && Emdc<4";
  TCut mh_strict_cut = "nemc>3  && S>0.05 && Eemc<2.5 && Emdc<4";
  TCut ee_base_cut = "nemc==2 && S<0.05 && Emdc<5 && Eemc>2.5";
  TCut ee_ext_cut = "(nemc==2 || nemc==3) && S<0.05 && Emdc<5 && Eemc>2.5";
  TCut ee_theta_cut  = "Sum$(sin(theta)<0.45)==mdc.ntrack";
  TCut mh_theta_cut  = "Sum$(sin(theta)>0.45)==mdc.ntrack";
  TCut gg_base_cut = "Etotal > 3.3 && Etotal < 4  && sqrt((Sum$(x)-2)**2 + Sum$(y)**2)<4 && abs(Sum$(z))<9";
  TCut gg_theta_cut  = "Sum$(sin(theta)>0.45)==2";
  TCut hp_cut = "Sum$(abs(hpz)<3)==2 && Sum$(hpr<0.25)==2";
  TCut rv_cut = "Sum$(rvxy[hpidx]<0.5)==2 && Sum$(abs(rvz[hpidx])<5)==2";
  TCut mh_ecut = "Sum$(E>=0.02)==nemc"; //energy deposition in calorimeter more then 20MeV
  /* this is standart cut */
  TCut  mh_cut = mh_base_cut  && rv_cut;
  TCut  ee_cut = ee_base_cut  && rv_cut && ee_theta_cut;
  TCut  gg_cut = gg_base_cut  && gg_theta_cut;
  /* this is strict cut */
  //mh_cut = mh_strict_cut  && rv_cut && mh_theta_cut && "pt100";
  /* bha bha ext cut */
  //ee_cut = ee_ext_cut  && rv_cut && ee_theta_cut;
  /* very strict cut */
  //mh_cut = mh_strict_cut && rv_cut && mh_theta_cut && "emc.ntrack==0";
  /* new test cut for E cut */
  //mh_cut = "Sum$(E>0.02)>2 && S>0.05 && Eemc<2.5 && Emdc<4"  && rv_cut && "pt50" && mh_theta_cut;
  //ee_cut = "Sum$(E>0.02)==2 && S<0.05 && Emdc<5 && Eemc>2.5"  && rv_cut && ee_theta_cut;
  //gg_cut = gg_base_cut  && gg_theta_cut && "Sum$(E>0.02)==2";

  mh_cut = "ngt >= 3 &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut;
  ee_cut = "ngt == 2 &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut;

  //mh_cut = "ngt > 3 &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut && mh_theta_cut && "pt100";
  //ee_cut = "ngt == 2 &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut;
  
  list<RunInfo_t> runinfo;
  make_runinfo(runinfo);
  vector <ScanPoint_t> pv;
  make_scan_points(pv);
  cout << setw(10) << "run #" << setw(20) << "multihadron" << setw(20) << "bhabha" << setw(20) << "gammagamma" << endl;
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
        }
      }
    }
  }
//temporary
  unsigned runidx=0;
  TGraphErrors * nchtr2_g = new TGraphErrors;
  TGraphErrors * nntr2_g = new TGraphErrors;
  for(unsigned run=20334; run!=20368;++run)
  {
    char buf[1024];
    sprintf(buf,"psip-%d.root",run);
    TFile file(buf);
    if(file.IsOpen())
    {
      TTree * mhadr = (TTree*)file.Get("mhadr");
      TTree * mdc = (TTree*)file.Get("mdc");
      TTree * emc = (TTree*)file.Get("emc");
      mhadr->AddFriend(mdc);
      mhadr->AddFriend(emc);
      set_alias(mhadr);
      TTree * gg = (TTree*)file.Get("gg");
      mhadr->Draw("Etotal",mh_cut,"goff");
      unsigned Nsignal = mhadr->GetSelectedRows();
      mhadr->Draw("Etotal",ee_cut,"goff");
      unsigned Nbhabha = mhadr->GetSelectedRows();
      //draw  number of charged tracks
      mhadr->Draw("mdc.ntrack>>hnchtr",mh_cut && mh_theta_cut && "emc.ntrack==0","goff");
      TH1F * hnchtr = (TH1F*)gDirectory->Get("hnchtr");
      hnchtr->Fit("gaus","IQ0");
      TF1 * nchtr_f = hnchtr->GetFunction("gaus");
      double nchtr=0, nntr=0, nchtr_rms=1, nntr_rms=1;
      nchtr=nchtr_f->GetParameter(1);
      double nchtr_sigma = nchtr_f->GetParError(1)*sqrt(nchtr_f->GetChisquare()/nchtr_f->GetNDF());
      //nchtr2_g->SetPoint(runidx,run,nchtr);
      //nchtr2_g->SetPointError(runidx,0,nchtr_sigma);
      
      //nchtr2_g->SetPoint(runidx,run,hnchtr->GetMean());
      //nchtr2_g->SetPointError(runidx,0,hnchtr->GetRMS()/sqrt(mhadr->GetSelectedRows()));
      nchtr = hnchtr->GetMean(); //number of charged tracks
      nchtr_rms = hnchtr->GetRMS()/sqrt(hnchtr->GetEntries()); //RMS of charged tracks.
      mhadr->Draw("emc.ntrack>>hnntr",mh_cut && mh_theta_cut,"goff");
      TH1F * hnntr = (TH1F*)gDirectory->Get("hnntr");
      hnntr->Fit("gaus","IQ0");
      TF1 * nntr_f = hnntr->GetFunction("gaus");
      nntr = hnntr->GetMean(); //number of neutral tracks
      nntr_rms = hnntr->GetRMS()/sqrt(hnntr->GetEntries()); //RMS of number of neutral tracks.
      
      gg->Draw("Etotal",gg_cut,"goff");
      unsigned Ngg = gg->GetSelectedRows();
      cout << setw(10) << run << setw(20) << Nsignal << setw(20)<< Nbhabha << setw(20)<< Ngg <<  setw(20)<< double(Nbhabha)/double(Ngg) << endl;

      for(unsigned pn=0; pn<13; pn++)
      {
        bool found=false;
        for(list<unsigned>::iterator i=pv[pn].runs.begin();i!=pv[pn].runs.end(); ++i)
        {
          if(*i==run) found=true;
        }
        if(found) 
        {
          pv[pn].Nh+=Nsignal;
          pv[pn].Nee+=Nbhabha;
          pv[pn].Ngg+=Ngg;
          pv[pn].Nchtr.add(nchtr,1./sq(nchtr_rms));
          pv[pn].Nntr.add(nntr,1./sq(nntr_rms));
          nchtr2_g->SetPoint(runidx,pv[pn].E,nchtr);
          nchtr2_g->SetPointError(runidx,pv[pn].Eerror,nchtr_rms);
          nntr2_g->SetPoint(runidx,pv[pn].E,nntr);
          nntr2_g->SetPointError(runidx,pv[pn].Eerror,nntr_rms);
          runidx++;
        }
      }
    }
  }



  cout.precision(8);
  cout << setw(3) << "#point" << setw(20) << "lum, nb^-1" << setw(20) << "energy, MeV" << setw(20)  << "error, MeV" << setw(20) << "signal (mhadr)" << setw(20) << "bhabha" << setw(20) << "gamma-gamma" << endl;
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
    os << setw(2) << pv[i].pn << setw(20) << pv[i].lum  << 
      setw(20)<< pv[i].E   << setw(20) << pv[i].Eerror << 
      setw(20)<< pv[i].Nh  << 
      setw(20)<< pv[i].Nee <<
      setw(20)<< pv[i].Ngg << endl;
    cout << os.str();
    scan12 << os.str();
    if(i<6) scan1 << os.str();
    else scan2 << os.str();
    
    nchtr_g->SetPoint(i,i, pv[i].Nchtr.average());
    nchtr_g->SetPointError(i, 0,pv[i].Nchtr.sigma());

    nntr_g->SetPoint(i,i, pv[i].Nntr.average());
    nntr_g->SetPointError(i, 0,pv[i].Nntr.sigma());
    bb_lum_g->SetPoint(i,i+1, pv[i].Nee/pv[i].lum);
    bb_lum_g->SetTitle("N_{ee} / L");
    bb_lum_g->GetXaxis()->SetTitle("point");
    bb_lum_g->SetPointError(i,0, sqrt(double(pv[i].Nee))/pv[i].lum);
    bb_gg_g->SetPoint(i,i+1, double(pv[i].Nee)/pv[i].Ngg);
    bb_gg_g->SetPointError(i,0, double(pv[i].Nee)/pv[i].Ngg*sqrt(1./pv[i].Nee+1./pv[i].Ngg));
    bb_gg_g->SetTitle("N_{ee} / N_{#gamma#gamma}");
    bb_gg_g->GetXaxis()->SetTitle("point");
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

  TCanvas * tr_c = new TCanvas("tracks_number","Number of tracks");
  tr_c->Divide(1,2);
  tr_c->cd(1);
  nchtr2_g->SetMarkerStyle(21);
  nchtr2_g->Draw("ap");
  tr_c->cd(2);
  nntr2_g->SetMarkerStyle(22);
  nntr2_g->Draw("ap");

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
    RLg[0]->SetPoint(i,pv[i].E, double(pv[i].Nh)/pv[i].lum);
    REEg[0]->SetPoint(i,pv[i].E, double(pv[i].Nh)/pv[i].Nee*Kee);
    RGGg[0]->SetPoint(i,pv[i].E, double(pv[i].Nh)/pv[i].Ngg*Kgg);
    RLg[0]->SetPointError(i,pv[i].Eerror, sqrt(double(pv[i].Nh))/pv[i].lum);
    REEg[0]->SetPointError(i,pv[i].Eerror,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
    RGGg[0]->SetPointError(i,pv[i].Eerror,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
  }
  for(unsigned i=Nscan1; i<pv.size(); i++)
  {
    RLg[1]->SetPoint(i-Nscan1,pv[i].E, double(pv[i].Nh)/pv[i].lum);
    REEg[1]->SetPoint(i-Nscan1,pv[i].E, double(pv[i].Nh)/pv[i].Nee*Kee);
    RGGg[1]->SetPoint(i-Nscan1,pv[i].E, double(pv[i].Nh)/pv[i].Ngg*Kgg);
     RLg[1]->SetPointError(i-Nscan1,pv[i].Eerror, sqrt(double(pv[i].Nh))/pv[i].lum);
    REEg[1]->SetPointError(i-Nscan1,pv[i].Eerror,double(pv[i].Nh)/pv[i].Nee*Kee*sqrt(1./pv[i].Nh+1./pv[i].Nee));
    RGGg[1]->SetPointError(i-Nscan1,pv[i].Eerror,double(pv[i].Nh)/pv[i].Ngg*Kgg*sqrt(1./pv[i].Nh+1./pv[i].Ngg));
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

