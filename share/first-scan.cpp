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
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TGraphErrors.h>

#include <iomanip>
#include <map>
#include <list>
#include <vector>
#include <iostream>
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
  list <int> runs; //list of runs
  double lum;
  double E;
  double Eerror;
  unsigned long Nh; //number of multihadron
  unsigned long Nee; //number of e+e- (bhabha)
  unsigned long Ngg;  //number of gamma gamma
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

void make_result(void)
{
  const char * signal_cut = "Sum$(E)>-100 && nchtrk>2  &&  S>0.05  && Etotal<2.5";
  //const char * bhabha_cut = "Sum$(E)>-100 && nchtrk==2 &&  S<0.05  && Etotal>2.5 && abs(cos(theta))>0.88";
  const char * bhabha_cut = "Sum$(E)>-100 && nchtrk==2 &&  S<0.05  && Etotal>2.5 &&  sin(theta[0])<0.45 && sin(theta[1])<0.45";
  const char * gg_cut =     "Sum$(E)>-100 && Etotal > 3 && Etotal < 4  && abs(z[1]+z[0])< 1 && abs(x[1]+x[0]-2.5)<5 && abs(y[1]+y[0]) < 5";
  list<RunInfo_t> runinfo;
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

  vector <ScanPoint_t> pv(13);
  pv[0].pn=1;   pv[0].E=3677.872; pv[0].lum=104.561 ;  pv[0].Eerror=0.241;  pv[0].runs.push_back(20334);   pv[0].runs.push_back(20335);
  pv[1].pn=1;   pv[1].E=3678.055; pv[1].lum=464.741 ;  pv[1].Eerror=0.146;  pv[1].runs.push_back(20339);   pv[1].runs.push_back(20336);
  pv[2].pn=3;   pv[2].E=3682.836; pv[2].lum=330.082 ;  pv[2].Eerror=0.141;  pv[2].runs.push_back(20340);   pv[2].runs.push_back(20341);
  pv[3].pn=5;   pv[3].E=3686.298; pv[3].lum=360.3861;  pv[3].Eerror=0.124;  pv[3].runs.push_back(20344);   pv[3].runs.push_back(20342); pv[3].runs.push_back(20343);
  pv[4].pn=7;   pv[4].E=3688.277; pv[4].lum=339.514 ;  pv[4].Eerror=0.192;  pv[4].runs.push_back(20346);   pv[4].runs.push_back(20347);
  pv[5].pn=9;   pv[5].E=3696.883; pv[5].lum=339.000 ;  pv[5].Eerror=0.176;  pv[5].runs.push_back(20350);   pv[5].runs.push_back(20348); pv[5].runs.push_back(20349);
  pv[6].pn=1;   pv[6].E=3676.428; pv[6].lum=230.509 ;  pv[6].Eerror=0.128;  pv[6].runs.push_back(20353);   pv[6].runs.push_back(20354);
  pv[7].pn=2;   pv[7].E=3681.665; pv[7].lum=239.3696;  pv[7].Eerror=0.131;  pv[7].runs.push_back(20357);   pv[7].runs.push_back(20355);
  pv[8].pn=4;   pv[8].E=3683.509; pv[8].lum=302.043 ;  pv[8].Eerror=0.122;  pv[8].runs.push_back(20358);   pv[8].runs.push_back(20359);
  pv[9].pn=5;   pv[9].E=3686.004; pv[9].lum=272.106 ;  pv[9].Eerror=0.126;  pv[9].runs.push_back(20361);   pv[9].runs.push_back(20360);
  pv[10].pn=6; pv[10].E=3687.093; pv[10].lum=283.544 ;  pv[10].Eerror=0.155; pv[10].runs.push_back(20362);  pv[10].runs.push_back(20363);
  pv[11].pn=8; pv[11].E=3690.152; pv[11].lum=325.140 ;  pv[11].Eerror=0.117; pv[11].runs.push_back(20365);  pv[11].runs.push_back(20364);
  pv[12].pn=9; pv[12].E=3693.086; pv[12].lum=354.908 ;  pv[12].Eerror=0.143; pv[12].runs.push_back(20366);  pv[12].runs.push_back(20367);
  cout << setw(10) << "Run #" << setw(20) << "Multihadron" << setw(20) << "Bhabha" << setw(20) << "GammaGamma" << endl;

  //reset luminosity in order to fill it from runinfo table.
  for(int point=0;point<pv.size();++point)
  {
    pv[point].lum=0;
    bool found=false;
    for(list<int>::iterator i=pv[point].runs.begin();i!=pv[point].runs.end(); ++i)
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
  for(unsigned run=20334; run!=20368;++run)
  {
    char buf[1024];
    sprintf(buf,"psip-%d.root",run);
    TFile file(buf);
    if(file.IsOpen())
    {
      TTree * mhadr = (TTree*)file.Get("mhadr");
      TTree * gg = (TTree*)file.Get("gg");
      mhadr->Draw("Etotal",signal_cut,"goff");
      unsigned Nsignal = mhadr->GetSelectedRows();
      mhadr->Draw("Etotal",bhabha_cut,"goff");
      unsigned Nbhabha = mhadr->GetSelectedRows();
      gg->Draw("Etotal",gg_cut,"goff");
      unsigned Ngg = gg->GetSelectedRows();
      cout << setw(10) << run << setw(20) << Nsignal << setw(20)<< Nbhabha << setw(20)<< Ngg <<  setw(20)<< double(Nbhabha)/double(Ngg) << endl;

      for(unsigned pn=0; pn<13; pn++)
      {
        bool found=false;
        for(list<int>::iterator i=pv[pn].runs.begin();i!=pv[pn].runs.end(); ++i)
        {
          if(*i==run) found=true;
        }
        if(found) 
        {
          pv[pn].Nh+=Nsignal;
          pv[pn].Nee+=Nbhabha;
          pv[pn].Ngg+=Ngg;
        }
      }
    }
  }

  cout.precision(8);
  cout << setw(3) << "#point" << setw(20) << "lum, nb^-1" << setw(20) << "energy, MeV" << setw(20)  << "error, MeV" << setw(20) << "signal (mhadr)" << setw(20) << "bhabha" << setw(20) << "gamma-gamma" << endl;
  TGraphErrors * result_g = new TGraphErrors;
  for(unsigned i=0; i<pv.size(); i++)
  {
    cout << setw(3) << pv[i].pn << setw(20) << pv[i].lum  << 
      setw(20)<< pv[i].E   << setw(20) << pv[i].Eerror << 
      setw(20)<< pv[i].Nh  << 
      setw(20)<< pv[i].Nee <<
      setw(20)<< pv[i].Ngg << endl;
    result_g->SetPoint(i,pv[i].E, double(pv[i].Nh)/pv[i].Nee);
    result_g->SetPointError(i,pv[i].Eerror, 0);
  }
  result_g->Draw("ap");
  result_g->SetMarkerStyle(21);
}

