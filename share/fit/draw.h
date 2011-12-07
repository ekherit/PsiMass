/*
 * =====================================================================================
 *
 *       Filename:  utils.h
 *
 *    Description:  Различные вспомогательные функции.
 *
 *        Version:  1.0
 *        Created:  05.12.2011 18:20:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_FITSCAN_UTILS_H
#define IBN_FITSCAN_UTILS_H

#include <vector>

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>

#include "FitOniumRCompactLib.hh"


using namespace std;


inline void draw_res(int resid, double sW, std::vector <double> & Es)
{
  TCanvas * psip_c  = new TCanvas;
  unsigned NPX=100; //points used to draw resonance curve
  TGraph * fun = new TGraph(NPX); //resonance curve
  TGraph * graph = new TGraph;   //layout points
  /* draw function */
  vector<double> par(4);  
  double M=0, Emax=0, Emin=0;
  par[0]=100; //background
  par[1]=1;  //efficiency
  par[2]=0;  //delta M/2, MeV
  par[3]=sW; //sigmaW, MeV
  switch (resid)
  {
    case _IdPsiPrime:
      M = _MPsiPrime;
      break;
    case _IdJPsi:
      M = _MJPsi;
      break;
  }
  Emax = M/2+sW*5;
  Emin = M/2-sW*5;
  for(int i=0;i<NPX;i++)
  {
    double E = Emin+(Emax-Emin)/NPX*i;
    fun->SetPoint(i, E, CrSOniumR(_MethodAzimov, resid, E , &par[0]));
  }

  for(int i=0;i<Es.size();i++)
  {
    graph->SetPoint(i,Es[i],CrSOniumR(_MethodAzimov, resid, Es[i], &par[0]));
        //cout << EPSI[i] << " " 
        //  << CrSOniumR(_MethodAzimov, _IdPsiPrime, EPSI[i], scen_psi2s_par) 
        //  << " " <<scen_psi2s->EvalPar(&EPSI[i], scen_psi2s_par) 
        //  << " " <<scen_psi2s->Eval(EPSI[i])
        //  << endl;
  }
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(3);
  graph->SetMarkerColor(kRed);
  fun->Draw("ac");
  fun->SetLineWidth(2);
  graph->Draw("p");
  fun->GetXaxis()->SetTitle("E, MeV");
}
#endif
