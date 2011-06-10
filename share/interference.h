/*
 * =====================================================================================
 *
 *       Filename:  interference.h
 *
 *    Description:  Discribe the interference of continuum and resonance in bhabha scattering
 *
 *        Version:  1.0
 *        Created:  05/27/2011 05:18:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_BHABHA_INTERFERENCE_H
#define IBN_BHABHA_INTERFERENCE_H

#include "utils.h"
#include <iostream>

const double MPDG=3686.09;
const double PI = 3.1415926535897;
const double ALPHA = 1./136.03599911;
const double GAMMAEE_PSI2S = 0.00235; //MeV
const double GAMMA_PSI2S = 0.304; // MeV
const double M_PSI2S  = 3686.09; //MeV
const double ME = 0.510998910; //MeV
const double HC2 = 0.389379304*1e6*1e6; //MeV^2*nbarn

#include <TF1.h>
inline double sigma(double *x,  double *p)
{
	double W = MPDG+x[0];
	if(W<=0) return 1e100;
	double G = p[3]; //width
  double den=sqrt(sq(W-MPDG)+sq(G/2.));
  //I must be accurate with calculating the phase
	double a = atan(-G/2./(W-MPDG));
  //if atan is negative it correspond angle > 90degre
  if(a<0) a+=PI;
	double beta = 4*ALPHA/PI*(log(W/ME)-0.5);
	double beta0 = 4*ALPHA/PI*(log((2*MPDG)/ME)-0.5);
  double bsin = PI*beta/sin(PI*beta);
  double bsin0 = PI*beta0/sin(PI*beta0);
	//resonance beheavioure
	double A=bsin/bsin0*pow(G/(2*den),1.-beta);
	//interference contribution max
	double C = p[0]*sq(MPDG/W); //continuum
	double I = p[1]*A*cos(a*(1-beta)); 
	double R = p[2]*A*sin(a*(1-beta)); 
	return C+I+R;
}

//energy spread
inline double spread(double dW,   double S)
{
	return 1./(sqrt(2*PI)*S)*exp(-0.5*sq(dW/S));
}





double ee_interference(double W, double theta);

double ee_interference(double W, double cos_min, double cos_max);

inline double ee_interference_vs_W(double *x, double *p)
{
    return  ee_interference(x[0], p[0], p[1]);
}



inline double cs_bhabha_with_interference(double *x, double *p)
{
  double W = x[0]+M_PSI2S;
  double QED = p[0];
  double INT = p[1];
  //this parameters should be fixed
  double cmin  = p[2];
  double cmax = p[3]; 
	double C = QED*sq(M_PSI2S/W); //continuum
  double I = INT*2.0*ee_interference(W, cmin, cmax); 
  /*
   *  The 2.0 must be becase integration over angle was made only for one side
   *  cmin and cmax here is absolute values of the cos theta.
   */
  return C+I;
}

inline double interference_correction(double W, double QED, double INT, double cmin, double cmax)
{
	double C = QED*sq(M_PSI2S/W); //continuum
  double I = INT*2.0*ee_interference(W, cmin, cmax); 
  return I/C;
}

inline double interference_correction(double *x, double *p)
{
  double W = x[0]+M_PSI2S;
  double QED = p[0];
  double INT = p[1];
  //this parameters should be fixed
  double cmin  = p[2];
  double cmax = p[3]; 
  return interference_correction(W,QED,INT,cmin,cmax);
	//double C = QED*sq(M_PSI2S/W); //continuum
  //double I = INT*2.0*ee_interference(W, cmin, cmax); 
  ///*
  // *  The 2.0 must be becase integration over angle was made only for one side
  // *  cmin and cmax here is absolute values of the cos theta.
  // */
  //return I/C;
}

inline double ee_interference2_vs_W(double *x, double *p)
{
    return  ee_interference(x[0], p[0]);
}

inline double ee_interference_vs_theta(double *x, double *p)
{
    return  ee_interference(p[0], x[0]);
}

double ee_interference(double W, double theta_min, double theta_max);

double ee_interference_spreaded(double W, double spread, double eff, double cmin, double cmax);

double ee_cs_qed_spreaded(double W, double spread, double qed);

inline double ee_interference_spreaded(double *x, double *p)
{
  double W = x[0]+M_PSI2S;
  double SPREAD=p[0];
  double EFF = p[1];
  double CMIN = p[2];
  double CMAX = p[3];
  return ee_interference_spreaded(W,SPREAD,EFF, CMIN,CMAX);
}


inline double ee_int_cor(double *x, double *p)
{
  double W = x[0]+M_PSI2S;
  double SPREAD=p[0];
  double QED = p[1];
  double EFF = p[2];
  double CMIN = p[3];
  double CMAX = p[4];
  return ee_interference_spreaded(W,SPREAD, EFF, CMIN, CMAX)/ee_cs_qed_spreaded(W,SPREAD,QED);
};


//subintegral 
inline double sigma_spread_sub(double *x,  double *P)
{
	double W=P[0];
	double S=P[1];
	double p[4];
	p[0] = P[2];
	p[1] = P[3];
	p[2] = P[4];
	p[3] = P[5];
	double res = sigma(x, p)*spread(W-x[0], S);
	return res;
};


//the bhabha with the interference and energy spread
inline double sigma_spread(double *x,  double *par)
{
	double p[6];
	p[0] = x[0]; //beam energy
	p[1] = par[0]; //spread
	p[2] = par[1]; //QED
	p[3] = par[2]; //Interf
	p[4] = par[3]; //Res
	p[5] = par[4]; //Gamma
	TF1  f("tmpbbfun", &sigma_spread_sub, -20, 20, 6);
	f.SetParameters(p);
	double res = f.Integral(p[0]-p[1]*5, p[0]+p[1]*5);
	//clog << "W="<<p[0]<< ", sig="<<p[1] << ",  qed="<<p[2]<<", int="<< p[3] << ", res=" << p[4]
	//	<< ", gam=" << p[5] << ",  sigma_spread=" << res << endl;
	return res;
}

//double BBIntCor(double W)
//{
//	//this number are taken from MC draw_bhabha
//	double SPREAD = 1.58; //MeV
//	double QED = 125.057; //nb
//	double INT = 10.7287; //nb
//	double RES = -4.34519; //nb
//	double GAMMA = 0.304; //MeV
//  double cor=0; //result of calculation
//	TF1 * sfun = new TF1("fbbcor_tmp",&sigma_spread,-10, 10, 5 );
//	double par[5];
//	par[0]=SPREAD; //Beam spread
//	par[1] = 0;
//	par[2]=INT/QED;//INT
//	par[3]=RES/QED;//RES
//	par[4]=GAMMA;//GAMMA
//	sfun->SetParameters(par);
//  cor = 1./(1.+sfun->Eval(W-MPDG));
//  delete sfun;
//	return cor;
//};

inline double BBIntCor(double W,double SPREAD=1.58)
{
  double QED = 124.7;
  double EFF = 0.5545;
  double CMIN = 0.86;
  double CMAX = 0.93;
  double cor = ee_interference_spreaded(W,SPREAD, EFF, CMIN, CMAX)/ee_cs_qed_spreaded(W,SPREAD,QED);
  return 1./(1.+cor);
}

inline double MhadrCor(double W)
{
  double k=5.164e-4;
  return  1.+ k*(W-M_PSI2S);
}

#endif
