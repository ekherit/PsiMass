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
const double MPDG=3686.09;
const double PI = 3.1415926535897;
const double ALPHA = 1./136.037;
const double GAMMAEE_PSI2S = 0.00235; //MeV
const double GAMMA_PSI2S = 0.304; // MeV
const double M_PSI2S  = 3686.09; //MeV
const double ME = 0.510998910; //MeV
const double HC2 = 0.389379304*1e6*1e3; //MeV^2*nbarn

#include <TF1.h>
double sigma(double *x,  double *p)
{
	double W = MPDG+x[0];
	if(W<=0) return 1e100;
	double G = p[3]; //width
  double den=sqrt(sq(W-MPDG)+sq(G/2.));
	//double a = atan(-G/2./(W-MPDG));
	double a = asin(G/(2*den));
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
double spread(double dW,   double S)
{
	return 1./(sqrt(2*PI)*S)*exp(-0.5*sq(dW/S));
}


//subintegral 
double sigma_spread_sub(double *x,  double *P)
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
double sigma_spread(double *x,  double *par)
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


double BBIntCor(double W)
{
	//this number are taken from MC draw_bhabha
	double SPREAD = 1.6; //MeV
	double QED = 124.505; //nb
	double INT = 0.89923; //nb
	double RES = 3.69557; //nb
	double GAMMA = 2.70303; //MeV
  double cor=0; //result of calculation
	TF1 * sfun = new TF1("fbbcor_tmp",&sigma_spread,-10, 10, 5 );
	double par[5];
	par[0]=SPREAD; //Beam spread
	par[1] = 0;
	par[2]=INT/QED;//INT
	par[3]=RES/QED;//RES
	par[4]=GAMMA;//GAMMA
	sfun->SetParameters(par);
  cor = 1./(1.+sfun->Eval(W-MPDG));
  delete sfun;
	return cor;
};

inline double ee_interference(double W, double theta);

double ee_interference_vs_W(double *x, double *p)
{
    return  ee_interference(x[0], p[0]);
}

double ee_interference_vs_theta(double *x, double *p)
{
    return  ee_interference(p[0], x[0]);
}

double ee_interference(double W, double theta_min, double theta_max);

#endif
