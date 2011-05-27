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
const double MPDG=3686.09;
const double ALPHA=1./136.037;
const double PI = 3.1415926535897;
const double ME = 0.51099; //Mass of the electron,  MeV

double sigma(double *x,  double *p)
{
	double W = MPDG+x[0];
	double G = p[3]; //width
	double a = atan(-G/2./(W-MPDG));
	double beta = 4*ALPHA/PI*(log(W/ME)-0.5);
	//resonance beheavioure
	double A=PI*beta/sin(PI*beta)*G/2./sqrt(sq(W-MPDG)+sq(G/2.));
	//interference contribution max
	double C = p[0]*sq(MPDG/W); //continuum
	double I = p[1]*A*cos(a*(1-beta)); 
	double R = p[2]*A*sin(a*(1-beta)); 
	return C+I+R;
}

//energy spread
double spread(double dW,   double S)
{
	const double PI=3.1415926535;
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
	TF1 * f = new TF1("tmp", &sigma_spread_sub, -20, 20, 6);
	f->SetParameters(p);
	double res = f->Integral(p[0]-p[1]*5, p[0]+p[1]*5);
	return res;
}


double BBIntCor(double W)
{
	//this number are taken from MC draw_bhabha
	double SPREAD = 1.6; //MeV
	double QED = 124.479; //nb
	double INT = 1.06325; //nb
	double RES = 3.82139; //nb
	double GAMMA = 2.66043; //MeV
	TF1 * f = new TF1("fbbcor",  &sigma_spread,  5);
	f->SetParameter(0,SPREAD);
	f->SetParameter(1, 0);//QED is zero because we are looking for the correction.
	f->SetParameter(2, INT/QED);
	f->SetParameter(3, RES/QED);
	f->SetParameter(4, GAMMA);
	double cor =  1./(1+f->Eval(W-MPDG));
	cout << W << " " << cor << endl;
	delete f;
	return cor;
};


#endif
