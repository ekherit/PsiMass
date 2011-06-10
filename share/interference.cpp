
/*
 * =====================================================================================
 *
 *       Filename:  interference.cpp
 *
 *    Description:  Calculate the interference in ee->ee process betwee resonance and continuum
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

#include "ibn/integral.h"

#include "interference.h"
#include <cmath>



#include <complex>



inline double ee_interference(double W, double theta)
{
  const double LWME = log(W/ME);
  double b  = 4*ALPHA/M_PI*(LWME-0.5);
  const double PI2=M_PI*M_PI;
  const double d = 0.75*b+ALPHA/M_PI*(PI2/3. - 0.5)+b*b*(37./96.-PI2/12.-LWME/36.);
  const double g = GAMMA_PSI2S/M_PSI2S;
  const double gee = GAMMAEE_PSI2S/M_PSI2S;
  double x = W/M_PSI2S;
  complex <double> den(2.*(-x+1.), - g);
  complex <double> F = M_PI*b/sin(M_PI*b)*pow(den, b-1.0);
  double R = 9./4.*gee*gee/g*(1.+sq(cos(theta)))*(1.+d)*F.imag();
  double I = - 1.5*ALPHA*gee*( 1.+sq(cos(theta)) - sq(1.+cos(theta))/(1.-cos(theta)))*F.real();
  double result = 1./sq(M_PSI2S)*(R + I);
  //cout << "Res = " << R << ", Inte = " << I <<  " res = " << res << endl;;
  return HC2*result; //result should be in nano barn
}

inline double ee_interference(double W, double cos_min, double cos_max)
{
  const double LWME = log(W/ME);
  double b  = 4*ALPHA/M_PI*(LWME-0.5);
  const double PI2=M_PI*M_PI;
  const double d = 0.75*b+ALPHA/M_PI*(PI2/3. - 0.5)+b*b*(37./96.-PI2/12.-LWME/36.);
  const double g = GAMMA_PSI2S/M_PSI2S;
  const double gee = GAMMAEE_PSI2S/M_PSI2S;
  double x = W/M_PSI2S;
  complex <double> den(2.*(-x+1.), - g);
  complex <double> F = M_PI*b/sin(M_PI*b)*pow(den, b-1.0);
  double R_angular_integral = 2.*M_PI*((cb(cos_max)-cb(cos_min))/3.+cos_max-cos_min);
  double R = 9./4.*gee*gee/g*(1.+d)*F.imag()*R_angular_integral;
  double I_angular_integral = 2.*M_PI*( 4.*log((1.-cos_max)/(1-cos_min)) + (cb(cos_max)-cb(cos_min))/3.+cos_max-cos_min +
      0.5*(sq(cos_max)-sq(cos_min)+6.*(cos_max-cos_min))); 
  double I = - 1.5*ALPHA*gee*I_angular_integral*F.real();
  double result = HC2/sq(M_PSI2S)*(R + I); //result should be in nano barn
  //cout << "W=" << W << " Res = " << R << ", Inte = " << I <<  " res = " << result << " Rang = " << R_angular_integral/2./M_PI << " Iang="<< I_angular_integral/2./M_PI << endl;
  //if(isnan(result)) return 1e300;
  return result; //result should be in nano barn
}


class eeint_sub
{
    double W;
    public:
    eeint_sub(double w) { W=w; }
    double operator()(double theta)
    {
      return ee_interference(W,theta)*sin(theta)*2.*M_PI;
    } 
};


class ee_cs_spreaded_sub
{
  const double W; //c.m. beam energy, MeV
  const double SPREAD; //beam energy spread in c.m.,MeV
  const double QED; //QED part in nb.
  const double EFF; //registration efficiency for interference
  const double CMIN; //minimum cosinus
  const double CMAX;//maximum cosinus
  public:
  ee_cs_spreaded_sub(double w, double s, double qed, double eff, double cmin, double cmax):
      W(w), SPREAD(s),QED(qed), EFF(eff), CMIN(cmin), CMAX(cmax){}
  double operator()(double w)
  {
    double C = QED*sq(M_PSI2S/w);
    double I = (EFF == 0) ? 0  : EFF*2.0*ee_interference(w, CMIN, CMAX); 
    return  (I+C)*spread(w-W,SPREAD);
  }
};

//class ee_interference_spreaded_sub
//{
//  const double W; //c.m. beam energy, MeV
//  const double SPREAD; //beam energy spread in c.m.,MeV
//  const double EFF; //registration efficiency for interference
//  const double CMIN; //minimum cosinus
//  const double CMAX;//maximum cosinus
//  const double QED; //QED part in nb.
//  public:
//  ee_interference_spreaded_sub(double w, double s, double eff, double cmin, double cmax, double qed=0):
//      W(w), SPREAD(s), EFF(eff), CMIN(cmin), CMAX(cmax), QED(qed){}
//  double operator()(double w)
//  {
//    double C = QED*sq(M_PSI2S/w);
//    double I = EFF*2.0*ee_interference(w, CMIN, CMAX); 
//    return  (I+C)*spread(w-W,SPREAD);
//  }
//};

double ee_interference_spreaded(double W, double spread, double int_eff, double cmin, double cmax)
{
  return ibn::dgaus(ee_cs_spreaded_sub(W,spread,0, int_eff,cmin,cmax),W-spread*5,W+spread*5,1e-12);
}

double ee_cs_qed_spreaded(double W, double spread, double qed)
{
  return ibn::dgaus(ee_cs_spreaded_sub(W,spread, qed, 0,0,0),W-spread*5,W+spread*5,1e-12);
}

