
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
  double I = - 1.5*ALPHA*gee*( 1.+sq(cos(theta)) - sq(1.+cos(theta))/(1.-cos(theta))*F.real());
  double res = 1./sq(M_PSI2S)*(R + I);
  //cout << "Res = " << R << ", Inte = " << I <<  " res = " << res << endl;;
  return HC2*res; //result should be in nano barn
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


double ee_interference(double W, double theta_min, double theta_max)
{
  return ibn::dgaus(eeint_sub(W), theta_min, theta_max,1e-12);
}
