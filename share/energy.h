/*
 * =====================================================================================
 *
 *       Filename:  energy.h
 *
 *    Description:  Energy manipulation.
 *
 *        Version:  1.0
 *        Created:  23.01.2012 17:07:19
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
#ifndef IBN_BEPC_ENERGY_H
#define IBN_BEPC_ENERGY_H

#include "utils.h" 
const double Me=0.510998918; //Electron mass, MeV
const double BEPC_ALPHA=0.022; //BEPC crossing angle

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

double dcm_energy(double Ee, double dEe, double Ep, double dEp)
{
    return sqrt(sq(dW_dE1(Ee,Ep)*dEe) + sq(dW_dE1(Ep,Ee)*dEp));
}
#endif
