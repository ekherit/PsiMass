/*
 * =====================================================================================
 *
 *       Filename:  RunInfo.h
 *
 *    Description:  Run information
 *
 *        Version:  1.0
 *        Created:  20.01.2012 10:58:42
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#ifndef IBN_BES_PSIP_MASS_H
#define IBN_BES_PSIP_MASS_H
#include "ibn/valer.h"
#include <time.h>
#include <vector>
struct EMS_t
{
    time_t begin_time;
    time_t end_time;
    double E; //beam energy, MeV
    double dE; 
    double S;//energy spread, MeV
    double dS;
};

struct RunInfo_t
{
  unsigned run;
  unsigned NBee;
  unsigned NEee;
  unsigned Nuu;
  unsigned Ngg;
  unsigned Nh;
  unsigned Nmh; //number of multihadron
  double SNR;
  double lum;
  double xsec;
  time_t begin_time;
  time_t end_time;
  double BEPC_Ee;
  double BEPC_Ep;
  double Ie1, Ie2;
  double Ip1, Ip2;
  double W;
  double dW;
  EMS_t e,p; //energy and error for electron (e) and positron p
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
    begin_time =0;
    end_time = 0;
    Nmh=0;
    e.E=0;
    p.E=0;
    e.dE=0;
    p.dE=0;
    e.S=0;
    p.S=0;
    e.dS=0;
    p.dS=0;
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

extern void read_run_info(const char * filename, RunInfo_t * RI, unsigned & RIsize);
extern void read_run_info(const char * filename, std::vector<RunInfo_t> & RI);

#endif
