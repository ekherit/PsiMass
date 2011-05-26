/*
 * =====================================================================================
 *
 *       Filename:  selection.h
 *
 *    Description:  Selection criteria are here
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
#ifndef IBN_SELECTION_H
#define IBN_SELECTION_H

#include <TCut.h>
void set_selection(int selection_version, TCut & mh_cut, TCut & ee_cut , TCut & gg_cut)
{
  TCut mh_base_cut; //base cut for signal
  TCut mh_strict_cut; //strict cut for signal
  TCut ee_base_cut; //base cut for bhabha
  TCut ee_ext_cut;
  TCut ee_theta_cut;
  TCut mh_theta_cut;
  TCut gg_base_cut;
  TCut gg_theta_cut;
  TCut hp_cut;
  TCut rv_cut;
  TCut ip_cut[3]; //interaction point cut
  TCut mdcEcut;
  TCut ggEcut;
  switch(selection_version)
  {
    case 3:
      mh_base_cut = "nemc>2  && S>0.05 && Eemc<2.5 && Emdc<4";
      mh_strict_cut = "nemc>3  && S>0.05 && Eemc<2.5 && Emdc<4";
      ee_base_cut = "nemc==2 && S<0.05 && Emdc<5 && Eemc>2.5";
      ee_ext_cut = "(nemc==2 || nemc==3) && S<0.05 && Emdc<5 && Eemc>2.5";
      ee_theta_cut  = "Sum$(sin(theta)<0.55)==mdc.ntrack";
      mh_theta_cut  = "Sum$(sin(theta)>0.45)==mdc.ntrack";
      gg_base_cut = "Etotal > 3.3 && Etotal < 4  && sqrt((Sum$(x)-2)**2 + Sum$(y)**2)<4 && abs(Sum$(z))<9";
      gg_theta_cut  = "Sum$(sin(theta)>0.45)==2";
      hp_cut = "Sum$(abs(hpz)<3)==2 && Sum$(hpr<0.25)==2";
      rv_cut = "Sum$(rvxy[hpidx]<0.5)==2 && Sum$(abs(rvz[hpidx])<5)==2";
      /* this is strict cut */
      //mh_cut = mh_strict_cut  && rv_cut && mh_theta_cut && "pt100";
      /* bha bha ext cut */
      //ee_cut = ee_ext_cut  && rv_cut && ee_theta_cut;
      /* very strict cut */
      //mh_cut = mh_strict_cut && rv_cut && mh_theta_cut && "emc.ntrack==0";
      mdcEcut = "Sum$(E>0.05)==nemc";
      ggEcut = "Sum$(E>0.05)==2";

      //mh_cut = mh_base_cut && rv_cut;
      //ee_cut = ee_base_cut && rv_cut &&  ee_theta_cut;
      //gg_cut = gg_base_cut  && gg_theta_cut && ggEcut;

      /* new test cut for E cut */
      //mh_cut = "Sum$(E>0.02)>=3 && S>=0.06 && Eemc<2.5 && Emdc<5"  && rv_cut;
      //ee_cut = "Sum$(E>0.02)==2 && S<=0.05 && Emdc<5 && Eemc>2.5"  && rv_cut && ee_theta_cut;
      gg_cut = gg_base_cut  && gg_theta_cut && "Sum$(E>0.02)==2";

      mh_cut = "ngt > 2 &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut && "nemc==ntrack";
      //mh_cut = "nemc > 2 &&  S>=0.06 && Eemc<2.5 && Emdc<5" && rv_cut && mdcEcut;
      //mh_cut = ("ngt > 2 &&  S>=0.06" || ( "ngt==2 && S>=0.06" && rv_cut && mh_theta_cut)) && "ngt_Eemc<2.5 && Emdc<5";
      ee_cut = "ngt == 2 &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut;
      //ee_cut = "nemc == 2 &&  S<=0.05 && Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut &&mdcEcut;

      //strict cut
      //mh_cut = "ngt >= 4  &&  S>=0.06 && ngt_Eemc<2.5 && Emdc<5" && rv_cut && mh_theta_cut && "pt100";
      //ee_cut = "ngt == 2  &&  S<=0.05 && ngt_Eemc>2.5 && Emdc<5" && ee_theta_cut && rv_cut && "emc.ntrack==0";
      //gg_cut = gg_base_cut  && gg_theta_cut;
      break;
    case 4:
      {
        mh_base_cut = "mdc.ntrack>2 && S>0.06 && Eemc<2.5 && Emdc<4 && pt50";
        mh_strict_cut = "mdc.ntrack>3 && S>0.06 && Eemc<2.5 && Emdc<4 && pt100";
        mh_theta_cut  = "Sum$(sin(theta)>0.45)==mdc.ntrack";
        ee_base_cut = "mdc.ntrack==2 && S<0.05 && Eemc>2.5 && Emdc<5";
        ee_theta_cut = "Sum$(sin(theta)<0.45)==mdc.ntrack";
        gg_base_cut = "Etotal > 3.3 && Etotal < 4  && sqrt((Sum$(x)-2)**2 + Sum$(y)**2)<4 && abs(Sum$(z))<9";
        gg_theta_cut = "Sum$(theta)>0.45";
        TCut gg_cos_cut = "acos(cos)>3.1";
        //TCut tof_cut = "Sum$(tof>0)==tof.tof.ntrack && Sum$(tof<10)==tof.tof.ntrack";
        //TCut tof_cut = "Sum$(t0>540 && t0<660)==Length$(t0)";
        TCut tof_cut = "Sum$(tof>1&&tof<6)==Length$(tof) && Sum$(t0>540 && t0<640)==Length$(t0)";

        for(int i=0;i<3;i++)
        {
          char buf[1024];
          double IP_R=0.5;//cm
          double IP_Z=5;//cm
          sprintf(buf,"rvxy[%d]<%f&&abs(rvz[%d])<%f", i, IP_R, i, IP_Z);
          ip_cut[i]=TCut(buf);
        }
        //mh_cut = mh_base_cut  && ip_cut[0] && ip_cut[1] && ip_cut[2]&& mh_theta_cut;
        mh_cut = mh_base_cut  && "Sum$(mdc.rvxy<0.5)==mdc.ntrack && Sum$(abs(mdc.rvz)<5)==mdc.ntrack" && mh_theta_cut;
        ee_cut = ee_base_cut  && ip_cut[0] && ip_cut[1] && ee_theta_cut && "q[0]*q[1]<0";
        gg_cut = gg_base_cut  && gg_theta_cut &&gg_cos_cut;
      }
      break;

    case 6:
      {
        TCut good_track = "Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack && Sum$(mdc.rvxy<1)==mdc.ntrack";
        mh_base_cut  =   "mdc.ntrack>3" && good_track;
        mh_cut = mh_base_cut;

        TCut ee_acol =   "abs(atheta)<0.03 && -0.06 < aphi && aphi<0.01";
        TCut ee_endcup = "abs(cos(mdc.theta[0]))>0.86 && abs(cos(mdc.theta[1]))>0.86";
        TCut ee_E("mdc.E[0]/Eb>0.8 && mdc.E[1]/Eb>0.8 && Emdc<5");
        ee_base_cut = "mdc.ntrack>1 && mdc.ntrack<4" && good_track && ee_endcup && ee_acol && ee_E;
        ee_cut = ee_base_cut  && "q[0]*q[1]<0";

        TCut gg_acol = "abs(atheta) < 0.05 &&  aphi>-0.08 && aphi<0.04";
        TCut gg_barel = "Sum$(cos(theta)<0.8)";
        TCut gg_E("E[0]/Eb && E[1]/Eb>0.8");
        gg_base_cut = gg_acol && gg_barel && gg_E;
        gg_cut = gg_base_cut;
      }
      break;
    case 7:
      {
        TCut good_track = "Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack && Sum$(mdc.rvxy<1)==mdc.ntrack";
        TCut mh_p = "Sum$(p<2.5)==mdc.ntrack";
        TCut mh_2good = "mdc.rvxy[0]<1 && mdc.rvxy[1]<1&& Sum$(abs(cos(mdc.theta))<0.93)==mdc.ntrack";
        TCut mh_2good2 = "mdc.rvxy[0]<1 && mdc.rvxy[1]<1&& abs(cos(mdc.theta[0]))<0.93 &&abs(cos(mdc.theta[1]))<0.93";
        mh_base_cut  =   "mdc.ntrack>3" && good_track && "S>0.06"&& mh_p;
        //mh_base_cut  =   "mdc.ntrack>2" && good_track; //week cut
        //mh_base_cut  =   "mdc.ntrack>4" && good_track && "S>0.06" && "Sum$(cos(theta)<1.8&&pt>0.05)==mdc.ntrack"; //strong cut
        //mh_base_cut = "mdc.ntrack>2" && "S>0.06" && mh_2good2;
        mh_cut = mh_base_cut;

        TCut ee_acol =   "abs(atheta)<0.03 && -0.06 < aphi && aphi<0.01";
        TCut ee_endcup = "abs(cos(mdc.theta[0]))>0.86 && abs(cos(mdc.theta[1]))>0.86";
        TCut ee_barrel = "abs(cos(mdc.theta[0]))<0.8 && abs(cos(mdc.theta[1]))<0.8";
        TCut ee_barrel2 = "abs(cos(mdc.theta[0]))>0.5 && abs(cos(mdc.theta[1]))>0.5";
        TCut ee_E("mdc.E[0]/Eb>0.8 && mdc.E[1]/Eb>0.8 && mdc.E[0]/Eb<1.2 && mdc.E[1]/Eb<1.2");
        TCut ee_p =  "mdc.p[0]<2.5 && mdc.p[1]<2.5 && mdc.p[0]/Eb>0.9 && mdc.p[1]/Eb>0.9";
        //ee_base_cut = "mdc.ntrack>1 && mdc.ntrack<4" && good_track && ee_barrel && ee_barrel2 && ee_acol && ee_E && ee_p;
        ee_base_cut = "mdc.ntrack>1 && mdc.ntrack<4" && good_track && ee_endcup && ee_acol && ee_E && ee_p;
        ee_cut = ee_base_cut  && "q[0]*q[1]<0";
        
        cout << "Bhabha cut selection: " << endl;
        cout << ee_cut << endl;

        TCut gg_n = "ngct==0 && ngt>1";
        TCut gg_acol = "abs(atheta) < 0.05 &&  aphi>-0.06 && aphi<0.02";
        TCut gg_barel = "abs(cos(theta[0]))<0.8 && abs(cos(theta[1]))<0.8";
        TCut gg_E = "E[0]/Eb>0.8 && E[1]/Eb>0.8 && E[0]/Eb<1.2 && E[1]/Eb<1.2";
        gg_base_cut = gg_acol && gg_barel && gg_E && gg_n;
        gg_cut = gg_base_cut;
      }
      break;
  }
}
#endif
