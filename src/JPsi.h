/*
 * =====================================================================================
 *
 *       Filename:  JPsi.h
 *
 *    Description:  Multihadron event selection for j/psi and psi prime resonance.
 *
 *        Version:  1.0
 *        Created:  04/27/2010 02:47:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#ifndef IBN_TAUEMU_H
#define IBN_TAUEMU_H
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"

#include <TMatrixD.h>
#include <vector>
#include <algorithm>

//#include "EventModel/EventModel.h"
//#include "EventModel/Event.h"
#include "EventModel/EventHeader.h"

class JPsi : public Algorithm 
{
	public:
  JPsi(const std::string& name, ISvcLocator* pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();  

	private:
  int prop_check_dedx;
	int USE_IPCUT; //use interection point cut
	int IPTRACKS; //tracks number from interection point
	int MIN_CHARGED_TRACKS; //minimum charged tracks in selection
	int MAX_TRACK_NUMBER; //minimum charged tracks in selection
	double DELTA_X, DELTA_Y, DELTA_Z; //interection point cut
	long int event_proceed;
	long int event_write;
	NTuple::Tuple * main_tuple;//main tuple 
	NTuple::Item<long> m_time; //time when events is writed (unixtime)
	NTuple::Item<double> Etotal; //Total energy deposition
	NTuple::Item<long> nchtrk; //number of charged tracks
	NTuple::Item<long> nneutrk;//number of neutral tracks
	NTuple::Item<long> ntrk; //total number of tracks.
	NTuple::Item<long> niptrk; //tracks from interection point
	NTuple::Item<long> m_nmdcemc; //good tracks with mdc and emc
	NTuple::Item<double> S1;
	NTuple::Item<double> S2;
	NTuple::Item<double> S3;
	NTuple::Item<double> m_S;
	NTuple::Item<double> m_cos_high_p;

  NTuple::Item<long> tr_idx;
	//charged track information
  NTuple::Array<double> m_E;
  NTuple::Array<double> m_p;
  NTuple::Array<double> m_px,m_py,m_pz;
	NTuple::Array<double> m_theta;
	NTuple::Array<double> m_phi;
  NTuple::Array<double> m_pt;
  NTuple::Array<double> m_M;
  NTuple::Array<double> m_q;
  NTuple::Array<double> m_x, m_y, m_z;
  NTuple::Array<double> m_X, m_Y, m_Z;
  NTuple::Array<double> m_ismu;


  /* ElecroMagnetic Calorimeter Information */
  struct EMC_t
  {
    NTuple::Item<long>    ntrack; //number of neutral tracks.
    NTuple::Array<long>   status; //status status=1: single seed cluster; status=2: splitted from multi-seeds cluster.
    NTuple::Array<long>   ncrstl; //Number of crystals in the shower
    NTuple::Array<long>   cellId; //Central crystal’s identifier
    NTuple::Array<long>   module; //module=0: east endcap;  module=1: barrel;  module=2: west endcap.
    NTuple::Array<double> x, y, z; //coordinates of claster
    NTuple::Array<double> theta, phi;  //angles
    NTuple::Array<double> E; // energy deposition
    NTuple::Array<double> dE; // energy deposition error
    NTuple::Item<double>  Etotal;
  };
	
  NTuple::Tuple * EMC_tuple;
  EMS_t emc;



	NTuple::Tuple * dedx_tuple;
  NTuple::Item<long> trdedx_idx;
	NTuple::Array<double> m_chie; //chi2_dEdx for electron
	NTuple::Array<double> m_chimu;//chi2_dEdx for muon
	NTuple::Array<double> m_chipi;//chi2_dEdx pion
	NTuple::Array<double> m_chik; //chi2_dEdx kaon
	NTuple::Array<double> m_chip; //chi2_dEdx proton
	NTuple::Array<double> m_ghit; //number of good de/dx hits (excluding overflow)
	NTuple::Array<double> m_thit; //nuber of total de/dx including overflow
  NTuple::Array<double> m_probPH; //most probable pulse height from trucated mean 
  NTuple::Array<double> m_normPH; //normalized pulse height
	void InitData(void);
	TMatrixD S; //sphericity tensor

	//gamma-gamma annihilation selection
	long gg_event_writed;
	NTuple::Tuple * gg_tuple;
  NTuple::Item<long> gg_nntrk;
	NTuple::Item<double> gg_cos; //aclolinearity
	NTuple::Item<double> gg_Etotal; //aclolinearity
	NTuple::Array<double> gg_x, gg_y, gg_z; //coordinat of claster
	NTuple::Array<double> gg_theta, gg_phi;//angles
	NTuple::Array<double> gg_E, gg_dE; //energy and error
	NTuple::Array<long> gg_n; //number of clasters
	NTuple::Array<long> gg_module; //module=0: east endcap;  module=1: barrel;  module=2: west endcap.
};

#endif
