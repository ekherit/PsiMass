/*
 * =====================================================================================
 *
 *       Filename:  TauEMU.h
 *
 *    Description:  Selection of tau decey into e mu
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
#define MAX_TRACK_NUMBER 9
#define MIN_CHARGED_TRACK 3

#include <vector>
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
	double prop_delta_x; //interaction  point
	double prop_delta_y;
	double prop_delta_z;
	long int event_proceed;
	long int event_write;
	NTuple::Tuple * main_tuple;//main tuple 
	NTuple::Item<double> Etotal; //Total energy deposition
	NTuple::Item<long> nchtrk; //number of charged tracks
	NTuple::Item<long> nneutrk;//number of neutral tracks
	NTuple::Item<long> ntrk; //total number of tracks.
	NTuple::Item<long> niptrk; //tracks from interection point

	//charged track information
	NTuple::Tuple * chtr_tuple; //charged grack tuple
  NTuple::Item<long> tr_idx;
  NTuple::Array<double> m_E;
  NTuple::Array<double> m_pt;
  NTuple::Array<double> m_M;
  NTuple::Array<double> m_q;
  NTuple::Array<double> m_x, m_y, m_z;
  NTuple::Array<double> m_ismu;
	//struct track_t
	//{
	//	NTuple::Item<double> pt; //transverse momentum
	//	NTuple::Item<double> E;  //Energy
	//	NTuple::Item<double> M; //Invariant mass
	//	NTuple::Item<double> q;//charge
	//	NTuple::Item<double> x, y, z; //vertex coordinate
	//	NTuple::Item<double> ismu;  //is it a muon track
	//};
	//std::vector <track_t> chtr;

	//de_dx infrmation temporary remove
	//NTuple::Tuple * dedx_tuple;
  //NTuple::Item<double> m_ptrk[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_chie[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_chimu[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_chipi[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_chik[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_chip[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_probPH[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_normPH[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_ghit[MAX_TRACK_NUMBER];
  //NTuple::Item<double> m_thit[MAX_TRACK_NUMBER];
	void InitData(void);
};

#endif
