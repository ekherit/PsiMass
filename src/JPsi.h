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

#include <TMatrixD.h>
#include <vector>
#include <algorithm>
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
	NTuple::Item<double> S1;
	NTuple::Item<double> S2;
	NTuple::Item<double> S3;
	NTuple::Item<double> m_S;
	NTuple::Item<long> m_Signal;

  NTuple::Item<long> tr_idx;
	//charged track information
	NTuple::Tuple * chtr_tuple; //charged grack tuple
  NTuple::Array<double> m_E;
  NTuple::Array<double> m_pt;
  NTuple::Array<double> m_M;
  NTuple::Array<double> m_q;
  NTuple::Array<double> m_x, m_y, m_z;
  NTuple::Array<double> m_ismu;

	void InitData(void);
	TMatrixD S; //sphericity tensor

	/*  BhaBha specific section section */
};

#endif
