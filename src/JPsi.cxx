/*
 * =====================================================================================
 *
 *       Filename:  JPsi.cxx
 *
 *    Description:  Multihadron event selection for j/psi and psi prime resonance.
 *
 *        Version:  1.0
 *        Created:  04/27/2010 02:50:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "VertexFit/IVertexDbSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/ISvcLocator.h"

#include "EventModel/EventModel.h"
#include "EventModel/Event.h"

#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "DstEvent/TofHitStatus.h"
#include "EventModel/EventHeader.h"

#include "TMath.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/NTuple.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/TwoVector.h"
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
   typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "JPsi.h"

#include <TMatrixDEigen.h>

#include "VertexFit/KinematicFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/Helix.h"
#include "ParticleID/ParticleID.h"

JPsi::JPsi(const std::string& name, ISvcLocator* pSvcLocator) :
  Algorithm(name, pSvcLocator)
{
  declareProperty("CheckDedx", prop_check_dedx = 1);
  declareProperty("DELTA_X", DELTA_X = 1.0); //cm?
  declareProperty("DELTA_Y", DELTA_Y = 1.0); //cm?
  declareProperty("DELTA_Z", DELTA_Z = 10.0); //cm?
  declareProperty("USE_IPCUT", USE_IPCUT=1); //to use interection point cut.
  declareProperty("IPTRACKS", IPTRACKS=2); //number of tracks from interection point
  declareProperty("MIN_CHARGED_TRACKS", MIN_CHARGED_TRACKS=2); //minimum number of charged tracks in selection
  declareProperty("MAX_TRACK_NUMBER", MAX_TRACK_NUMBER=20); //maximum number of charged tracks
	S.ResizeTo(3, 3);
}


StatusCode JPsi::initialize(void)
{
  MsgStream log(msgSvc(), name());

  log << MSG::INFO << "in initialize()" << endmsg;
	event_proceed=0;
	event_write = 0;
	gg_event_writed=0;
  
  StatusCode status;
	NTuplePtr my_nt(ntupleSvc(), "FILE1/mhadr");
	if(my_nt) main_tuple=my_nt;
	else
	{
		main_tuple = ntupleSvc()->book("FILE1/mhadr", CLID_ColumnWiseTuple, "Multihadron tree plus bhabha");
		if(main_tuple)
		{
			//common
			status=main_tuple->addItem("t", m_time);
			status=main_tuple->addItem("Etotal", Etotal);
			status=main_tuple->addItem("nchtrk", nchtrk);
			status=main_tuple->addItem("nneutrk", nneutrk);
			status=main_tuple->addItem("ntrk", ntrk);
			status=main_tuple->addItem("niptrk", niptrk);
			/*  sphericity part */
			status=main_tuple->addItem("S1", S1);
			status=main_tuple->addItem("S2", S2);
			status=main_tuple->addItem("S3", S3);
			status=main_tuple->addItem("S", m_S);
			status=main_tuple->addItem("signal", m_Signal);
			//charged tracks
			status=main_tuple->addItem("coshp", m_cos_high_p);
      status = main_tuple->addItem ("nchtr", tr_idx, 0, MAX_TRACK_NUMBER);
      status = main_tuple->addIndexedItem ("E", tr_idx, m_E );
      status = main_tuple->addIndexedItem ("p", tr_idx, m_p );
      status = main_tuple->addIndexedItem ("px", tr_idx, m_px );
      status = main_tuple->addIndexedItem ("py", tr_idx, m_py );
      status = main_tuple->addIndexedItem ("pz", tr_idx, m_pz );
      status = main_tuple->addIndexedItem ("pt", tr_idx, m_pt );
      status = main_tuple->addIndexedItem ("theta", tr_idx, m_theta );
      status = main_tuple->addIndexedItem ("phi", tr_idx, m_phi );
      status = main_tuple->addIndexedItem ("M", tr_idx, m_M );
      status = main_tuple->addIndexedItem ("q", tr_idx, m_q );
      status = main_tuple->addIndexedItem ("x", tr_idx, m_x );
      status = main_tuple->addIndexedItem ("y", tr_idx, m_y );
      status = main_tuple->addIndexedItem ("z", tr_idx, m_z );
      status = main_tuple->addIndexedItem ("ismu", tr_idx, m_ismu );
		}
		else
		{
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(main_tuple) << endmsg;
      return StatusCode::FAILURE;
		}
	}
	NTuplePtr nt2(ntupleSvc(), "FILE1/dedx");
	if(nt2) dedx_tuple=nt2;
	else
	{
		dedx_tuple = ntupleSvc()->book("FILE1/dedx", CLID_ColumnWiseTuple, "dedx information");
		if(dedx_tuple)
		{
      status = dedx_tuple->addItem ("nchtr", trdedx_idx, 0, MAX_TRACK_NUMBER);
      status = dedx_tuple->addIndexedItem ("chie", trdedx_idx, m_chie );
      status = dedx_tuple->addIndexedItem ("chimu",trdedx_idx, m_chimu );
      status = dedx_tuple->addIndexedItem ("chipi", trdedx_idx, m_chipi );
      status = dedx_tuple->addIndexedItem ("chik", trdedx_idx, m_chik );
      status = dedx_tuple->addIndexedItem ("chip", trdedx_idx, m_chip );
      status = dedx_tuple->addIndexedItem ("ghit", trdedx_idx, m_ghit );
      status = dedx_tuple->addIndexedItem ("thit", trdedx_idx, m_thit );
      status = dedx_tuple->addIndexedItem ("probPH", trdedx_idx, m_probPH );
      status = dedx_tuple->addIndexedItem ("normPH", trdedx_idx, m_normPH );
		}
		else
		{
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(dedx_tuple) << endmsg;
      return StatusCode::FAILURE;
		}
	}
	NTuplePtr nt_gg(ntupleSvc(), "FILE1/gg");
	if(nt_gg) gg_tuple=nt_gg;
	else
	{
		gg_tuple = ntupleSvc()->book("FILE1/gg", CLID_ColumnWiseTuple, "gamma-gamma annihilation");
		if(gg_tuple)
		{
      status = gg_tuple->addItem ("nneu", gg_nntrk, 0, 2);
      status = gg_tuple->addItem ("cos", gg_cos);
      status = gg_tuple->addItem ("Etotal", gg_Etotal);
      status = gg_tuple->addIndexedItem ("x", gg_nntrk, gg_x );
      status = gg_tuple->addIndexedItem ("y", gg_nntrk, gg_y );
      status = gg_tuple->addIndexedItem ("z", gg_nntrk, gg_z );
      status = gg_tuple->addIndexedItem ("theta", gg_nntrk, gg_theta );
      status = gg_tuple->addIndexedItem ("phi", gg_nntrk, gg_phi );
      status = gg_tuple->addIndexedItem ("E", gg_nntrk, gg_E );
      status = gg_tuple->addIndexedItem ("dE", gg_nntrk, gg_dE );
      status = gg_tuple->addIndexedItem ("model", gg_nntrk, gg_module );
      status = gg_tuple->addIndexedItem ("n", gg_nntrk, gg_n );
		}
		else
		{
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(gg_tuple) << endmsg;
      return StatusCode::FAILURE;
		}
	}
	return StatusCode::SUCCESS;
}


StatusCode JPsi::execute()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  int runNo=eventHeader->runNumber();
  int event=eventHeader->eventNumber();
	m_time = eventHeader->time();
	if(event_proceed%1000==0)
	{
		log << MSG::DEBUG <<"run, evtnum = "
			<< runNo << " , "
			<< event <<endreq;
		std::cout << "proceed event " << event_proceed << std::endl;
		time_t t=eventHeader->time();
		cout << t << " "  << ctime(&t) << endl;
	}

	/*  Get information about reconstructed events */
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	std::vector<HepLorentzVector> p(MAX_TRACK_NUMBER);
	std::vector<Hep3Vector> p3(MAX_TRACK_NUMBER);
	InitData();

	nchtrk = evtRecEvent->totalCharged();
	nneutrk= evtRecEvent->totalNeutral();

	/************    Multihadron event and BhaBha selection ****************/
	/*  the selection is based on charged tracks */
	if(MIN_CHARGED_TRACKS<=nchtrk && nchtrk <=MAX_TRACK_NUMBER)
	{
		if(evtRecEvent->totalCharged()<3) 
			m_Signal=0;
		else 
			m_Signal=1;
		ntrk=nchtrk+nneutrk;
		double p2sum=0;
		double  Eh[2]={0, 0};
		Hep3Vector ph[2];//Two high momentum
		unsigned has_mdc_emc=0;
		/*  loop over charged track */
		for(int i = 0; i < evtRecEvent->totalCharged(); i++)
		{
			tr_idx=i+1;
			trdedx_idx=i+1;
			EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
			if(!(*itTrk)->isMdcTrackValid() || !(*itTrk)->isEmcShowerValid() ) continue;
			has_mdc_emc++;
			RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
			RecEmcShower *emcTrk = (*itTrk)->emcShower(); //Electro Magnet Calorimeer
			p[i]=HepLorentzVector(1, 1, 1, 1);
			p[i].setTheta(mdcTrk->theta());
			p[i].setPhi(mdcTrk->phi());
			p[i].setRho(mdcTrk->p());
			p[i].setE(emcTrk->energy());
			m_pt[i]=mdcTrk->p()*sin(mdcTrk->theta());
			m_E[i]=emcTrk->energy();
			m_p[i]=mdcTrk->p();
			m_px[i]=mdcTrk->px();
			m_py[i]=mdcTrk->py();
			m_pz[i]=mdcTrk->pz();
			m_theta[i]=mdcTrk->theta();
			m_phi[i]=mdcTrk->phi();
			m_M[i]=p[i].m();
			m_q[i]=mdcTrk->charge();
			m_x[i]=mdcTrk->x();
			m_y[i]=mdcTrk->y();
			m_z[i]=mdcTrk->z();
			m_ismu[i]=(*itTrk)->isMucTrackValid();
			Etotal+=m_E[i];
			/* Calculate sphericity tensor */
			for(int i=0;i<3;i++)
				for(int j=0;j<3;j++)
				{
					S[i][j]+=mdcTrk->p3()[i]*mdcTrk->p3()[j];
				}
			p2sum+=mdcTrk->p()*mdcTrk->p();
			if(fabs(m_x[i])<DELTA_X && fabs(m_y[i]) < DELTA_Y && fabs(m_z[i])< DELTA_Z) niptrk++;
			/* find two high energy track */
			if(emcTrk->energy() >= Eh[0])
			{
				Eh[0]=emcTrk->energy();
				ph[0]=mdcTrk->p3();
			}
			else
			{
				if(emcTrk->energy() >= Eh[1])
				{
					ph[1]=mdcTrk->p3();
					Eh[1]=emcTrk->energy();
				}
			}
			/* dEdx information */
			if(prop_check_dedx == 1 && (*itTrk)->isMdcDedxValid())
			{
				RecMdcDedx* dedxTrk = (*itTrk)->mdcDedx();
				m_chie[i] = dedxTrk->chiE();
				m_chimu[i] = dedxTrk->chiMu();
				m_chipi[i] = dedxTrk->chiPi();
				m_chik[i] = dedxTrk->chiK();
				m_chip[i] = dedxTrk->chiP();
				m_ghit[i] = dedxTrk->numGoodHits();
				m_thit[i] = dedxTrk->numTotalHits();
				m_probPH[i] = dedxTrk->probPH();
				m_normPH[i] = dedxTrk->normPH();
			}
		}
		if(has_mdc_emc<2) return StatusCode::SUCCESS; //at list two tracks must have drift and shower.

		//Two tracks from interaction points. The same condion for BhaBha and for multihadron
		if(USE_IPCUT && niptrk <IPTRACKS) return StatusCode::SUCCESS;

		/*  calculate angles of high energy tracks */
		double tmp = ph[0].mag()*ph[1].mag()<=0 ? -10 : (ph[0].dot(ph[1]))/(ph[0].mag()*ph[1].mag());
		m_cos_high_p=tmp;


		//normalize sphericity tensor
		for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
				S[i][j]/=p2sum;
		TMatrixDEigen Stmp(S);
		const TVectorD & eval = Stmp.GetEigenValuesRe();
		std::vector<double> v(3);
		for(int i=0;i<3;i++) v[i]=eval[i];
		std::sort(v.begin(), v.end());
		S1=v[0];
		S2=v[1];
		S3=v[2];
		m_S = 1.5*(v[0]+v[1]);
		if(!(v[0]<=v[1] && v[1]<=v[2])) //test the order of eigenvalues
		{
			cerr << "Bad sphericity" << endl;
			exit(1);
		}
		/* now fill the data */
		main_tuple->write();
		dedx_tuple->write();
		event_write++;
	}
	else 
	{
		//gamma gamma selection only two neutral tracks
		if(nneutrk==2 && nchtrk==0)
		{
			double r[2];
			for(int track = 0; track<nneutrk; track++)
			{
				gg_nntrk=track+1;
				EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + track;
				if((*itTrk)->isEmcShowerValid())
				{
					RecEmcShower *emcTrk = (*itTrk)->emcShower();
					gg_x[track] = emcTrk->x();
					gg_y[track] = emcTrk->y();
					gg_z[track] = emcTrk->z();
			 		r[track] = sqrt(gg_x[track]*gg_x[track] + gg_y[track]*gg_y[track] + gg_z[track]*gg_z[track]);
					gg_theta[track] = emcTrk->theta();
					gg_phi[track] = emcTrk->phi();
					gg_E[track]  =  emcTrk->energy();
					gg_dE[track] =  emcTrk->dE();
					gg_module[track] = emcTrk->module();
					gg_n[track] = emcTrk->numHits();
				}
			}
			//calculate colliniarity
			gg_cos = (gg_x[0]*gg_x[1] + gg_y[0]*gg_y[1] + gg_z[0]*gg_z[1])/(r[0]*r[1]);
			gg_Etotal = gg_E[0]+gg_E[1];
			gg_tuple->write();
			gg_event_writed++;
		}
	}

	event_proceed++;

	// part for ee->gg annihilation
	// big angles, two neutral track,  no charged.
	
  return StatusCode::SUCCESS;
}

StatusCode JPsi::finalize()
{
	std::cout << "Event proceed: " << event_proceed << std::endl;
	std::cout << "Event selected: " << event_write << std::endl;
	std::cout << "Selection efficiency: " << event_write/double(event_proceed) << std::endl;
  return StatusCode::SUCCESS;
}


void JPsi::InitData(void)
{
	Etotal=0;
	ntrk=0;
	nchtrk=0;
	nneutrk=0;
	niptrk=0;
	m_Signal=-1;
	m_cos_high_p=-10;
	for(unsigned i=0;i<MAX_TRACK_NUMBER;++i)
	{
		m_E[i]=-999;
		m_pt[i]=-999;
		m_p[i]=-999;
		m_px[i]=-999;
		m_py[i]=-999;
		m_pz[i]=-999;
		m_theta[i]=-999;
		m_phi[i]=-999;
		m_M[i]=-999;
		m_q[i]=-999;
		m_x[i]=-999;
		m_y[i]=-999;
		m_z[i]=-999;
		m_ismu[i]=-999;
		//dedx information
		m_chie[i] = -999;
		m_chimu[i] = -999;
		m_chipi[i] = -999;
		m_chik[i] = -999;
		m_chip[i] = -999;
		m_ghit[i] = -999;
		m_thit[i] = -999;
		m_probPH[i] = -999;
		m_normPH[i] = -999;
	}
	//sphericity initialization
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			S[i][j]=0;
	for(int i=0; i<2;i++)
	{
		gg_x [i]=-9999;
		gg_y [i]=-9999;
		gg_z [i]=-9999;
		gg_theta [i]=-9999;
		gg_phi [i]=-9999;
		gg_E [i]=-9999;
		gg_dE [i]=-9999;
		gg_module [i]=-9999;
		gg_n[i]=-9999;
	}
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
