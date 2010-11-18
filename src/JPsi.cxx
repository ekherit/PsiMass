/*
 * =====================================================================================
 *
 *       Filename:  TauEMU.cxx
 *
 *    Description:  My bha bha implementation algorithm
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
  declareProperty("Delta_x", prop_delta_x = 1.0); //cm?
  declareProperty("Delta_y", prop_delta_y = 1.0); //cm?
  declareProperty("Delta_z", prop_delta_z = 10.0); //cm?
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
  
  StatusCode status;
	NTuplePtr my_nt(ntupleSvc(), "FILE1/chtr");
	if(my_nt) chtr_tuple=my_nt;
	else
	{
		chtr_tuple = ntupleSvc()->book("FILE1/chtr", CLID_ColumnWiseTuple, "Charged rack");
		if(chtr_tuple)
		{
      status = chtr_tuple->addItem ("nchtr", tr_idx, 0, MAX_TRACK_NUMBER);
      status = chtr_tuple->addIndexedItem ("E", tr_idx, m_E );
      status = chtr_tuple->addIndexedItem ("pt", tr_idx, m_pt );
      status = chtr_tuple->addIndexedItem ("M", tr_idx, m_M );
      status = chtr_tuple->addIndexedItem ("q", tr_idx, m_q );
      status = chtr_tuple->addIndexedItem ("x", tr_idx, m_x );
      status = chtr_tuple->addIndexedItem ("y", tr_idx, m_y );
      status = chtr_tuple->addIndexedItem ("z", tr_idx, m_z );
      status = chtr_tuple->addIndexedItem ("ismu", tr_idx, m_ismu );
		}
		else
		{
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(chtr_tuple) << endmsg;
      return StatusCode::FAILURE;
		}
	}
	NTuplePtr mainnt(ntupleSvc(), "FILE1/main");
	if(mainnt) main_tuple=mainnt;
	else
	{
		main_tuple = ntupleSvc()->book("FILE1/main", CLID_ColumnWiseTuple, "Event information");
		if(main_tuple)
		{
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
		}
		else
		{
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(main_tuple) << endmsg;
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
  log << MSG::DEBUG <<"run, evtnum = "
      << runNo << " , "
      << event <<endreq;
	if(event_proceed%1000==0)
		std::cout << "proceed event " << event_proceed << std::endl;
	event_proceed++;

	/*  Get information about reconstructed events */
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);
	std::vector<HepLorentzVector> p(MAX_TRACK_NUMBER);
	std::vector<Hep3Vector> p3(MAX_TRACK_NUMBER);
	InitData();
	/* Multihadron selecion and Bhabha*/
	if(evtRecEvent->totalCharged()<MIN_CHARGED_TRACKS || evtRecEvent->totalCharged()>=MAX_TRACK_NUMBER) return StatusCode::SUCCESS;
	if(evtRecEvent->totalCharged()<3) 
		m_Signal=0;
	else 
		m_Signal=1;
	nchtrk = evtRecEvent->totalCharged();
	nneutrk= evtRecEvent->totalCharged();
	ntrk=nchtrk+nneutrk;
	double p2sum=0;
	/*  loop over charged track */
  for(int i = 0; i < evtRecEvent->totalCharged(); i++)
	{
		tr_idx=i;
    EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
    if(!(*itTrk)->isMdcTrackValid()) continue;
    if(!(*itTrk)->isEmcShowerValid()) continue;
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
    RecEmcShower *emcTrk = (*itTrk)->emcShower(); //Electro Magnet Calorimeer
		p[i]=HepLorentzVector(1, 1, 1, 1);
		p[i].setTheta(mdcTrk->theta());
		p[i].setPhi(mdcTrk->phi());
		p[i].setRho(mdcTrk->p());
		p[i].setE(emcTrk->energy());
		m_pt[i]=mdcTrk->p()*sin(mdcTrk->theta());
		m_E[i]=emcTrk->energy();
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
		if(fabs(m_x[i])<prop_delta_x && fabs(m_y[i]) < prop_delta_y && fabs(m_z[i])< prop_delta_z) niptrk++;
	}

	//Two tracks from interaction points.
	//The same condion for BhaBha and for multihadron
	if(USE_IPCUT && niptrk <IPTRACKS) return StatusCode::SUCCESS;

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
	chtr_tuple->write();
	main_tuple->write();
	event_write++;
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
	for(unsigned i=0;i<MAX_TRACK_NUMBER;++i)
	{
		m_E[i]=-999;
		m_pt[i]=-999;
		m_M[i]=-999;
		m_q[i]=-999;
		m_x[i]=-999;
		m_y[i]=-999;
		m_z[i]=-999;
		m_ismu[i]=0;
	}
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			S[i][j]=0;
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
