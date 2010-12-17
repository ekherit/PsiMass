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


const double PI_MESON_MASS=0.13957018; //GeV

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
    declareProperty("MAX_TRACK_NUMBER", MAX_TRACK_NUMBER=30); //maximum number of charged tracks
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
            status=main_tuple->addItem("ntrack", m_ntrack);
            status=main_tuple->addItem("nchtrk", m_nchtr);
            status=main_tuple->addItem("nneutrk", m_nneutr);
            status=main_tuple->addItem("Etotal", m_Etotal);
            status=main_tuple->addItem("Eemc", m_Eemc);
        }
        else
        {
            log << MSG::ERROR << "    Cannot book N-tuple:" << long(main_tuple) << endmsg;
            return StatusCode::FAILURE;
        }
    }

    /*  MDC ntuple intialization */
    NTuplePtr nt_mdc(ntupleSvc(), "FILE1/mdc");
    if(nt_mdc) mdc_tuple=nt_mdc;
    else
    {
        mdc_tuple = ntupleSvc()->book("FILE1/mdc", CLID_ColumnWiseTuple, "Charged  tracks");
        if(mdc_tuple)
        {
            status = mdc_tuple->addItem ("ntrack", mdc.ntrack, 0, MAX_TRACK_NUMBER);
            status = mdc_tuple->addItem ("Emdc", mdc.Emdc);
            status = mdc_tuple->addItem ("Eemc", mdc.Eemc);
            status = mdc_tuple->addItem ("nemc", mdc.nemc);
            status = mdc_tuple->addItem ("nip", mdc.nip);
            status = mdc_tuple->addItem ("idx1", mdc.idx1);
            status = mdc_tuple->addItem ("idx2", mdc.idx2);
            status = mdc_tuple->addItem ("hp_cos", mdc.hp_cos);
            status = mdc_tuple->addItem ("pt50", mdc.pt50);
            status = mdc_tuple->addItem ("pt100", mdc.pt100);
            status = mdc_tuple->addIndexedItem ("p", mdc.ntrack, mdc.p);
            status = mdc_tuple->addIndexedItem ("pt", mdc.ntrack, mdc.pt);
            status = mdc_tuple->addIndexedItem ("px", mdc.ntrack, mdc.px);
            status = mdc_tuple->addIndexedItem ("py", mdc.ntrack, mdc.py);
            status = mdc_tuple->addIndexedItem ("pz", mdc.ntrack, mdc.pz);
            status = mdc_tuple->addIndexedItem ("theta", mdc.ntrack, mdc.theta);
            status = mdc_tuple->addIndexedItem ("phi", mdc.ntrack, mdc.phi);
            status = mdc_tuple->addIndexedItem ("x", mdc.ntrack, mdc.x);
            status = mdc_tuple->addIndexedItem ("y", mdc.ntrack, mdc.y);
            status = mdc_tuple->addIndexedItem ("z", mdc.ntrack, mdc.z);
            status = mdc_tuple->addIndexedItem ("q", mdc.ntrack, mdc.q);
            status = mdc_tuple->addIndexedItem ("emc", mdc.ntrack, mdc.isemc);
            status = mdc_tuple->addIndexedItem ("E", mdc.ntrack, mdc.E);
            status = mdc_tuple->addIndexedItem ("dE", mdc.ntrack, mdc.dE);
            status = mdc_tuple->addIndexedItem ("ncrstl", mdc.ntrack, mdc.ncrstl);
            status = mdc_tuple->addIndexedItem ("cellId", mdc.ntrack, mdc.cellId);
            status = mdc_tuple->addIndexedItem ("status", mdc.ntrack, mdc.status);
            status = mdc_tuple->addIndexedItem ("module", mdc.ntrack, mdc.module);
            status = mdc_tuple->addIndexedItem ("M", mdc.ntrack, mdc.M);
            status = mdc_tuple->addIndexedItem ("mu", mdc.ntrack, mdc.ismu);
            status = mdc_tuple->addIndexedItem ("X", mdc.ntrack, mdc.X);
            status = mdc_tuple->addIndexedItem ("Y", mdc.ntrack, mdc.Y);
            status = mdc_tuple->addIndexedItem ("Z", mdc.ntrack, mdc.Z);

            /*  sphericity part */
            status = mdc_tuple->addItem("S", mdc.S);
        }
        else
        {
            log << MSG::ERROR << "    Cannot book N-tuple:" << long(mdc_tuple) << endmsg;
            return StatusCode::FAILURE;
        }
    }

    NTuplePtr nt_emc(ntupleSvc(), "FILE1/emc");
    if(nt_emc) emc_tuple=nt_emc;
    else
    {
        emc_tuple = ntupleSvc()->book("FILE1/emc", CLID_ColumnWiseTuple, "Netutral track");
        if(emc_tuple)
        {
            status = emc_tuple->addItem ("ntrack", emc.ntrack, 0, MAX_TRACK_NUMBER);
            status = emc_tuple->addItem ("Etotal", emc.Etotal);
            status = emc_tuple->addIndexedItem ("model", emc.ntrack, emc.module );
            status = emc_tuple->addIndexedItem ("status", emc.ntrack, emc.status );
            status = emc_tuple->addIndexedItem ("ncrstl", emc.ntrack, emc.ncrstl );
            status = emc_tuple->addIndexedItem ("cellId", emc.ntrack, emc.cellId );
            status = emc_tuple->addIndexedItem ("x", emc.ntrack, emc.x );
            status = emc_tuple->addIndexedItem ("y", emc.ntrack, emc.y );
            status = emc_tuple->addIndexedItem ("z", emc.ntrack, emc.z );
            status = emc_tuple->addIndexedItem ("E", emc.ntrack, emc.E );
            status = emc_tuple->addIndexedItem ("dE", emc.ntrack, emc.dE );
            status = emc_tuple->addIndexedItem ("theta", emc.ntrack, emc.theta );
            status = emc_tuple->addIndexedItem ("phi", emc.ntrack, emc.phi );
        }
        else
        {
            log << MSG::ERROR << "    Cannot book N-tuple:" << long(emc_tuple) << endmsg;
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
            status = dedx_tuple->addIndexedItem ("pid", trdedx_idx, m_pid );
            status = dedx_tuple->addIndexedItem ("chie", trdedx_idx, m_chie );
            status = dedx_tuple->addIndexedItem ("chimu",trdedx_idx, m_chimu );
            status = dedx_tuple->addIndexedItem ("chipi", trdedx_idx, m_chipi );
            status = dedx_tuple->addIndexedItem ("chik", trdedx_idx, m_chik );
            status = dedx_tuple->addIndexedItem ("chip", trdedx_idx, m_chip );
            status = dedx_tuple->addIndexedItem ("ghit", trdedx_idx, m_ghit );
            status = dedx_tuple->addIndexedItem ("thit", trdedx_idx, m_thit );
            status = dedx_tuple->addIndexedItem ("probPH", trdedx_idx, m_probPH );
            status = dedx_tuple->addIndexedItem ("normPH", trdedx_idx, m_normPH );
            status = dedx_tuple->addIndexedItem ("dedx_e", trdedx_idx, m_dedx_e );
            status = dedx_tuple->addIndexedItem ("dedx_mu", trdedx_idx, m_dedx_mu );
            status = dedx_tuple->addIndexedItem ("dedx_pi", trdedx_idx, m_dedx_pi );
            status = dedx_tuple->addIndexedItem ("dedx_K", trdedx_idx, m_dedx_K );
            status = dedx_tuple->addIndexedItem ("dedx_p", trdedx_idx, m_dedx_p );
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
    time_t t=eventHeader->time();
    if(event_proceed%1000==0)
    {
        std::cout << "proceed event " << event_proceed << std::endl;
    }
    event_proceed++;
    if(event_proceed<45000) return StatusCode::SUCCESS;

    /*  Get information about reconstructed events */
    SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(), EventModel::EvtRec::EvtRecEvent);
    SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),  EventModel::EvtRec::EvtRecTrackCol);

    InitData(evtRecEvent->totalCharged(), evtRecEvent->totalNeutral());

    /************    Multihadron event and BhaBha selection ****************/
    /*  the selection is based on charged tracks */
    if(MIN_CHARGED_TRACKS<=evtRecEvent->totalCharged() && evtRecEvent->totalCharged() <=MAX_TRACK_NUMBER)
    {
        double p2sum=0;
        double  Eh[2]={0, 0};
        Hep3Vector ph[2];//Two high momentum
        /*  loop over charged track */
        //mdc.ntrack=evtRecEvent->totalCharged();
        mdc.ntrack=0;
        bool ispt50=true;
        bool ispt100=true;
        for(int i = 0; i < evtRecEvent->totalCharged(); i++)
        {
            EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + i;
            mdc.ntrack=i+1;
            if(!(*itTrk)->isMdcTrackValid()) continue; 
            RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();  //main drift chambe
            mdc.p[i]=mdcTrk->p();
            mdc.pt[i]=mdcTrk->p()*sin(mdcTrk->theta());
            ispt50 = ispt50 && mdc.pt[i]>0.05;
            ispt100 = ispt100 && mdc.pt[i]>0.1;
            mdc.px[i]=mdcTrk->px();
            mdc.py[i]=mdcTrk->py();
            mdc.pz[i]=mdcTrk->pz();
            mdc.theta[i]=mdcTrk->theta();
            mdc.phi[i]=mdcTrk->phi();
            mdc.q[i]=mdcTrk->charge();
            mdc.x[i]=mdcTrk->x();
            mdc.y[i]=mdcTrk->y();
            mdc.z[i]=mdcTrk->z();
            mdc.X[i]=mdcTrk->getVX0();
            mdc.Y[i]=mdcTrk->getVY0();
            mdc.Z[i]=mdcTrk->getVZ0();
            mdc.Emdc+=sqrt(mdc.p[i]*mdc.p[i]+PI_MESON_MASS*PI_MESON_MASS);
            /* Calculate sphericity tensor */
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                {
                    S[i][j]+=mdcTrk->p3()[i]*mdcTrk->p3()[j];
                }
            p2sum+=mdcTrk->p()*mdcTrk->p();

            /*  If we have EMC information for this charged track */
            mdc.isemc[i]=(*itTrk)->isEmcShowerValid();
            if(mdc.isemc[i]) 
            {
                mdc.nemc++;//increase number of emc clasters
                RecEmcShower *emcTrk = (*itTrk)->emcShower(); //Electro Magnet Calorimeer

                mdc.E[i]=emcTrk->energy();
                mdc.dE[i]=emcTrk->dE();
                mdc.ncrstl[i]=emcTrk->numHits();
                mdc.status[i]=emcTrk->status();
                mdc.cellId[i]=emcTrk->cellId();
                mdc.module[i]=emcTrk->module();
                mdc.Eemc+=mdc.E[i];

                HepLorentzVector P(mdc.px[i], mdc.py[i], mdc.pz[i], mdc.E[i]);
                mdc.M[i]=P.m();


                /* find two high energy track */
                if(emcTrk->energy() >= Eh[0])
                {
                    //if we'v found energy more then highest finded energy
                    //then we should save old value to lower one.
                    long tmp=mdc.idx1;
                    mdc.idx2=tmp;
                    Eh[1]=Eh[0];
                    ph[1]=ph[0];
                    mdc.idx1=i;
                    Eh[0]=emcTrk->energy();
                    ph[0]=mdcTrk->p3();
                }
                else
                {
                    if(emcTrk->energy() >= Eh[1])
                    {
                        mdc.idx2=i;
                        ph[1]=mdcTrk->p3();
                        Eh[1]=emcTrk->energy();
                    }
                }
            }


            /* Check muon system information for this track */
            mdc.ismu[i]=(*itTrk)->isMucTrackValid();

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
                m_dedx_e[i] = dedxTrk->getDedxExpect(0);
                m_dedx_mu[i] = dedxTrk->getDedxExpect(1);
                m_dedx_pi[i] = dedxTrk->getDedxExpect(2);
                m_dedx_K[i] = dedxTrk->getDedxExpect(3);
                m_dedx_p[i] = dedxTrk->getDedxExpect(4);
                m_pid[i]=dedxTrk->particleId();
            }
        }

        mdc.pt50 = ispt50;
        mdc.pt100 = ispt100;

        /* Use data at least two charged track with signal in EMC */
        if(mdc.nemc<2) return StatusCode::SUCCESS;

        //Two tracks from interaction points. The same condion for BhaBha and for multihadron
        if( USE_IPCUT 
                && fabs(mdc.x[mdc.idx1]) > DELTA_X && fabs(mdc.y[mdc.idx1]) > DELTA_Y && fabs(mdc.z[mdc.idx1]) > DELTA_Z
                && fabs(mdc.x[mdc.idx2]) > DELTA_X && fabs(mdc.y[mdc.idx2]) > DELTA_Y && fabs(mdc.z[mdc.idx2]) > DELTA_Z
          ) return StatusCode::SUCCESS;
        /*  calculate angles of high energy tracks */
        double tmp = ph[0].mag()*ph[1].mag()<=0 ? -10 : (ph[0].dot(ph[1]))/(ph[0].mag()*ph[1].mag());
        mdc.hp_cos=tmp;


        //normalize sphericity tensor
        for(int i=0;i<3;i++)
            for(int j=0;j<3;j++)
                S[i][j]/=p2sum;

        TMatrixDEigen Stmp(S);
        const TVectorD & eval = Stmp.GetEigenValuesRe();
        std::vector<double> v(3);
        for(int i=0;i<3;i++) v[i]=eval[i];
        std::sort(v.begin(), v.end());
        mdc.S = 1.5*(v[0]+v[1]);
        if(!(v[0]<=v[1] && v[1]<=v[2])) //test the order of eigenvalues
        {
            cerr << "Bad sphericity" << endl;
            exit(1);
        }

        /*  fill data for neutral tracks */
        emc.ntrack=evtRecEvent->totalNeutral();
        int track=0; //index for neutral tracks
        emc.Etotal=0;
        cout << event_proceed << " ncharged=" << evtRecEvent->totalCharged() << " nneutral=" << evtRecEvent->totalNeutral() << endl;
        for(int idx = evtRecEvent->totalCharged(); idx<evtRecEvent->totalTracks(); idx++, track++)
        {
            cout << "idx = " << idx << " track=" << track << endl;
            EvtRecTrackIterator itTrk=evtRecTrkCol->begin() + track;
            //emc.ntrack=track+1;
            if(!(*itTrk)->isEmcShowerValid()) continue;
            RecEmcShower *emcTrk = (*itTrk)->emcShower();
            emc.status[track] = emcTrk->status();
            emc.ncrstl[track] = emcTrk->numHits();
            emc.cellId[track] = emcTrk->numHits();
            emc.module[track] = emcTrk->module();
            emc.x[track] = emcTrk->x();
            emc.y[track] = emcTrk->y();
            emc.z[track] = emcTrk->z();
            emc.theta[track] = emcTrk->theta();
            emc.phi[track] = emcTrk->phi();
            emc.E[track]  =  emcTrk->energy();
            emc.dE[track] =  emcTrk->dE();
            emc.Etotal+=emcTrk->energy();
        }

        cout << "After ECM hit" << endl;

        m_nchtr=evtRecEvent->totalCharged();
        m_nneutr=evtRecEvent->totalNeutral();
        m_ntrack=evtRecEvent->totalCharged()+evtRecEvent->totalNeutral();
        m_Etotal = emc.Etotal+mdc.Emdc;
        m_Eemc = emc.Etotal+mdc.Eemc;
        m_time = eventHeader->time();
        cout << "Before tuple write" << endl;

        /* now fill the data */
        main_tuple->write();
        cout << "Before dedx tuple write" << endl;
        dedx_tuple->write();
        cout << "Before mdc tuple write" << endl;
        mdc_tuple->write();
        cout << "Before emc tuple write" << endl;
        emc_tuple->write();
        event_write++;
    }
    else 
    {
        //gamma gamma selection only two neutral tracks
        // part for ee->gg annihilation
        // big angles, two neutral track,  no charged.
        if(evtRecEvent->totalNeutral()==2 && evtRecEvent->totalCharged()==0)
        {
            double r[2];
            for(int track = 0; track<evtRecEvent->totalNeutral(); track++)
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
            if(gg_Etotal>1) 
            {
                gg_tuple->write();
                gg_event_writed++;
            }
        }
    }
    return StatusCode::SUCCESS;
}

StatusCode JPsi::finalize()
{
    std::cout << "Event proceed: " << event_proceed << std::endl;
    std::cout << "Event selected: " << event_write << std::endl;
    std::cout << "Selection efficiency: " << event_write/double(event_proceed) << std::endl;
    return StatusCode::SUCCESS;
}


void JPsi::InitData(long nchtrack, long nneutrack)
{
    m_ntrack=nchtrack+nneutrack;
    m_nchtr=nchtrack;
    m_nneutr=nneutrack;
    m_Etotal=0;
    m_Eemc=0;
    //mdc track informaion init
    mdc.nemc=0;
    mdc.nip=0;
    mdc.Eemc=0;
    mdc.Emdc=0;
    mdc.S=0;
    mdc.idx1=-1000;
    mdc.idx2=-1000;
    mdc.hp_cos=-1000;
    mdc.pt50=-1000;
    mdc.pt100=-1000;
    mdc.ntrack=0;
    for(int i=0;i<MAX_TRACK_NUMBER; i++)
    {
        mdc.p[i]=-1000;
        mdc.px[i]=-1000;
        mdc.py[i]=-1000;
        mdc.pz[i]=-1000;
        mdc.pt[i]=-1000;
        mdc.x[i]=-1000;
        mdc.y[i]=-1000;
        mdc.z[i]=-1000;
        mdc.theta[i]=-1000;
        mdc.phi[i]=-1000;
        mdc.q[i]=-1000;
        mdc.isemc[i]=-1000;
        mdc.ncrstl[i]=-1000;
        mdc.cellId[i]=-1000;
        mdc.status[i]=-1000;
        mdc.module[i]=-1000;
        mdc.E[i]=-1000;
        mdc.dE[i]=-1000;
        mdc.M[i]=-1000;
        mdc.ismu[i]=-1000;
        mdc.X[i]=-1000;
        mdc.Y[i]=-1000;
        mdc.Z[i]=-1000;

        //dedx information
        m_pid[i]=-1000;
        m_chie[i] = -1000;
        m_chimu[i] = -1000;
        m_chipi[i] = -1000;
        m_chik[i] = -1000;
        m_chip[i] = -1000;
        m_ghit[i] = -1000;
        m_thit[i] = -1000;
        m_probPH[i] = -1000;
        m_normPH[i] = -1000;
        m_dedx_e[i]=-1000;
        m_dedx_mu[i]=-1000;
        m_dedx_pi[i]=-1000;
        m_dedx_K[i]=-1000;
        m_dedx_p[i]=-1000;
    }
    //sphericity initialization
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            S[i][j]=0;

    // emc information init.
    emc.ntrack=0;
    emc.Etotal=0;
    emc.ntrack=nneutrack;
    for(int i=0;i<nneutrack;i++)
    {
        emc.status[i]=-1000;
        emc.ncrstl[i]=-1000;
        emc.cellId[i]=-1000;
        emc.module[i]=-1000;
        emc.E[i]=-1000;
        emc.dE[i]=-1000;
        emc.x[i]=-1000;
        emc.y[i]=-1000;
        emc.z[i]=-1000;
        emc.theta[i]=-1000;
        emc.phi[i]=-1000;
    }
    for(int i=0; i<2;i++)
    {
        gg_x [i]=-1000;
        gg_y [i]=-1000;
        gg_z [i]=-1000;
        gg_theta [i]=-1000;
        gg_phi [i]=-1000;
        gg_E [i]=-1000;
        gg_dE [i]=-1000;
        gg_module [i]=-1000;
        gg_n[i]=-1000;
    }
}
// for particle id look /ihepbatch/bes/alex/workarea/Analysis/Physics/PsiPrime/G2MuMuAlg-00-00-01/PipiJpsiAlg/src
