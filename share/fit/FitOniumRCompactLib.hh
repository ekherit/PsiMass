//--------------------------------------------------------------------------
// File and Version Information:
// FitOniumRCompactLib.hh
//  Description:
//  Library for fitting onium resonances (short version)
// Environment:
//      Software developed for the KEDR Detector at the VEPP-4M.
//      Branch for BES analysis
// Author List:
//      Korneliy Todyshev               Originator
// Copyright Information:
//      Copyright (C) 2000               KEDR
//
//------------------------------------------------------------------------
#ifndef  FitOniumRCompactLib
#define  FitOniumRCompactLib
#include<TMath.h>
#include<iostream>
#include<stdlib.h>
#include<string>
#include<assert.h>
#include<stdio.h>
#include<fstream>
#define HALF_KNOT_COUNT 21
#define _me         0.5109982
#define _mmu      105.658367
#define _ConvConstant      389379323.0e+3 // nb*MeV^2
#define _alpha 1./137.036
#define _3part 0.33333333333
#define _4part 0.25
#define _6part 0.66666666667
#define _24part 0.0416666667
#define _MPsiPrime          3686.111
#define _MJPsi              3096.917
#define _GeeJPsi            5.55e-3
#define _GeePsiPrime        2.38e-3 //average PDG
#define _GtotJPsi           93.2e-3  //
#define _GtotPsiPrime       304e-3  //
#define _BllJPsi            0.0594
#define _BllPsiPrime        0.00743
#define _BllPsiDPrime       9.8e-6
#define _IdJPsi                0
#define _IdPsiPrime         1
#define _MethodSimple         0
#define _MethodAzimov       1
#define _MethodDIntegral    2
#define   ScaleEGr          2.0    

Double_t K_FuncRInterfPsiP(Double_t W, Double_t* parf);
Double_t K_FuncRInterfJPsi(Double_t W, Double_t* parf);
Double_t CrossSBhabhaPP(Double_t Eb,Double_t* par);
Double_t FCrSPPrimeAzimov(Double_t* Eb,Double_t* par);

Double_t CrSOniumR(Int_t Method,Int_t Id,Double_t Eb,Double_t* par);
Double_t CrSOniumRAzimov(Int_t Id,Double_t Eb,Double_t* par);
Double_t HANDLEDGAUSS(Double_t F(Double_t W,Double_t* parf),Double_t A,Double_t B,Double_t* par, Double_t eps);
Double_t HANDLEDGAUSSOLD(Double_t F(Double_t W,Double_t* parf),Double_t A,Double_t B,Double_t* par, Double_t eps);
Double_t FCrossSection(Double_t* Eb,Double_t* par);


void SeparatePointsPartNew(Int_t npini,Double_t* q,Int_t* npfin,Int_t* NpUse,Double_t* SortArray,Double_t Diff);
void SumPointsSimple(Int_t npini,Int_t npfin,Int_t* UsePoints,Double_t* v,Double_t* vn);
void SumPointsByQuant(Int_t npini,Int_t npfin,Int_t* UsePoints,Double_t* v,Double_t* q,Double_t*s,Double_t* vn,Double_t* sn,bool ZeroOpt);
int GetNumRows(char *FileName,int npar);
void FillArrayFromFile(char* FileName,Double_t** Array,int npar,int nps);
void FillArrayFromFile(char* FileName,Double_t** Array,int npar,int nparNew,int nps);
Int_t compar(const void* a,const void* b);
Int_t comparD(Int_t UpDown,Double_t a,Double_t b);
Int_t comparDRows(Int_t UpDown,Double_t* a,Double_t* b,Int_t n);
void swapD(Double_t& a,Double_t& b);
void swapI(Int_t& a,Int_t& b);
void swapDRows(Double_t* a,Double_t* b,Int_t n);

#define  idRbg       0		  
#define  idReff      1
#define  idRM        2
#define  idRSw       3
#define  idRGee      4
#define  idRTauEff   5
#define  idRFreeGee  6
#define  idRNP       7


#define  idbeta     0
#define  idM        1
#define  idW        2
#define  idSw       3
#define  idGt       4
#define  ideff      5
#define  idGee      6
#define  idDeltaE   7
#define  idDeltaEMB 8
#define  idLsqlog   9
#define  idSQB      10
#define  idsqaphi   11
#define  idXe       12
#define  idBhadr    13
#define  idNPar     14
#endif
