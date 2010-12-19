#ifndef JpsiLibrary
#define JpsiLibrary
#include <TMath.h>
#define nsequent 300
#define me  0.510999
#define alfa 1./137.07
#define MUpsilon1S         9460.30
#define GeeUpsilon1S       1.314e-3
#define GtotUpsilon1S      0.053
#define MPsiPrime          3686.111
#define MJPsi              3096.917
#define MPsiDoublePrime    3770
#define GeeJPsi            5.3963e-3
#define GeePsiPrime        2.12155e-3
#define GeePsiDoublePrime  0.26e-3
#define BllJPsi            0.0593
#define BllPsiPrime        0.0090
#define MDc      1869.3
#define MD0      1864.5
Double_t FuncR(Double_t W, Double_t* parf);
Double_t FuncRDoublePrime(Double_t W, Double_t* parf);
Double_t FuncRUpsilon(Double_t W, Double_t* parf);
Double_t HANDLE_DGAUSS(Double_t F(Double_t W,Double_t* parf),Double_t A,Double_t B,Double_t* par, Double_t eps);
long double HANDLE_DGAUSS2(long double F(long double W,long double* parf),long double A,long double B,long double* par, long double eps);

Double_t GrIntPsiPrime(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par);

Double_t M_h_f_(Double_t* a,Double_t *b,Double_t *x);
Double_t M_h_f_d(Double_t* a,Double_t *b,Double_t *x);
Double_t fsn_(Double_t* t,Double_t* z,Int_t n,Int_t s);
Double_t fsn_d(Double_t* t,Double_t* z,Int_t n,Int_t s);
Double_t Fzt_Aprox(Double_t* z,Double_t *t);
Double_t Fzt_(Double_t* z,Double_t *t);
Double_t Fzt_1(Double_t z,Double_t *t);
Double_t Fzt_d(Double_t* z,Double_t *t);
Double_t xsecnbsim(Double_t* Eb,Double_t* par);
Double_t xs_bhabha(Double_t Eb,Double_t iscan,Double_t* par);
Double_t xsbhabha_(Double_t Eb,Double_t* par);
Double_t xsbhabha(Double_t Eb,Double_t iscan,Double_t* par);
Double_t xsecnbScansSeparate(Double_t Eb,Double_t* par);
Double_t xsecnbYpsilon(Double_t Eb,Double_t* par);
Double_t GrIntYpsilon(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par);
Double_t GrInt(Double_t Eb,Double_t* par);
Double_t GrIntNew(Double_t Eb,Double_t epsilon_, Double_t RangeE,Double_t* par);
Double_t RPsiDoublePrime(Double_t Eb,Double_t epsilon_,Double_t RangeE,Double_t* par);
Double_t RPsi(Double_t* Eb,Double_t* par);
//Double_t HANDLE_DGAUSS(Double_t F(Double_t W,Double_t* parf,Double_t A,Double_t B,Double_t* par, Double_t eps);
Double_t xsecnbScansSeparatePsiPrime(Double_t Eb,Double_t* par);
Double_t xsecnbScansSeparateBeta(Double_t Eb,Double_t* par);
Double_t xsecnbScansSeparateErf(Double_t Eb,Double_t* par);
Double_t xsecnbFit(Double_t* Eb,Double_t* par);
Double_t xsecnbFitPsiPrime(Double_t* Eb,Double_t* par);
Double_t xsecnbFitNew(Double_t* Eb,Double_t* par);
Double_t xsecnbFitNew_(Double_t* Eb,Double_t* par);
Double_t xsecnb(Double_t* Eb,Double_t* par);
Double_t xsecnbEI(Double_t Eb,Double_t I,Double_t* par);
Double_t xsecnbEIFix(Double_t* x,Double_t* par);
Double_t xsecnb_gaus(Double_t* Eb,Double_t* par);
Double_t xsecnb_d_p3(Double_t* Eb,Double_t* par);
Double_t xsecnb_0(Double_t* Eb,Double_t* par);
Double_t xsecnb_d(Double_t* Eb,Double_t* par);
Double_t SEvsI(Double_t icur,Double_t *par);
#endif
