#ifndef MathLibraryFitTools
#define MathLibraryFitTools
#include <TMath.h>
#define sqr_2pi 0.3989422804
class TF1;

class TGraphErrors;

Double_t AGauss(Double_t *x, Double_t *par);
Double_t dEdxResolutionN(Double_t* x,Double_t* p);
Double_t dEdxResolutionL(Double_t* x,Double_t* p);
Double_t DistExp(Double_t* x,Double_t* par);
Double_t DistF(Double_t* x,Double_t* par);
Double_t DistFtest(Double_t* x,Double_t* par);
Double_t DistFN(Double_t* x,Double_t* par);
Double_t DistF9(Double_t* x,Double_t* p);
Double_t PolN(Double_t* t,Double_t* p,Int_t N);
Double_t powerN(Double_t* x, Double_t* p);
Double_t powerL(Double_t* x, Double_t* p);
Double_t powerBB(Double_t* x, Double_t* p);
Double_t powerP(Double_t* x, Double_t* p);


//Double_t Pol3(Double_t* x,Double_t* p);
void FitPoln(TF1** poln, Double_t* param, TGraphErrors* dEvsDist,Double_t* x,Double_t* y, Double_t* sigmay,Int_t npoints,Int_t* degr,Int_t NumPar);
Double_t ThetaF(Double_t* x,Double_t* par);
Double_t ThetaF9test(Double_t* x,Double_t* p);
Double_t TimeF(Double_t* x,Double_t* par);
Double_t TimeFtest(Double_t* x,Double_t* par);
Double_t ThetaF9(Double_t* x,Double_t* par);
Double_t AlphaF(Double_t* x,Double_t* par);
Double_t AlphaF3(Double_t* x,Double_t* p);
Double_t AlphaF5(Double_t* x,Double_t* p);
Double_t AlphaFtest(Double_t* x,Double_t* p);
Double_t RF6(Double_t* x,Double_t* p);
Double_t RF4(Double_t* x,Double_t* p);
Double_t RF5(Double_t* x,Double_t* p);
Double_t BetheBloch(Double_t* x,Double_t* par);
Double_t BetheBloch2(Double_t* x,Double_t* par);
Double_t BetheBloch_new(Double_t* x,Double_t* par);
Double_t BetheBloch_new1(Double_t* x,Double_t* par);
Double_t test1(Double_t x,Double_t* par);
Double_t test(Double_t x,Double_t* par);
Int_t compar(const void* a,const void* b);
Int_t comparD(Int_t UpDown,Double_t a,Double_t b);
Int_t comparDRows(Int_t UpDown,Double_t* a,Double_t* b,Int_t n);
void swapD(Double_t& a,Double_t& b);
void swapI(Int_t& a,Int_t& b);
void swapDRows(Double_t* a,Double_t* b,Int_t n);
Double_t ZF(Double_t* x, Double_t* par);
Double_t powerN(Double_t* x, Double_t* p);
Double_t powerL(Double_t* x, Double_t* p);
Double_t powerBB(Double_t* x, Double_t* p);
Double_t powerP(Double_t* x, Double_t* p);

inline Double_t sq(Double_t x){return (x*x);}
inline Double_t cub(Double_t x){return (x*x*x);}
inline Double_t sq4(Double_t x){return (x*x*x*x);}
inline Double_t sq5(Double_t x){return sq(x)*cub(x);}
inline Double_t sq6(Double_t x){return cub(x)*cub(x);}
inline Double_t sq7(Double_t x){return sq4(x)*cub(x);}
inline Double_t sq8(Double_t x){return sq4(x)*sq4(x);}
inline Double_t sq9(Double_t x){return sq4(x)*sq5(x);}
inline Double_t sq10(Double_t x){return sq5(x)*sq5(x);}
inline Double_t sq11(Double_t x){return sq5(x)*sq6(x);}

inline Double_t Min(Double_t a,Double_t b){return TMath::Min(a,b);}
inline Double_t Max(Double_t a,Double_t b){return TMath::Max(a,b);}

inline Double_t T0(){return 1;}
inline Double_t T1(Double_t x){return x;}
inline Double_t T2(Double_t x){return (2*sq(x)-1);}
inline Double_t T3(Double_t x){return (4*cub(x)-3*x);}
inline Double_t T4(Double_t x){return (8*sq4(x)-8*sq(x)+1);}
inline Double_t T5(Double_t x){return (16*sq5(x)-20*cub(x)+5*x);}
inline Double_t T6(Double_t x){return (32*sq6(x)-48*sq4(x)+18*sq(x)-1);}
inline Double_t T7(Double_t x){return (64*sq7(x)-112*sq5(x)+56*cub(x)-7*x);}
inline Double_t T8(Double_t x){return (128*sq8(x)-256*sq6(x)+160*sq4(x)-32*sq(x)+1);}
inline Double_t T9(Double_t x){return (256*sq9(x)-576*sq7(x)+432*sq5(x)-120*cub(x)+9*x);}
inline Double_t T10(Double_t x){return (512*sq10(x)-1280*sq8(x)+1120*sq6(x)-400*sq4(x)+50*sq(x)-1);}
inline Double_t T11(Double_t x){return (1024*sq11(x)-2816*sq9(x)+2816*sq7(x)-1232*sq5(x)+220*cub(x)-11*x);}
inline Double_t T12(Double_t x){return (2*x*T11(x)-T10(x));}
inline Double_t T13(Double_t x){return (2*x*T12(x)-T11(x));}
inline Double_t T14(Double_t x){return (2*x*T13(x)-T12(x));}
inline Double_t T15(Double_t x){return (2*x*T14(x)-T13(x));}
inline Double_t T16(Double_t x){return (2*x*T15(x)-T14(x));}
inline Double_t T17(Double_t x){return (2*x*T16(x)-T15(x));}
inline Double_t T18(Double_t x){return (2*x*T17(x)-T16(x));}
inline Double_t T19(Double_t x){return (2*x*T18(x)-T17(x));}
inline Double_t T20(Double_t x){return (2*x*T19(x)-T18(x));}
inline Double_t Hermite0(){return 1;}
inline Double_t Hermite1(Double_t x){return  2*x;}
inline Double_t Hermite2(Double_t x){return  4*sq(x)- 2;}
inline Double_t Hermite3(Double_t x){return  8*cub(x)-12*x;}
inline Double_t Hermite4(Double_t x){return  16*sq4(x)-48*sq(x)+12;}
inline Double_t Hermite5(Double_t x){return  32*sq5(x)-160*cub(x)+120*x;}
inline Double_t Hermite6(Double_t x){return  64*sq6(x)-480*sq4(x)+720*sq(x)-120;}
inline Double_t Hermite7(Double_t x){return  128*sq7(x)-1344*sq5(x)+3360*cub(x)-1680*x;}
inline Double_t Hermite8(Double_t x){return  256*sq8(x)-3584*sq6(x)+13440*sq4(x)-13440*sq(x)+1680;}
inline Double_t Hermite9(Double_t x){return  512*sq9(x)-9216*sq7(x)+48384*sq5(x)-80640*cub(x)+30240*x;}
inline Double_t Hermite10(Double_t x){return 1024*sq10(x)-23040*sq8(x)+161280*sq6(x)-403200*sq4(x)+302400*sq(x)-30240;}
inline Double_t Hermite11(Double_t x){return 2*(x*Hermite10(x)-10*Hermite9(x));}
inline Double_t Hermite12(Double_t x){return 2*(x*Hermite11(x)-11*Hermite10(x));}
inline Double_t Hermite13(Double_t x){return 2*(x*Hermite12(x)-12*Hermite11(x));}
inline Double_t Hermite14(Double_t x){return 2*(x*Hermite13(x)-13*Hermite12(x));}
inline Double_t Hermite15(Double_t x){return 2*(x*Hermite14(x)-14*Hermite13(x));}
inline Double_t Hermite16(Double_t x){return 2*(x*Hermite15(x)-15*Hermite14(x));}
Double_t factorial(Double_t* N);
Double_t PartPuasson(Double_t Nexpect,Double_t Nobserv);
Double_t PartMultinom(Double_t Nexpect,Double_t Nobserv);
Double_t chebysh_exp(Double_t* x,Double_t* par);
Double_t power(Double_t* x, Double_t* par);
Double_t linear(Double_t* x, Double_t* par);
Double_t linear_(Double_t x, Double_t* par);
Double_t linear1(Double_t* x, Double_t* par);
Double_t linear2(Double_t* x, Double_t* par);
Double_t linear0(Double_t* x, Double_t* par);
Double_t along(Double_t* x, Double_t* par);
Double_t alongz(Double_t* x, Double_t* par);
Double_t alongzA(Double_t* x, Double_t* par);
Double_t standard_landau(Double_t* x, Double_t* par);
Double_t standard_Landau(Double_t* x, Double_t* par);
Double_t landauFF(Double_t* x, Double_t* par);
Double_t de_loss(Double_t* x, Double_t* par);
Double_t de_loss_fit(Double_t x, Double_t* par);
Double_t ChargeR(Double_t ChargeI, Double_t* par);
Double_t dE_Real(Double_t ChargeI, Double_t* par);
Double_t Charge_R(Double_t ChargeI, Double_t* par);
Double_t signum(Double_t x);
void  ScalePoints(Double_t* a,Double_t xmin,Double_t xmax,Double_t xrmin,Double_t xrmax,Double_t xshift,Int_t nleft,Int_t nr,Int_t nright);
void devide_l(Int_t n,Double_t* a,Double_t xmin,Double_t xmax,Double_t xbet,Int_t nmin);
Double_t out_l(Int_t i,Double_t* a);
void aver_value(Int_t dim,Int_t *dnp,Double_t* v,Double_t* s,Double_t* vn,Double_t* sn);
void aver_value_scan(Int_t npe,Int_t isc,Double_t** v,Double_t** s,Double_t* vn,Double_t* sn);
void aver_value_cut(Int_t nps,Int_t ns,Double_t* v,Double_t* s,Double_t* vn,Double_t* sn);
Double_t dist_exp(Double_t* x, Double_t* par);
Float_t define_bins(Int_t nslices,const Int_t* scaleweight,
	Float_t xmin,Float_t xmax,Float_t* bincenter,Int_t* defplace);
Float_t define_bins_(Int_t nslices,Float_t xmin,Float_t xmax,Float_t* bincenter);
Double_t BB(Double_t* x,Double_t* par);
Double_t BBfix(Double_t x);
Double_t BB3(Double_t* x,Double_t* par);
// BaBar functions
Double_t ChargeCorrCoupal(Double_t chargein,Double_t sindip, Double_t* CoupalPar);
Double_t ChargeCorrCheby(Double_t chargein,Double_t sindip, Double_t* ChebyPar);
Double_t ChebyP4(Double_t* sindip, Double_t* par);
Double_t ChebyP3(Double_t* z, Double_t* par);
Float_t ChebyP3Fix(Float_t z, Float_t* par);
Float_t ChebyP4Fix(Float_t sindip, Double_t* par);
Double_t ChargeAfterSat(Double_t sindip,Double_t rhoraw, Double_t* par);
Double_t ChargeAfterSat_(Double_t* x, Double_t* par);
Double_t ChargeAfterSatFix(Double_t* x, Double_t* par);
Double_t CalibrationCharge1(Double_t* x, Double_t* p);
Double_t CalibrationCharge2(Double_t* x, Double_t* p);
Double_t CalibrationCharge3(Double_t* x, Double_t* p);
Double_t CalibrationChargeGlobal(Double_t* x, Double_t* p);
Double_t squareF(Double_t* x, Double_t* p);
Double_t XtC(Double_t* x, Double_t* par);
Double_t dEdxResolutionM(Double_t* x,Double_t* par);
Double_t TMFCut(Double_t* x, Double_t* par);


#endif
