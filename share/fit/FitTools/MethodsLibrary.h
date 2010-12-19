#ifndef MethodsLibrary
#define MethodsLibrary
#include "FitTools/MathLibrary.h"

bool ApproxRoot(Double_t* r,Double_t F(Double_t t,Double_t* p),Double_t t0,Double_t t1,Double_t* par,Int_t it,Double_t eps)
{
    bool root=false;
    Int_t iq=0;
    Int_t ic=0;
    Double_t tmin,tmax,dt;
    Double_t tf,tin=t1,tin0=t0,tcmax,tcmin;
    Double_t fin,fin0,ff;
   // tmin=Min(tin,tin0);
  //  tmax=Max(tin,tin0);
    tmin=(tin < tin0) ? tin : tin0;
    tmax=(tin > tin0) ? tin : tin0;
    tcmax=tmax;
    tcmin=tmin;
  
    dt=fabs(tmax-tmin);
    fin=F(tin,par);
    fin0=F(tin0,par);
    while(iq==0)
    {
        tf=tin-(tin-tin0)/(fin-fin0)*fin;
        if(fabs(tf-tin)<eps||ic>it
      //     ||(Max(tin,tin0)>tmax&&fabs(tin-tin0)>dt)
       //    ||(Min(tin,tin0)<tmin&&fabs(tin-tin0)>dt)) iq=1;
         ||(tcmax>tmax&&fabs(tin-tin0)>dt)
           ||(tcmin<tmin&&fabs(tin-tin0)>dt)) iq=1;
        
        else
        {
          tin0=tin;
          tin=tf;
          fin0=fin;
          fin=F(tf,par);
          ic++;
          tcmin=(tin < tin0) ? tin : tin0;
          tcmax=(tin > tin0) ? tin : tin0;
   

        }
    }
    if(fabs(tf-tin)<eps)
    {
        root=true;
        r[0]=tf;
    }
    return root;
}

void MultiRoot(Double_t* roots,Int_t nr,Double_t F(Double_t t,Double_t* p),Double_t t0,Double_t t1,Double_t* par,Int_t it,Double_t eps)
{
    Double_t r=0;
    bool     root;
    Double_t step;
    Int_t iq=0;
    Int_t cr=0;
    step=(t1-t0)/nr;
    for(Int_t j=0;j<nr;j++)
    {
        root=ApproxRoot(&r,F,t0+step*j,t1+step*(j+1),par,it,eps);
        roots[cr]=0;
        if(root==true)
        {
            iq=0;
           // cout<<"root:"<<r<<"x:"<<t0+step*j<<endl;
            if(cr>=1)
            {
                for(Int_t i=0;i<cr;i++)
                {
                    if(roots[i]!=0&&fabs(roots[i]-r)<eps) iq=1;
                }
            }
            if(iq==0)
            {
                roots[cr]=r;
                cr++;
            }
        }

    }
}
/*
bool ApproxRootLineCircle(Double_t* r,Double_t F(Double_t t,Double_t* p),
                          Double_t t0,Double_t t1,Double_t* par,
                          Int_t it,Double_t eps)
{
    bool root=false;
    Int_t iq=0;
    Int_t ic=0;
    Double_t tmin,tmax,dt;
    Double_t tf,tin=t1,tin0=t0;
    Double_t fin,fin0,ff;
    tmin=Min(tin,tin0);
    tmax=Max(tin,tin0);
    dt=fabs(tmax-tmin);
    fin=F(tin,par);
    fin0=F(tin0,par);
    while(iq==0)
    {
        tf=tin-(tin-tin0)/(fin-fin0)*fin;
        if(fabs(tf-tin)<eps||ic>it
           ||(Max(tin,tin0)>tmax&&fabs(tin-tin0)>dt)
           ||(Min(tin,tin0)<tmin&&fabs(tin-tin0)>dt)) iq=1;
        else
        {
          tin0=tin;
          tin=tf;
          fin0=fin;
          fin=F(tf,par);
          ic++;
        }
    }
    if(fabs(tf-tin)<eps)
    {
        root=true;
        r[0]=tf;
    }
    return root;
}

void MultiRoot(Double_t* roots,Int_t nr,Double_t F(Double_t t,Double_t* p),Double_t t0,Double_t t1,Double_t* par,Int_t it,Double_t eps)
{
    Double_t r=0;
    bool     root;
    Double_t step;
    Int_t iq=0;
    Int_t cr=0;
    step=(t1-t0)/nr;
    for(Int_t j=0;j<nr;j++)
    {
        root=ApproxRoot(&r,F,t0+step*j,t1+step*(j+1),par,it,eps);
        roots[cr]=0;
        if(root==true)
        {
            iq=0;
            cout<<"root:"<<r<<"x:"<<t0+step*j<<endl;
            if(cr>=1)
            {
                for(Int_t i=0;i<cr;i++)
                {
                    if(roots[i]!=0&&fabs(roots[i]-r)<eps) iq=1;
                }
            }
            if(iq==0)
            {
                roots[cr]=r;
                cr++;
            }
        }

    }
    }
*/
#endif



