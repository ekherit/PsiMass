#ifndef ParseFile
#define ParseFile
#include<TMath.h>
#include<iostream.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<stdio.h>
#include<fstream.h>

int GetNumRows(char *FileName,int npar)
{
    int counter=-1;
    Double_t Spool;
    ifstream test(FileName,ios::in);//|ios::floatfield);
    if(!test){
      cout<<"there is no file:"<<FileName<<endl;
      counter=0;
      test.close();
    }
    else
   {
    while(test.get()!=EOF) // find the end off file  runs.par
    {
    	counter++;             // number of lines
	for(int j=0;j<npar;j++)
	{
		test>>Spool;
	}
    }
    test.close();
   } 
    return counter;
}

void FillArrayFromFile(char* FileName,Double_t** Array,int npar,int nps)
{
     ifstream readingfile(FileName,ios::in);
     if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
     else
     {
     for(int i=0;i<nps;i++)
     {
     	Array[i]=new Double_t [npar];
     	for(int j=0;j<npar;j++)
	{
	    readingfile>>Array[i][j];
	}
     }
     cout<<endl;
     }
     readingfile.close();
}

void FillArrayFromFile(char* FileName,Double_t** Array,int npar,int nparNew,int nps)
{
    ifstream readingfile(FileName,ios::in);
    
     if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
     else
     {
     for(int i=0;i<nps;i++)
     {
     	Array[i]=new Double_t [nparNew];
     	for(int j=0;j<npar;j++)
	{
	    readingfile>>Array[i][j];
	}
     }
     cout<<endl;
     }
     readingfile.close();
}

void FillArrayMHFromFile(char* FileName,Double_t** Array,int npar,int nps)
{
     Double_t t;
     ifstream readingfile(FileName,ios::in);
     //        if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
      if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
     else
     {
  
     for(int i=0;i<nps;i++)
     {
     	Array[i]=new Double_t [npar];
     	for(int j=0;j<npar;j++)
	{
		readingfile>>Array[i][j];
	}
	t= Array[i][2];
	Array[i][2]= Array[i][3];
	Array[i][3]= Array[i][2];
     }
     cout<<endl;
     }
     readingfile.close();
}
void ReadSkeleton(char* FileName,Int_t nps,Double_t Array[][2])
{
    ifstream readingfile(FileName,ios::in);
    if(!readingfile) cout<<"there is no file:"<<FileName<<endl;
   
  
     for(Int_t i=0;i<nps;i++)
     {
     	for(Int_t j=0;j<2;j++)
	{
		readingfile>>Array[i][j];
	}
     }
     readingfile.close();
}
#endif

