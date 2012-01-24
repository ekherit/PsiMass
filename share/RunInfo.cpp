/*
 * =====================================================================================
 *
 *       Filename:  RunInfo.cpp
 *
 *    Description:  Utils to read runinfo files.
 *
 *        Version:  1.0
 *        Created:  20.01.2012 11:05:08
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */

#include "RunInfo.h"
#include "energy.h"
#include <vector>
#include <fstream>
#include <iostream>



class time_string
{
  tm t;
  public:
  char str[20];
  time_t unixtime;
  time_string(void)
  {
    str[0]=0;
    unixtime=0;
  }
  time_string(const char * s)
  {
    set_unixtime(s);
  }

  void set_unixtime(const char *s)
  {
    strptime(s,"%Y-%m-%d %T",&t);
    unixtime=mktime(&t);
  }
  void set_unixtime(void)
  {
    set_unixtime(str);
  }

};


inline std::istream & operator>>(std::istream & is, time_string & t)
{
  do
  {
    t.str[0]=is.get();
  } while(!isdigit(t.str[0]));
  is.get(&t.str[1],19);
  t.str[19]=0;
  t.set_unixtime();
  return is;
}



//void read_run_info(const char * filename, std::vector<RunInfo_t> & RI)
void read_run_info(const char * filename, RunInfo_t * RI, unsigned & RIsize)
{
  std::ifstream file(filename);
  if(!file)
  {
    std::cerr << "Unable to open file " << filename << endl;
    return;
  }
  RunInfo_t ri;
  time_string ts1, ts2, emst1, emst2;
  time_string trec; //BES3 record time to data base 
  std::vector<RunInfo_t> RI2;
  while(file)
  {
    unsigned tmprun;
    file >> tmprun;
    if(file.eof()) break;
    ri.run=tmprun;
    file >> ri.Nmh >> ri.lum >> ts1 >> ts2;
    file >> ri.BEPC_Ee >> ri.BEPC_Ep >>  ri.Ip1 >> ri.Ip2 >> ri.Ie1 >> ri.Ie2 >> trec;
    ri.begin_time = ts1.unixtime;
    ri.end_time = ts2.unixtime;
    file  >> ri.e.E >> ri.e.dE >> ri.e.S >> ri.e.dS >> emst1 >> emst2;
    ri.e.begin_time = emst1.unixtime;
    ri.e.end_time = emst2.unixtime;
    file  >> ri.p.E >> ri.p.dE >> ri.p.S >> ri.p.dS >> emst1 >> emst2;
    ri.p.begin_time = emst1.unixtime;
    ri.p.end_time = emst2.unixtime;
    ri.W = cm_energy(ri.e.E,ri.p.E);
    ri.dW=dcm_energy(ri.e.E,ri.e.dE,ri.p.E,ri.p.dE);
    file.ignore(1000,'\n');
    RI2.push_back(ri);
    std::cout.precision(14);
    std::cout << ri.run << " " << ri.lum << " " << ts1.unixtime << " " << ts2.unixtime << " ";
    std::cout << ri.e.E << " " << ri.p.E << std::endl;
  }
  file.close();
  RIsize = RI2.size();
  RI = new RunInfo_t[RIsize];
  for(unsigned i=0;i<RI2.size();++i)
  {
    RI[i]=RI2[i];
  }
  std::cout << "Finish reading file" << std::endl;
  cout << "RI="<<RI << endl;
}

void read_run_info(const char * filename, std::vector<RunInfo_t> & RI)
{
  std::ifstream file(filename);
  if(!file)
  {
    std::cerr << "Unable to open file " << filename << endl;
    return;
  }
  RunInfo_t ri;
  time_string ts1, ts2, emst1, emst2;
  time_string trec; //BES3 record time to data base 
  std::vector<RunInfo_t> RI2;
  while(file)
  {
    unsigned tmprun;
    file >> tmprun;
    if(file.eof()) break;
    ri.run=tmprun;
    file >> ri.Nmh >> ri.lum >> ts1 >> ts2;
    file >> ri.BEPC_Ee >> ri.BEPC_Ep >>  ri.Ip1 >> ri.Ip2 >> ri.Ie1 >> ri.Ie2 >> trec;
    ri.begin_time = ts1.unixtime;
    ri.end_time = ts2.unixtime;
    file  >> ri.e.E >> ri.e.dE >> ri.e.S >> ri.e.dS >> emst1 >> emst2;
    ri.e.begin_time = emst1.unixtime;
    ri.e.end_time = emst2.unixtime;
    cout << emst1.str << " " << emst2.str << endl;
    cout << emst1.unixtime << " " << emst2.unixtime << endl;

    file  >> ri.p.E >> ri.p.dE >> ri.p.S >> ri.p.dS >> emst1 >> emst2;
    ri.p.begin_time = emst1.unixtime;
    ri.p.end_time = emst2.unixtime;
    ri.W = cm_energy(ri.e.E,ri.p.E);
    ri.dW=dcm_energy(ri.e.E,ri.e.dE,ri.p.E,ri.p.dE);
    file.ignore(1000,'\n');
    RI2.push_back(ri);
    std::cout.precision(14);
    std::cout << ri.run << " " << ri.lum << " " << ts1.unixtime << " " << ts2.unixtime << " ";
    std::cout << ri.e.E << " " << ri.p.E << " ";
    std::cout << ri.e.S << " " << ri.p.S << " ";
    std::cout << ri.W << " " << ri.dW << " ";
    std::cout<< std::endl;
  }
  file.close();
  RI.resize(RI2.size());
  for(unsigned i=0;i<RI2.size();++i)
  {
    RI[i]=RI2[i];
  }
  std::cout << "Finish reading file" << std::endl;
  cout << "RI="<<RI.size() << endl;
}
