/*
 * =====================================================================================
 *
 *       Filename:  RunInfo.C 
 *
 *    Description:  Test for function working with RunInfo.
 *
 *        Version:  1.0
 *        Created:  06.12.2010 18:47:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
{
#include "RunInfo.h"
//#include <vector>
  gSystem->CompileMacro("RunInfo.cpp","kO");
  //std::vector<RunInfo_t> RI;
  RunInfo_t * RI;
  unsigned RISIZE;
  read_run_info("psip-2011-run-info.txt", RI, RISIZE);
  cout << "test" << endl;
}
