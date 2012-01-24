/*
 * =====================================================================================
 *
 *       Filename:  mhsel.C
 *
 *    Description:  Select multihadron events
 *
 *        Version:  1.0
 *        Created:  24.01.2012 16:39:27
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan B. Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics, Novosibirsk, Russia
 *
 * =====================================================================================
 */
{
  gSystem->CompileMacro("interference.cpp","kO");
  gSystem->CompileMacro("RunInfo.cpp","kO");
  gSystem->CompileMacro("mhsel.cpp","kO");
  draw_energy_vs_time();
  //make_result();
}

