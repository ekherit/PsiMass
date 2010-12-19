/*
 * =====================================================================================
 *
 *       Filename:  first-scan.C
 *
 *    Description:  Root script for make result of first psi prime scan.
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
  gSystem->CompileMacro("first-scan.cpp","kO");
  make_result();
}
