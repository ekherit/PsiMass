/*
 * =====================================================================================
 *
 *       Filename:  draw.C
 *
 *    Description:  Draw reslt of selection
 *
 *        Version:  1.0
 *        Created:  11/18/2010 10:26:54 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */
{
	TFile f("jpsi.root");
	TTree * chtr = (TTree*)f.Get("chtr");
	TTree * main = (TTree*)f.Get("main");
	chtr->AddFriend(main);


	const char * scut1 = "signal==1";
	const char * bcut1 = "signal==0";
	TCanvas * cE  = new TCanvas;
	chtr->Draw("Etotal", "");
	chtr->Draw("Etotal>>hEs", scut1);
	hEs->SetLineColor(kRed);
	hEs->SetTitle("Total charged track energy deposition");
	hEs->GetXaxis()->SetTitle("E [GeV]");
	hEs->SetLineWidth(2);
	chtr->Draw("Etotal>>hEb",bcut1,"same");
	hEb->SetLineColor(kBlue);
	TLegend * lE = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE->AddEntry(hEs, "Signal (ncharged>2)", "lp");
	lE->AddEntry(hEb, "Bhabha (ncharded=2)", "lp");
	lE->Draw();

	TCanvas * cS = new TCanvas; //sphericity
	chtr->Draw("S", "");
	chtr->Draw("S>>hSs",scut1);
	hSs->SetLineColor(kRed);
	hSs->SetTitle("Sphericity");
	hSs->GetXaxis()->SetTitle("S");
	hSs->SetLineWidth(2);
	chtr->Draw("S>>hSb", bcut1,"same");
	hSb->SetLineColor(kBlue);
	TLegend * lS = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS->AddEntry(hSs, "Signal (ncharged>2)", "lp");
	lS->AddEntry(hSb, "Bhabha (ncharded=2)", "lp");
	lS->Draw();

	//Now apply some cuts to see what happen
	const char * scut2 = "signal==1&&S>0.1";
	const char * bcut2 = "signal==0&&S<0.1";
	TCanvas * cE2  = new TCanvas;
	chtr->Draw("Etotal", "");
	chtr->Draw("Etotal>>hE2s", scut2);
	hE2s->SetLineColor(kRed);
	hE2s->SetTitle("Total charged track energy deposition");
	hE2s->GetXaxis()->SetTitle("E [GeV]");
	hE2s->SetLineWidth(2);
	chtr->Draw("Etotal>>hE2b",bcut2,"same");
	hE2b->SetLineColor(kBlue);
	TLegend * lE2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE2->AddEntry(hE2s, "Signal (ncharged>2&&S>0.1)", "lp");
	lE2->AddEntry(hE2b, "Bhabha (ncharded=2&&S<0.1)", "lp");
	lE2->Draw();

	const char * scut22 = "signal==1&&Etotal>1&&Etotal<2.5";
	const char * bcut22 = "signal==0&&Etotal>2.5";
	TCanvas * cS2 = new TCanvas; //sphericity
	chtr->Draw("S", "");
	chtr->Draw("S>>hS2s",scut22);
	hS2s->SetLineColor(kRed);
	hS2s->SetTitle("Sphericity");
	hS2s->GetXaxis()->SetTitle("S");
	hS2s->SetLineWidth(2);
	chtr->Draw("S>>hS2b", bcut22,"same");
	hS2b->SetLineColor(kBlue);
	TLegend * lS2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS2->AddEntry(hS2s, "Signal (ncharged>2&&Etotal>1&&Etotal<2.5)", "lp");
	lS2->AddEntry(hS2b, "Bhabha (ncharded=2&&Etotal>2.5)", "lp");
	lS2->Draw();
}
