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
	TFile f("proceed.root");
	TTree * tree = (TTree*)f.Get("mhadr");
	TTree * dedx = (TTree*)f.Get("dedx");
	tree->AddFriend(dedx);
	TCanvas * c = new TCanvas;
	c->Divide(2, 2);


	const char * scut1 = "signal==1";
	const char * bcut1 = "signal==0";
	//TCanvas * cE  = new TCanvas;
	c->cd(1);
	tree->Draw("Etotal", "");
	tree->Draw("Etotal>>hEs", scut1);
	hEs->SetLineColor(kRed);
	hEs->SetTitle("Total charged track energy deposition");
	hEs->GetXaxis()->SetTitle("E [GeV]");
	hEs->SetLineWidth(2);
	tree->Draw("Etotal>>hEb",bcut1,"same");
	hEb->SetLineColor(kBlue);
	TLegend * lE = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE->AddEntry(hEs, "Signal (ncharged>2)", "lp");
	lE->AddEntry(hEb, "Bhabha (ncharded=2)", "lp");
	lE->Draw();

	c->cd(3);
	//TCanvas * cS = new TCanvas; //sphericity
	tree->Draw("S", "");
	tree->Draw("S>>hSs",scut1);
	hSs->SetLineColor(kRed);
	hSs->SetTitle("Sphericity");
	hSs->GetXaxis()->SetTitle("S");
	hSs->SetLineWidth(2);
	tree->Draw("S>>hSb", bcut1,"same");
	hSb->SetLineColor(kBlue);
	TLegend * lS = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS->AddEntry(hSs, "Signal (ncharged>2)", "lp");
	lS->AddEntry(hSb, "Bhabha (ncharded=2)", "lp");
	lS->Draw();

	//Now apply some cuts to see what happen
	const char * scut2 = "signal==1&&S>0.1";
	const char * bcut2 = "signal==0&&S<0.1";
	//TCanvas * cE2  = new TCanvas;
	c->cd(2);
	tree->Draw("Etotal", "");
	tree->Draw("Etotal>>hE2s", scut2);
	hE2s->SetLineColor(kRed);
	hE2s->SetTitle("Total charged track energy deposition");
	hE2s->GetXaxis()->SetTitle("E [GeV]");
	hE2s->SetLineWidth(2);
	tree->Draw("Etotal>>hE2b",bcut2,"same");
	hE2b->SetLineColor(kBlue);
	TLegend * lE2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE2->AddEntry(hE2s, "Signal (ncharged>2&&S>0.1)", "lp");
	lE2->AddEntry(hE2b, "Bhabha (ncharded=2&&S<0.1)", "lp");
	lE2->Draw();

	const char * scut22 = "signal==1&&Etotal>1&&Etotal<2.5";
	const char * bcut22 = "signal==0&&Etotal>2.5";
	//TCanvas * cS2 = new TCanvas; //sphericity
	c->cd(4);
	tree->Draw("S", "");
	tree->Draw("S>>hS2s",scut22);
	hS2s->SetLineColor(kRed);
	hS2s->SetTitle("Sphericity");
	hS2s->GetXaxis()->SetTitle("S");
	hS2s->SetLineWidth(2);
	tree->Draw("S>>hS2b", bcut22,"same");
	hS2b->SetLineColor(kBlue);
	TLegend * lS2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS2->AddEntry(hS2s, "Signal (ncharged>2&&Etotal>1&&Etotal<2.5)", "lp");
	lS2->AddEntry(hS2b, "Bhabha (ncharded=2&&Etotal>2.5)", "lp");
	lS2->Draw();

	TCanvas * ccos = new TCanvas;
	tree->Draw("coshp>>hCs", "signal==1&&coshp>-5");
	hCs->SetLineColor(kRed);
	hCs->SetLineWidth(2);
	hCs->GetXaxis()->SetTitle("cos #phi");
	tree->Draw("coshp>>hCb", "signal==0&&coshp>-5", "same");
	TLegend * lC = new TLegend(0.6, 0.8, 1.0, 1.0);
	lC->AddEntry(hCs, "Signal (ncharged>2)", "lp");
	lC->AddEntry(hCb, "Bhabha (ncharded=2)", "lp");
	lC->Draw();
}
