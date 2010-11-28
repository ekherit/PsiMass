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
#include <TTree.h>
void draw(TTree * tree,  const char * title)
{
	TCanvas * c = new TCanvas;
	c->SetTitle(title);
	c->Divide(2, 2);
	c->SetBorderMode(0);
	c->SetFillColor(0);


	const char * scut1 = "signal==1";
	const char * bcut1 = "signal==0";
	//TCanvas * cE  = new TCanvas;
	c->cd(1);
	gPad->SetFillColor(0);
	tree->Draw("Etotal", "");
	tree->Draw("Etotal>>hEs", scut1);
	hEs->SetLineColor(kRed);
	hEs->SetTitle("Total charged track energy deposition");
	hEs->GetXaxis()->SetTitle("E [GeV]");
	hEs->SetLineWidth(2);
	tree->Draw("Etotal>>hEb",bcut1,"same");
	hEb->SetLineColor(kBlue);
	TLegend * lE = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE->AddEntry(hEs, "N_{q}>2", "lp");
	lE->AddEntry(hEb, "N_{q}=2", "lp");
	lE->Draw();

	c->cd(3);
	gPad->SetFillColor(0);
	gPad->SetLogy();
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
	lS->AddEntry(hSs, "N_{q}>2", "lp");
	lS->AddEntry(hSb, "N_{q}=2", "lp");
	lS->Draw();

	//Now apply some cuts to see what happen
	const char * scut2 = "signal==1&&S>0.1";
	const char * bcut2 = "signal==0&&S<0.1";
	//TCanvas * cE2  = new TCanvas;
	c->cd(2);
	gPad->SetFillColor(0);
	tree->Draw("Etotal", "");
	tree->Draw("Etotal>>hE2s", scut2);
	hE2s->SetLineColor(kRed);
	hE2s->SetTitle("Total charged track energy deposition");
	hE2s->GetXaxis()->SetTitle("E [GeV]");
	hE2s->SetLineWidth(2);
	tree->Draw("Etotal>>hE2b",bcut2,"same");
	hE2b->SetLineColor(kBlue);
	TLegend * lE2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE2->AddEntry(hE2s, "N_{q}>2, S>0.1", "lp");
	lE2->AddEntry(hE2b, "N_{q}=2, S<0.1", "lp");
	lE2->Draw();

	const char * scut22 = "signal==1&&Etotal>0.5&&Etotal<2.5";
	const char * bcut22 = "signal==0&&Etotal>2.5";
	//TCanvas * cS2 = new TCanvas; //sphericity
	c->cd(4);
	gPad->SetFillColor(0);
	gPad->SetLogy();
	tree->Draw("S", "");
	tree->Draw("S>>hS2s",scut22);
	hS2s->SetLineColor(kRed);
	hS2s->SetTitle("Sphericity");
	hS2s->GetXaxis()->SetTitle("S");
	hS2s->SetLineWidth(2);
	tree->Draw("S>>hS2b", bcut22,"same");
	hS2b->SetLineColor(kBlue);
	TLegend * lS2 = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS2->AddEntry(hS2s, "N_{q}>2, E_{tot} #in (0.5, 2.5)GeV", "lp");
	lS2->AddEntry(hS2b, "N_{q}=2, E_{tot}>2.5 GeV", "lp");
	lS2->Draw();

	TCanvas * ccos = new TCanvas;
	ccos->SetTitle(title);
	tree->Draw("coshp>>hCs", "signal==1&&coshp>-5");
	hCs->SetLineColor(kRed);
	hCs->SetLineWidth(2);
	hCs->GetXaxis()->SetTitle("cos #phi");
	tree->Draw("coshp>>hCb", "signal==0&&coshp>-5", "same");
	TLegend * lC = new TLegend(0.6, 0.8, 1.0, 1.0);
	lC->AddEntry(hCs, "Signal (ncharged>2)", "lp");
	lC->AddEntry(hCb, "Bhabha (ncharded=2)", "lp");
	lC->Draw();

	TCanvas * c  = new TCanvas;
	c->SetTitle(title);
	tree->Draw("Etotal",  "fabs(chimu[0])<1.5");

	//Calculate result.
	//const char * main_signal_cut = "nchtrk>2&&Etotal>0.5&&Etotal<2.5&&S>0.1";
	const char * main_signal_cut = "nchtrk>2&&S>0.1&&Etotal>0.5&&Etotal<2.5";
	const char * main_bhabha_cut = "nchtrk==2&&S<0.1&&Etotal>2.5";
	TCanvas * rc = new TCanvas;
	rc->SetTitle(title);
	tree->Draw("Etotal>>Rsh(200, 0, 4.0)", main_signal_cut);
	unsigned Nsignal = tree->GetSelectedRows();
	Rsh->SetLineColor(kRed);
	Rsh->SetLineWidth(2);
	tree->Draw("Etotal>>Rbh(200, 0, 4.0)", main_bhabha_cut, "same");

	TCanvas * cnch = new TCanvas;
	cnch->SetTitle(title);
	cnch->SetLogy();
	cnch->SetFillColor(0);
	cnch->SetBorderMode(0);
	tree->Draw("nchtrk>>hNs(20, 0, 20)", "Etotal>0.5&&Etotal<2.5&&S>0.1");
	hNs->SetLineColor(kRed);
	hNs->SetLineWidth(2);
	hNs->SetTitle(title);
	hNs->GetXaxis()->SetTitle("number of charged tracks");
	tree->Draw("nchtrk>>hNb(20, 0, 20)", "Etotal>2.5&&S<0.1", "same");
	double maxb = hNb->GetMaximum();
	cout <<"Maximum = "<<maxb << endl;
	hNs->SetMaximum(maxb);
	cnch->Update();
	TLegend * lnch = new TLegend(0.6, 0.8, 1.0, 1.0);
	lnch->AddEntry(hNs, "0.5 < E_{tot} < 2.5 GeV, S>0.1", "lp");
	lnch->AddEntry(hNb, " E_{tot} >2.5 GeV, S<0.1", "lp");
	lnch->Draw();


	unsigned Nbhabha = tree->GetSelectedRows();
	double R=(double)Nsignal/(double)Nbhabha;
	cout << title << " : " << " signal = " << Nsignal << ",  bhabha = " << Nbhabha << ",  signal/bhabha = " << R << endl;
	
}


void draw_charged_cut(TTree * tree,  const char * title)
{
	TCanvas * c = new TCanvas(title, title, 640, 480/2*3*3/2);
	c->SetTitle(title);
	c->Divide(1, 2);
	c->SetBorderMode(0);
	c->SetFillColor(0);

	const char * scut1 = "signal==1";
	const char * bcut1 = "signal==0";
	c->cd(1);
	gPad->SetFillColor(0);
	gPad->SetBorderMode(0);
	tree->Draw("Etotal", "");
	tree->Draw("Etotal>>hEs", scut1);
	hEs->SetLineColor(kRed);
	hEs->SetTitle("Total charged track energy deposition");
	hEs->GetXaxis()->SetTitle("E [GeV]");
	hEs->SetLineWidth(3);
	double hEsmax = hEs->GetMaximum();
	tree->Draw("Etotal>>hEb",bcut1,"same");
	hEb->SetLineColor(kBlue);
	hEb->SetLineWidth(3);
	double hEbmax = hEb->GetMaximum();
	hEs->SetMaximum(max(hEsmax, hEbmax));
	TLegend * lE = new TLegend(0.6, 0.8, 1.0, 1.0);
	lE->AddEntry(hEs, "N_{q}>2", "lp");
	lE->AddEntry(hEb, "N_{q}=2", "lp");
	lE->Draw();

	c->cd(2);
	gPad->SetFillColor(0);
	gPad->SetBorderMode(0);
	gPad->SetLogy();
	tree->Draw("S", "");
	tree->Draw("S>>hSs",scut1);
	hSs->SetLineColor(kRed);
	hSs->SetTitle("Sphericity");
	hSs->GetXaxis()->SetTitle("S");
	hSs->SetLineWidth(3);
	double hSsmax = hSs->GetMaximum();
	tree->Draw("S>>hSb", bcut1,"same");
	hSb->SetLineColor(kBlue);
	hSb->SetLineWidth(3);
	double hSbmax = hSb->GetMaximum();
	hSs->SetMaximum(max(hSbmax, hSsmax));


	TLegend * lS = new TLegend(0.6, 0.8, 1.0, 1.0);
	lS->AddEntry(hSs, "N_{q}>2", "lp");
	lS->AddEntry(hSb, "N_{q}=2", "lp");
	lS->Draw();
}


void draw_angle_plot(TTree * tree,  const char * title)
{
	TCanvas * c = new TCanvas;
	c->SetTitle(title);
	c->SetBorderMode(0);
	c->SetFillColor(0); tree->Draw("theta:Etotal>>hsang", "fabs(theta)<10","lego2");
	hsang->SetTitle("Scattering angle");
	hsang->GetXaxis()->SetTitle("E [GeV]");
	hsang->GetYaxis()->SetTitle("#theta");
}
void draw_angle_plot2(TTree * tree,  const char * title)
{
	TCanvas * c = new TCanvas;
	c->SetLogz();
	c->SetTitle(title);
	c->SetBorderMode(0);
	c->SetFillColor(0); tree->Draw("theta:S>>hsangS", "fabs(theta)<10","lego2");
	hsangS->SetTitle("Scattering angle");
	hsangS->GetXaxis()->SetTitle("S");
	hsangS->GetYaxis()->SetTitle("#theta");
}

void draw_ES_plot(TTree * tree,  const char * title)
{
	TCanvas * c = new TCanvas;
	c->SetLogz();
	c->SetTitle(title);
	c->SetBorderMode(0);
	c->SetFillColor(0); 
	tree->Draw("Etotal:S>>hES", "fabs(E)<10&&fabs(S)<10","lego2");
	hES->SetTitle("Total energy and sphericity");
	hES->GetXaxis()->SetTitle("S");
	hES->GetYaxis()->SetTitle("E");
}


void draw_dalits(TTree * tree, const char * title)
{
	new TCanvas;
	tree->SetAlias("M01","(E[0]-E[1])**2-(px[0]-px[1])**2-(py[0]-py[1])**2-(pz[0]-pz[1])**2");
	tree->SetAlias("M02","(E[0]-E[2])**2-(px[0]-px[2])**2-(py[0]-py[2])**2-(pz[0]-pz[2])**2");
	tree->SetAlias("M12","(E[1]-E[2])**2-(px[1]-px[2])**2-(py[1]-py[2])**2-(pz[1]-pz[2])**2");
	tree->Draw("E[2]", "nchtrk==3&&fabs(M01)<100","lego2");
}

void draw(const char * file, const char * title)
{
}

void draw(void)
{
	//signal
	TFile * psip_file = new TFile("psip.root");
	TTree * psip = (TTree*)psip_file->Get("mhadr");
	psip->AddFriend("dedx");
	draw(psip, "psi prime data");
	draw_charged_cut(psip, "#Psi^{#prime} data");
	draw_angle_plot(psip, "psip");
	draw_angle_plot2(psip, "psip");
	draw_ES_plot(psip, "psip");
  //draw_dalits(psip,"psip");


	////monte carlo
	TFile * mc_file = new TFile("mcpsip.root");
	TTree * mc = (TTree*)mc_file->Get("mhadr");
	mc->AddFriend("dedx");

	draw_charged_cut(mc,  "psi prime Monte Carlo");
	draw_angle_plot(mc, "psip MC");
	draw_angle_plot2(mc, "psip MC");

	//continuum
	TFile * cont_file = new TFile( "con365.root");
	TTree * cont = (TTree*)cont_file->Get("mhadr");
	cont->AddFriend("dedx");
	draw_charged_cut(cont,  "Continuum 3.65 GeV");
	draw_angle_plot(cont, "Continum");
	
}
