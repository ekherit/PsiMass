/*
 * =====================================================================================
 *
 *       Filename:  extract.C
 *
 *    Description:  Extract the ee data from files
 *
 *        Version:  1.0
 *        Created:  05/25/2011 11:23:02 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ivan Nikolaev (), I.B.Nikolaev@inp.nsk.su
 *        Company:  Budker Institute of Nuclear Physics
 *
 * =====================================================================================
 */

#include <iomanip>
#include "selection.h"
#include "interference.h"
#include "utils.h"

const unsigned JOB_NUMBER=4;
//string file_prefix="bbyg_gg_";
//string file_prefix="bhwide_";
//string file_prefix="bbyg_ee_";
string file_prefix="bbyg_geom_ee_";

void FixParNoInterference(TF1 * f)
{
	f->FixParameter(1, 0); 
	f->FixParameter(2, 0); 
	f->FixParameter(3, 1); //From PDG table
}

void draw_bhabha(void)
{
	unsigned Elist_size=15;
	double Elist[15]={1839.0, 
                  1841.4, 
                  1843.2, 
                  1844.1, 
                  1848.4, 
                  1838.2, 
                  1840.8, 
                  1841.8, 
                  1843.0, 
                  1843.5, 
                  1845.1, 
                  1846.5, 
									1842.0, 1842.5,  1842.75}; //this is extra points.
	unsigned N0[1024];
	string Estr[1024];
	//energy in GeV
	double E[1024];
	//initial cross section and error
	double CR0[1024], CR0err[1024];

	//extracting the CrossSection
	cout << "Reading cross section" << endl;
	TGraphErrors * cr0g = new TGraphErrors;
	for(unsigned i=0; i<Elist_size;++i)
	{
		double sum=0;
		double sum2=0;
		double wsum=0;
		unsigned job=1;
		N0[i]=0;
		TGraphErrors graph(0);
		for(unsigned job=1;job<JOB_NUMBER+1;++job)
		{
			char name[1024];
			sprintf(name, "%s%4.1f_%d", file_prefix.c_str(), Elist[i], job);
			if(Elist[i]==1842.75)
			sprintf(name, "%s%4.2f_%d", file_prefix.c_str(), Elist[i], job);
			char cros_name[1024];
			sprintf(cros_name, "%s/CrossSection.txt",name);
			double cross_section;
			double cross_section_error;
			double w;//weight
			ifstream f(cros_name);
			if(!f) 
			{
				char cros2_name[1024];
				sprintf(cros2_name, "%s/fort.16",name);
				ifstream f2(cros2_name);
				if(!f2)
				{
					cout << "Unable to find cross section" << endl;
					continue;
				}
				string line;
				while(getline(f2, line))
				{
					size_t found = line.find("Accepted total          NEVGEN");
					if(found!=string::npos) 
					{
						string tmp;
						istringstream is(line.c_str());
						unsigned event_number;
						is >> tmp >> event_number;
						N0[i]+=event_number;
					}
					found = line.find("Xsec M.C. [pb]");
					if(found!=string::npos) 
					{
						string tmp;
						istringstream is(line.c_str());
						is >> tmp >> cross_section;
						cross_section/=1e3; //cross section in bhwide fort.16 file in pb.
						getline(f2, line); //read the relative error
						istringstream is2(line.c_str());
						is2 >> tmp >> cross_section_error;
						cross_section_error*=cross_section;
					}
				}
			}
			else 
			{
				string line;
				while(getline(f, line))
				{
					size_t found = line.find("EVENT NUMBER");
					if(found!=string::npos) 
					{
						istringstream is(line.c_str());
						unsigned event_number;
						string tmp;
						is >> tmp >> tmp >> event_number;
						N0[i]+=event_number;
					}
					found = line.find("UNWEIGHTED CROSS SECTION");
					if(found!=string::npos) 
					{
						f >> cross_section >> line >> cross_section_error;
					}
				}
			}
			graph.SetPoint(job-1, job-1,  cross_section);
			graph.SetPointError(job-1, 0,  cross_section_error);
			w = 1./sq(cross_section_error);
			sum+=cross_section*w;
			sum2+=sq(cross_section)*w;
			wsum+=w;
			cout << Elist[i] <<  " " << job << ": " << cross_section << "+-" << cross_section_error << endl;
		}
		double cross_section = sum/wsum;
		double cross_section_error = sqrt(sum2/wsum - sq(cross_section));
		graph.Fit("pol0", "Qgoff");
		TF1 * fun = graph.GetFunction("pol0");
		cross_section = fun->GetParameter(0);
		cross_section_error = fun->GetParError(0);
		CR0[i] = cross_section;
		CR0err[i] = cross_section_error;
		cout << name << ": " << cross_section<< "+-"<< cross_section_error << " nb,   N0="<< N0[i] <<  endl;
		cr0g->SetPoint(i, Elist[i]*2-MPDG, CR0[i]);
		cr0g->SetPointError(i, 0, CR0err[i]);
	}
	gStyle->SetOptFit();
	TF1 * fun_cr0 = new TF1("fun_cr0",&sigma, -10, 10, 4);
	fun_cr0->SetParameter(0, 900);//nb
	fun_cr0->SetParameter(1, 30);//nb
	fun_cr0->SetParameter(2, 30);//nb
	fun_cr0->SetParameter(3, 0.3);//MeV
	//fun_cr0->FixParameter(3, 0.304); //From PDG table
	fun_cr0->SetParName(0, "QED");
	fun_cr0->SetParName(1, "INT");
	fun_cr0->SetParName(2, "RES");
	fun_cr0->SetParName(3, "Gamma");
	//FixParNoInterference(fun_cr0);
	TCanvas * cr0c = new TCanvas("cr0c", "Total crossection from generator");
	cr0g->SetMarkerStyle(21);
	cr0g->Draw("ap");
	cr0g->GetXaxis()->SetTitle("W-M_{#psi},  MeV");
	cr0g->GetYaxis()->SetTitle("#sigma_{ee},  nb");
	TF1 * cr0_sfun = new TF1("sfun_cr0",&sigma_spread,-10, 10, 5 );
	double par_cr0[5];
	par_cr0[0]=1; //some spread.
	par_cr0[1]=900;
	par_cr0[2]=0;
	par_cr0[3]=0;
	par_cr0[4]=0.304; //psip width
	cr0_sfun->SetParameters(par_cr0);
	cr0_sfun->SetParLimits(0,0.1, 3);
	//cr0_sfun->SetParLimits(2,-100, 100);
	//cr0_sfun->SetParLimits(3,-100, 100);
	cr0_sfun->FixParameter(4, 0.304); //From PDG table
	cr0g->Fit("fun_cr0");
	//cr0g->Fit("sfun_cr0");
	return;

	TCut mh_cut,  ee_cut,  gg_cut;
	set_selection(7, mh_cut,  ee_cut,  gg_cut);
	
	TGraphErrors * sigma_g = new TGraphErrors;
	TGraphErrors * ggsigma_g = new TGraphErrors;
	TGraphErrors * ngg_g = new TGraphErrors; //Number of selected gg events
	TGraphErrors * nmh_g = new TGraphErrors; //Number of selected gg events
	//open root file and procceed number of events
	cout << "Calculating number of events" << endl;
	cout << setw(12) << "Energy, MeV" << setw(10) << "Nmh" << setw(10)  
		<< "Nee" << setw(10) << "Ngg" 
		<< setw(15) << "See, nb" << setw(15) << "Sgg, nb" << endl;
	for(unsigned i=0; i< Elist_size; i++)
	{
		TChain * mdc = new TChain("mdc", "mdc");
		TChain * gg = new TChain("gg", "gg");
		char estr[1024];
		sprintf(estr, "%f*1", Elist[i]/1.0e3);
		mdc->SetAlias("Eb", estr);
		gg->SetAlias("Eb", estr);
		for(unsigned job=1; job<JOB_NUMBER+1; ++job)
		{
			char dir_name[1024];
			sprintf(dir_name, "%s%4.1f_%d", file_prefix.c_str(), Elist[i], job);
			if(Elist[i]==1842.75)
			sprintf(dir_name, "%s%4.2f_%d", file_prefix.c_str(), Elist[i], job);
			char file_name[1024];
			sprintf(file_name, "%s/%s.root", dir_name, dir_name);
			//cout << dir_name << " " << file_name<< endl;
			int rcmdc = mdc->AddFile(file_name);
			int rcgg = gg->AddFile(file_name);
			//if(rcgg==0 || rcmdc==0) continue;
		}
		if(mdc->GetEntries()==0) continue;
		mdc->Draw("ntrack", ee_cut, "goff");
		double Nee=mdc->GetSelectedRows();
		double dNee = sqrt(Nee*(1.-Nee/N0[i]));
		double sigma_ee = Nee/N0[i]*CR0[i];
		double dsigma_ee = sqrt(sq(dNee*CR0[i]/N0[i]) + sq(Nee/N0[i]*CR0err[i]));
		sigma_g->SetPoint(i, Elist[i]*2-MPDG, sigma_ee);
		sigma_g->SetPointError(i, 0, dsigma_ee);

		mdc->Draw("ntrack", mh_cut, "goff");
		unsigned Nmh=mdc->GetSelectedRows();
		double dNmh = sqrt(Nmh*(1.-Nmh/N0[i]));
		nmh_g->SetPoint(i, Elist[i]*2-MPDG,Nmh);
		nmh_g->SetPointError(i, 0,dNmh);

		gg->Draw("ntrack", gg_cut, "goff");
		double Ngg=gg->GetSelectedRows();
		double dNgg = sqrt(Ngg*(1.-Ngg/N0[i])); 
		ngg_g->SetPoint(i, Elist[i]*2-MPDG, Ngg); 
		ngg_g->SetPointError(i, 0, dNgg); 
		double sigma_gg = Ngg/N0[i]*CR0[i];
		double dsigma_gg = sqrt(sq(dNgg*CR0[i]/N0[i]) + sq(Ngg/N0[i]*CR0err[i]));
		ggsigma_g->SetPoint(i, Elist[i]*2-MPDG, sigma_gg);
		ggsigma_g->SetPointError(i, 0, dsigma_gg);
		cout << setw(12) << Elist[i] << setw(10) << Nmh << setw(10)  << Nee << setw(10) << Ngg << "+-"<< dNgg
			<< setw(15) << sigma_ee << setw(15) << sigma_gg << " +- " << dsigma_gg << endl;
	}
	TCanvas * sigmac = new TCanvas("sigma_ee_canvas", "#sigma_{ee}^{vis}");
	sigma_g->SetMarkerStyle(21);
	sigma_g->Draw("ap");
	sigma_g->GetXaxis()->SetTitle("W-M_{#psi},  MeV");
	sigma_g->GetYaxis()->SetTitle("#sigma_{vis},  nb");
	//char fit_formula[1024];
	//sprintf(fit_formula, "[0]*pow(%f/(%f+x), 2)+2*[0]*[1]/sqrt(x**2+[2]**2)*cos()", MPDG, MPDG);
	////TF1 * fun = new TF1("fun", "[0]*pow(3686./x, 2)");
	//TF1 * fun = new TF1("fun", fit_formula);
	fun->SetParameter(0, 120);
	TF1 * fun_sigma = new TF1("fun_sigma",&sigma, -10, 10, 4);
	fun_sigma->SetParName(0, "QED");
	fun_sigma->SetParName(1, "INT");
	fun_sigma->SetParName(2, "RES");
	fun_sigma->SetParName(3, "Gamma");
	fun_sigma->SetParameter(0, 100);//nb
	fun_sigma->SetParameter(1, 10);//nb
	fun_sigma->SetParameter(2, 10);//nb
	fun_sigma->SetParameter(3, 0.1);//MeV
	fun_sigma->FixParameter(3, 0.304);
	fun_sigma->SetNpx(1000);
	//FixParNoInterference(fun_sigma);
	sigma_g->Fit("fun_sigma", "E");

	TF1 * sfun = new TF1("sfun",&sigma_spread,-10, 10, 5 );
	double par[5];
	par[0]=1.58; //Beam spread
	par[1]=fun_sigma->GetParameter(0);//QED
	double qed = fun_sigma->GetParameter(0);
	par[1] = qed;
	par[2]=fun_sigma->GetParameter(1);//INT
	par[3]=fun_sigma->GetParameter(2);//RES
	par[4]=fun_sigma->GetParameter(3);//GAMMA
	sfun->SetParameters(par);
	sfun->Draw("same");
	//sfun->Draw("");
	cout << "Parameters of correction is for spread " << par[0] <<  " MeV" <<  endl;
	cout << "QED: " << fun_sigma->GetParameter(0) << endl;
	cout << "INT: " << fun_sigma->GetParameter(1) << endl;
	cout << "RES: " << fun_sigma->GetParameter(2) << endl;
	cout << "GAM: " << fun_sigma->GetParameter(3) << endl;
	/* number of gamma gamma and multihadronic events 
	TCanvas * nggc = new TCanvas("ngg_canvas", "Number of gg events");
	ngg_g->SetMarkerStyle(21);
	ngg_g->Draw("ap");
	ngg_g->GetXaxis()->SetTitle("W-M_{#psi},  MeV");
	ngg_g->GetYaxis()->SetTitle("number of #gamma#gamma events");
	TCanvas * nmhc = new TCanvas("nmh_canvas", "Number of multihadronic events");
	nmh_g->SetMarkerStyle(21);
	nmh_g->Draw("ap");
	nmh_g->GetXaxis()->SetTitle("W-M_{#psi},  MeV");
	nmh_g->GetYaxis()->SetTitle("number of multihadronic  events"); */
	TCanvas * ggsigmac = new TCanvas("sigma_gg_canvas", "#sigma_{gg}^{vis}");
	ggsigma_g->SetMarkerStyle(21);
	ggsigma_g->Draw("ap");
	ggsigma_g->GetXaxis()->SetTitle("W-M_{#psi},  MeV");
	ggsigma_g->GetYaxis()->SetTitle("#sigma_{vis},  nb");
	TF1 * fun_sigma = new TF1("fun_ggsigma",&sigma, -10, 10, 4);
	fun_ggsigma->SetParameter(0, 100);//nb
	fun_ggsigma->SetParName(0, "QED");
	fun_ggsigma->SetParName(1, "INT");
	fun_ggsigma->SetParName(2, "RES");
	fun_ggsigma->SetParName(3, "Gamma");
	fun_ggsigma->SetParameter(1, 10);//nb
	fun_ggsigma->SetParameter(2, 0);//nb
	fun_ggsigma->SetParameter(3, 0);//MeV
	//No resonance contribution here
	fun_ggsigma->FixParameter(1, 0); 
	fun_ggsigma->FixParameter(2, 0); 
	fun_ggsigma->FixParameter(3, 0); //From PDG table
	ggsigma_g->Fit("fun_ggsigma");
};
