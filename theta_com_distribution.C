//This script is used to plot the theta com (center of mass) distributions for ee (moller).
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void theta_com_distribution(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	TGaxis::SetMaxDigits(3);
	
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";//directory where the rootfiles exist
	TString file[] = {"main_sm_moller/main_sm_moller*.root"};
	int nfile = sizeof(file)/sizeof(*file);
	
	TH1F* h_thcom[nfile];//histograms for total th_thcom distributions
	TH1F* h_thcom_c[nfile];//histograms for closed sector th_thcom distributions
	TH1F* h_thcom_o[nfile];//histograms for open sector th_thcom distribution
	TH1F* h_thcom_t[nfile];//histograms for transition sector th_thcom distributions
	int color[] = {1,1,1};//colors for total th_thcom distributions
	int coloro[] = {2,2,2};//colors for open sector th_thcom distributions
	int colorc[] = {4,4,4};//colors for closed sector th_thcom distributions
	int colort[] = {3,3,3};//colors for transition sector th_thcom distributions
	TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph)
	TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sectors
	TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sectors
	TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sectors
	TString weight;//define conditional weight inside the loop below

	TString sim[] = {"ee","ep-el","ep-inel"};

	for(int ifile=0;ifile<nfile;ifile++){
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 500;
	double x_min = 0;
	double x_max = 180;
	double bin_width = (x_max-x_min)/nbin;
	h_thcom[ifile] = new TH1F(Form("h_thcom[%d]",ifile),Form("%s #theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_thcom_o[ifile] = new TH1F(Form("h_thcom_o[%d]",ifile),Form("%s #theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_thcom_c[ifile] = new TH1F(Form("h_thcom_c[%d]",ifile),Form("%s #theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_thcom_t[ifile] = new TH1F(Form("h_thcom_t[%d]",ifile),Form("%s #theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_thcom[ifile]->SetLineColor(color[ifile]);
	h_thcom_o[ifile]->SetLineColor(coloro[ifile]);
	h_thcom_c[ifile]->SetLineColor(colorc[ifile]);
	h_thcom_t[ifile]->SetLineColor(colort[ifile]);
	if(ifile==0)
	weight = "1e-9*rate/2.0";//moller rate needs to be divided by 2.0
	else
	weight = "1e-9*rate";
	T->Draw(Form("ev.thcom/deg>>h_thcom[%d]",ifile),Form("%s*(%s)",weight.Data(),ring5Cut.Data()),"goff");
	T->Draw(Form("ev.thcom/deg>>h_thcom_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("ev.thcom/deg>>h_thcom_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("ev.thcom/deg>>h_thcom_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	}
	
	TCanvas* c1[nfile];
	TLegend* leg = new TLegend(0.70,0.65,0.9,0.9);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	TLegendEntry* leg1[4];
	for(int ifile=0;ifile<nfile;ifile++){
	c1[ifile] = new TCanvas(Form("c1[%d]",ifile),Form("%s",sim[ifile].Data()));
	h_thcom[ifile]->Draw();
	h_thcom_o[ifile]->Draw("same");
	h_thcom_c[ifile]->Draw("same");
	h_thcom_t[ifile]->Draw("same");

	leg1[0] = leg->AddEntry(h_thcom[ifile],Form("total %s",sim[ifile].Data()),"le");
	leg1[1] = leg->AddEntry(h_thcom_o[ifile],"open","le");
	leg1[2] = leg->AddEntry(h_thcom_c[ifile],"close","le");
	leg1[3] = leg->AddEntry(h_thcom_t[ifile],"transition","le");
	leg->SetTextSize(0.05);
	leg1[0]->SetTextColor(color[ifile]);
	leg1[1]->SetTextColor(coloro[ifile]);
	leg1[2]->SetTextColor(colorc[ifile]);
	leg1[3]->SetTextColor(colort[ifile]);
	leg->Draw();
	c1[ifile]->SaveAs(Form("../temp/plot%d_thcom.pdf",ifile));
	}
	TCanvas* c2[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	c2[ifile] = new TCanvas(Form("c2[%d]",ifile),Form("%s logy",sim[ifile].Data()));
	gPad->SetLogy();
	h_thcom[ifile]->Draw();
	h_thcom_o[ifile]->Draw("same");
	h_thcom_c[ifile]->Draw("same");
	h_thcom_t[ifile]->Draw("same");

	leg->Draw();
	c2[ifile]->SaveAs(Form("../temp/plot%d_thcom_logy.pdf",ifile));
	}
	TCanvas* c3[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	c3[ifile] = new TCanvas(Form("c3[%d]",ifile),Form("%s hist",sim[ifile].Data()));
	h_thcom[ifile]->Draw("hist");
	h_thcom_o[ifile]->Draw("hist same");
	h_thcom_c[ifile]->Draw("hist same");
	h_thcom_t[ifile]->Draw("hist same");
	leg->Draw();
	c3[ifile]->SaveAs(Form("../temp/plot%d_thcom_hist.pdf",ifile));
	}
	TCanvas* c4[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	c4[ifile] = new TCanvas(Form("c4[%d]",ifile),Form("%s logy hist",sim[ifile].Data()));
	gPad->SetLogy();
	h_thcom[ifile]->Draw("hist");
	h_thcom_o[ifile]->Draw("hist same");
	h_thcom_c[ifile]->Draw("hist same");
	h_thcom_t[ifile]->Draw("hist same");
	leg->Draw();
	c4[ifile]->SaveAs(Form("../temp/plot%d_thcom_logy_hist.pdf",ifile));
	}
	gSystem->Exec(Form("pdfunite ../temp/plot*_thcom*.pdf ../plots/thcom_ee.pdf"));
	gSystem->Exec(Form("rm -rf ../temp/plot*_thcom*.pdf"));
}
