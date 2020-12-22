//This script is used to plot the theta lab distribution for ee (moller).
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void theta_lab_distribution(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	TGaxis::SetMaxDigits(3);
	
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";//directory where the rootfiles exist
	TString file = "main_sm_moller/main_sm_moller_*.root";
	
	TH1F* h_thlab;//histograms for total th_thlab distribution
	TH1F* h_thlab_c;//histograms for closed sector th_thlab distribution
	TH1F* h_thlab_o;//histograms for open sector th_thlab distribution
	TH1F* h_thlab_t;//histograms for transition sector th_thlab distribution

	TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph)
	TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sectors
	TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sectors
	TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sectors
	TString weight;//define conditional weight inside the loop below

	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file.Data()));

	int nbin = 500;
	double x_min = 0;
	double x_max = 5;
	double bin_width = (x_max-x_min)/nbin;
	h_thlab = new TH1F("h_thlab",Form("#theta_{lab} distribution (rate weighted);#theta_{lab} (deg);rate (GHz)/%.2fdeg",bin_width),nbin,x_min,x_max);
	h_thlab_o = new TH1F("h_thlab_o",Form("#theta_{lab} distribution (rate weighted);#theta_{lab} (deg);rate (GHz)/%.2fdeg",bin_width),nbin,x_min,x_max);
	h_thlab_c = new TH1F("h_thlab_c",Form("#theta_{lab} distribution (rate weighted);#theta_{lab} (deg);rate (GHz)/%.2fdeg",bin_width),nbin,x_min,x_max);
	h_thlab_t = new TH1F("h_thlab_t",Form("#theta_{lab} distribution (rate weighted);#theta_{lab} (deg);rate (GHz)/%.2fdeg",bin_width),nbin,x_min,x_max);
	h_thlab->SetLineColor(1);
	h_thlab_o->SetLineColor(2);
	h_thlab_c->SetLineColor(4);
	h_thlab_t->SetLineColor(3);
	weight = "1e-9*rate/2.0";//moller rate needs to be divided by 2.0
	T->Draw("part.th/deg>>h_thlab",Form("%s*(%s)",weight.Data(),ring5Cut.Data()),"goff");
	T->Draw("part.th/deg>>h_thlab_o",Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw("part.th/deg>>h_thlab_c",Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw("part.th/deg>>h_thlab_t",Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	
	TLegend* leg = new TLegend(0.70,0.65,0.9,0.9);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetFillColor(0);
	TLegendEntry* leg1[4];
	TCanvas* c1 = new TCanvas("c1","ee");
	h_thlab->Draw();
	h_thlab_o->Draw("same");
	h_thlab_c->Draw("same");
	h_thlab_t->Draw("same");

	leg1[0] = leg->AddEntry(h_thlab,"total ee","le");
	leg1[1] = leg->AddEntry(h_thlab_o,"open","le");
	leg1[2] = leg->AddEntry(h_thlab_c,"close","le");
	leg1[3] = leg->AddEntry(h_thlab_t,"transition","le");
	leg->SetTextSize(0.05);
	leg1[0]->SetTextColor(1);
	leg1[1]->SetTextColor(2);
	leg1[2]->SetTextColor(4);
	leg1[3]->SetTextColor(3);
	leg->Draw();
	c1->SaveAs("../temp/plot_thlab.pdf");

	TCanvas* c2 = new TCanvas("c2","ee logy");
	gPad->SetLogy();
	h_thlab->Draw();
	h_thlab_o->Draw("same");
	h_thlab_c->Draw("same");
	h_thlab_t->Draw("same");
	leg->Draw();
	c2->SaveAs("../temp/plot_thlab_logy.pdf");
	
	TCanvas* c3 = new TCanvas("c3","ee hist");
	h_thlab->Draw("hist");
	h_thlab_o->Draw("hist same");
	h_thlab_c->Draw("hist same");
	h_thlab_t->Draw("hist same");
	leg->Draw();
	c3->SaveAs("../temp/plot_thlab_hist.pdf");
	
	TCanvas* c4 = new TCanvas("c4","ee logy hist");
	gPad->SetLogy();
	h_thlab->Draw("hist");
	h_thlab_o->Draw("hist same");
	h_thlab_c->Draw("hist same");
	h_thlab_t->Draw("hist same");
	leg->Draw();
	c4->SaveAs("../temp/plot_thlab_logy_hist.pdf");

	gSystem->Exec(Form("pdfunite ../temp/plot*_thlab*.pdf ../plots/thlab_ee.pdf"));
	gSystem->Exec(Form("rm -rf ../temp/plot*_thlab*.pdf"));
}
