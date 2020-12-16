//This script is used to plot transverse distribution in the detector acceptance
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ring5_sm_transverse_distribution(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	TGaxis::SetMaxDigits(3);

	Double_t ring5_rmin = 885;
	Double_t ring5_rmax = 1045;
	Double_t sm_rmin = 995;
	Double_t sm_rmax = 1155;
	int color[] = {1,2,4};
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/";
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH2F* h_main_transverse[nfile];
	TH2F* h_sm_transverse[nfile];

	for(int ifile=0;ifile<nfile;ifile++){
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbinx = 420;
	int nbiny = 420;
	int x_min = -1400;
	int x_max = 1400;
	int y_min = -1400;
	int y_max = 1400;

	h_main_transverse[ifile] = new TH2F(Form("h_main_transverse[%d]",ifile), "transverse on main det;-x (mm);y (mm)",nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_sm_transverse[ifile] = new TH2F(Form("h_sm_transverse[%d]",ifile), "transverse on sm det;- x (mm);y (mm)",nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_main_transverse[ifile]->SetMarkerColor(color[ifile]);
	h_sm_transverse[ifile]->SetMarkerColor(color[ifile]);

	T->Draw(Form("hit.y:hit.x>>h_main_transverse[%d]",ifile),"hit.det==28&&hit.e>1000&&hit.r>500","goff");
	T->Draw(Form("hit.y:hit.x>>h_sm_transverse[%d]",ifile),"hit.det==53&&hit.e>1000&&hit.r>500","goff");
	}

	TCanvas* c1 = new TCanvas("c1","main transverse",600,600);
	h_main_transverse[2]->Draw();
	h_main_transverse[1]->Draw("same");
	h_main_transverse[0]->Draw("same");
	TLegend* legMT = new TLegend(0.75,0.75,0.9,0.9);
	legMT->SetBorderSize(0);
	legMT->SetFillColor(0);
	legMT->SetFillStyle(0);
	TLegendEntry* leg3[nfile];;
	leg3[0] = legMT->AddEntry(h_main_transverse[0],"ee","le");
	leg3[1] = legMT->AddEntry(h_main_transverse[1],"ep_el","le");
	leg3[2] = legMT->AddEntry(h_main_transverse[2],"ep_inel","le");
	legMT->SetTextSize(0.05);
	for(int ifile=0;ifile<nfile;ifile++){
	leg3[ifile]->SetTextColor(color[ifile]);
	}
	legMT->Draw();
	TArc* mainIn = new TArc(0,0,ring5_rmin,0,360);
	TArc* mainOut = new TArc(0,0,ring5_rmax,0,360);
	mainIn->SetLineColor(6);
	mainOut->SetLineColor(6);
	mainIn->SetLineWidth(2);
	mainIn->SetFillColor(0);
	mainIn->SetFillStyle(0);
	mainIn->SetLineStyle(7);
	mainOut->SetLineWidth(2);
	mainOut->SetFillColor(0);
	mainOut->SetFillStyle(0);
	mainOut->SetLineStyle(7);
	mainIn->Draw();
	mainOut->Draw();
	c1->SaveAs(Form("../temp/plot1.pdf"));

	TCanvas* c2 = new TCanvas("c2","sm transverse",600,600);
	h_sm_transverse[2]->Draw();
	h_sm_transverse[1]->Draw("same");
	h_sm_transverse[0]->Draw("same");
	TLegend* legST = new TLegend(0.75,0.75,0.9,0.9);
	legST->SetBorderSize(0);
	legST->SetFillColor(0);
	legST->SetFillStyle(0);
	TLegendEntry* leg4[nfile];;
	leg4[0] = legST->AddEntry(h_sm_transverse[0],"ee","le");
	leg4[1] = legST->AddEntry(h_sm_transverse[1],"ep_el","le");
	leg4[2] = legST->AddEntry(h_sm_transverse[2],"ep_inel","le");
	legST->SetTextSize(0.05);
	for(int ifile=0;ifile<nfile;ifile++){
	leg4[ifile]->SetTextColor(color[ifile]);
	}
	legST->Draw();
	TArc* smIn = new TArc(0,0,sm_rmin,0,360);
	TArc* smOut = new TArc(0,0,sm_rmax,0,360);
	smIn->SetLineColor(6);
	smOut->SetLineColor(6);
	smIn->SetLineWidth(2);
	smIn->SetFillColor(0);
	smIn->SetFillStyle(0);
	smIn->SetLineStyle(7);
	smOut->SetLineWidth(2);
	smOut->SetFillColor(0);
	smOut->SetFillStyle(0);
	smOut->SetLineStyle(7);
	smIn->Draw();
	smOut->Draw();
	c2->SaveAs(Form("../temp/plot2.pdf"));
	gSystem->Exec(Form("pdfunite ../temp/plot* ../plots/ring5_sm_transverse_distribution.pdf"));
	gSystem->Exec(Form("rm -rf ../temp/plot*"));
}
