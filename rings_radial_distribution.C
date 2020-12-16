//This script is used for testing the ep background in the various rings of main detector (rate and asymmetry weighted).
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void rings_radial_distribution(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	TGaxis::SetMaxDigits(3);

	Double_t ring_border[] = {640,680,720,770,885,1045,1145};
	TString ring[] = {"ring 1","ring 2","ring 3","ring 4","ring 5","ring 6"};
	int ringN = sizeof(ring)/sizeof(*ring);
	int color[] = {1,2,4};
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/";
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH1F* h_main[nfile];//histogram for main (ring5) detector

	TString ring5_Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//cut for ring5
	TString weight;//define weight within the loop below

	for(int ifile=0;ifile<nfile;ifile++){
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 500;
	int x_min = 500;
	int x_max = 1500;
	double bin_width = (x_max-x_min)/nbin;

	h_main[ifile] = new TH1F(Form("h_main[%d]",ifile),Form("hit distribution on main det (rate*Asym weighted);hit_r (mm);rate*Asym (GHz*ppb)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main[ifile]->SetLineColor(color[ifile]);
	h_main[ifile]->Sumw2();
	
	if(ifile==0)
	weight = "-1*A*1.e-9*rate/2.0";
	else
	weight = "-1*A*1.e-9*rate";
	T->Draw(Form("hit.r>>h_main[%d]",ifile),Form("%s*(%s)",weight.Data(),ring5_Cut.Data()),"goff");
	}

	TLegend* leg = new TLegend(0.73,0.65,0.9,0.9);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	TLegendEntry* leg1[nfile];
	leg1[0] = leg->AddEntry(h_main[0],"ee","le");
	leg1[1] = leg->AddEntry(h_main[1],"ep_el","le");
	leg1[2] = leg->AddEntry(h_main[2],"ep_inel","le");
	
	leg->SetTextSize(0.05);
	for(int ifile;ifile<nfile;ifile++){
	leg1[ifile]->SetTextColor(color[ifile]);
	}
	TLine* line[ringN];
	for(int iring=0;iring<ringN+1;iring++){
	line[iring] = new TLine(ring_border[iring],0.0,ring_border[iring],h_main[1]->GetMaximum());
	line[iring]->SetLineWidth(2);
	line[iring]->SetLineColor(kRed-2);
	}
	TLatex latex;
	latex.SetTextSize(0.04);
	latex.SetTextColor(kRed-2);
	latex.SetTextAngle(90);

	TCanvas* c1 = new TCanvas("c1","Main");
	gPad->SetLogy(0);
	h_main[1]->Draw();
	h_main[0]->Draw("same");
	h_main[2]->Draw("same");
	leg->Draw();
	for(int iring=0;iring<ringN+1;iring++)
	line[iring]->Draw();
	for(int iring=0;iring<ringN;iring++)
	latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.95*h_main[1]->GetMaximum(),Form("%s",ring[iring].Data()));
	c1->SaveAs(Form("../temp/plot0_Main.pdf"));
	
	TCanvas* c2 = new TCanvas("c2","Main logy");
	gPad->SetLogy();
	h_main[1]->Draw();
	h_main[0]->Draw("same");
	h_main[2]->Draw("same");
	leg->Draw();
	for(int iring=0;iring<ringN+1;iring++)
	line[iring]->Draw();
	for(int iring=0;iring<ringN;iring++)
	latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.5*h_main[1]->GetMaximum(),Form("%s",ring[iring].Data()));
	c2->SaveAs(Form("../temp/plot1_Main.pdf"));

	TCanvas* c3 = new TCanvas("c3","Main hist");
	gPad->SetLogy(0);
	h_main[1]->Draw("hist");
	h_main[0]->Draw("hist same");
	h_main[2]->Draw("hist same");
	leg->Draw();
	for(int iring=0;iring<ringN+1;iring++)
	line[iring]->Draw();
	for(int iring=0;iring<ringN;iring++)
	latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.90*h_main[1]->GetMaximum(),Form("%s",ring[iring].Data()));
	c3->SaveAs(Form("../temp/plot2_Main.pdf"));

	TCanvas* c4 = new TCanvas("c4","Main hist logy");
	gPad->SetLogy(1);
	h_main[1]->Draw("hist");
	h_main[0]->Draw("hist same");
	h_main[2]->Draw("hist same");
	leg->Draw();
	for(int iring=0;iring<ringN+1;iring++)
	line[iring]->Draw();
	for(int iring=0;iring<ringN;iring++)
	latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.5*h_main[1]->GetMaximum(),Form("%s",ring[iring].Data()));
	c4->SaveAs(Form("../temp/plot3_Main.pdf"));

	gSystem->Exec(Form("pdfunite ../temp/plot*_Main.pdf ../plots/rings_radial_distribution.C"));
	gSystem->Exec(Form("rm -rf ../temp/plot*_Main.pdf"));

	double ee_int, ep_el_int, ep_inel_int;
	double ee_accept_int[ringN], ep_el_accept_int[ringN],ep_inel_accept_int[ringN], ee_percent[ringN], ep_el_percent[ringN], ep_inel_percent[ringN];
	ee_int = h_main[0]->Integral();
	ep_el_int = h_main[1]->Integral();
	ep_inel_int = h_main[2]->Integral();
	for(int iring=0;iring<ringN;iring++){
	ee_accept_int[iring] = h_main[0]->Integral(h_main[0]->FindBin(ring_border[iring]),h_main[0]->FindBin(ring_border[iring+1]));
	ep_el_accept_int[iring] = h_main[1]->Integral(h_main[1]->FindBin(ring_border[iring]),h_main[1]->FindBin(ring_border[iring+1]));
	ep_inel_accept_int[iring] = h_main[2]->Integral(h_main[2]->FindBin(ring_border[iring]),h_main[2]->FindBin(ring_border[iring+1]));
	ee_percent[iring] = ee_accept_int[iring]/ee_int*100.;
	ep_el_percent[iring] = ep_el_accept_int[iring]/ep_el_int*100.;
	ep_inel_percent[iring] = ep_inel_accept_int[iring]/ep_inel_int*100.;
	cout<<ring[iring].Data()<<"\t ee: "<<ee_percent[iring]<<"\t ep_el: "<<ep_el_percent[iring]<<"\t ep_inel: "<<ep_inel_percent[iring]<<endl;
	}
}
