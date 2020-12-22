//This script is used for testing the ep background in the detector acceptance (rate and asymmetry weighted).
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
//offset is used to scan the radial position of the quartz
void ring5_sm_radial_distribution(Double_t offset=0.0){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	TGaxis::SetMaxDigits(3);

	Double_t ring5_rmin = 885+offset;
	Double_t ring5_rmax = 1045+offset;
	Double_t sm_rmin = 995+offset;
	Double_t sm_rmax = 1155+offset;
	int color[] = {1,2,4};
	int colorQ[] = {3,6,7};
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH1F* h_main[nfile];//histogram for main (ring5) detector
	TH1F* h_sm[nfile];//histogram for showermax

	TString ring5_Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//cut for ring5
	TString sm_Cut = "(hit.det==53&&hit.e>1000&&hit.r>500)";//cut for showermax
	TString weight;//define weight within the loop below

	for(int ifile=0;ifile<nfile;ifile++){
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 500;
	int x_min = 500;
	int x_max = 1500;
	double bin_width = (x_max-x_min)/nbin;

	h_main[ifile] = new TH1F(Form("h_main[%d]",ifile),Form("hit distribution on main det (rate*Asym weighted), scan %.1f;hit_r (mm);rate*Asym (GHz*ppb)/%.1fmm",offset,bin_width),nbin,x_min,x_max);
	h_sm[ifile] = new TH1F(Form("h_sm[%d]",ifile),Form("hit distribution on sm det (rate*Asym weighted), scan %.1f;hit_r (mm);rate*Asym (GHz*ppb)/%.1fmm",offset,bin_width),nbin,x_min,x_max);
	h_main[ifile]->SetLineColor(color[ifile]);
	h_sm[ifile]->SetLineColor(color[ifile]);
	h_main[ifile]->Sumw2();
	h_sm[ifile]->Sumw2();
	
	if(ifile==0)
	weight = "-1*A*1.e-9*rate/2.0";
	else
	weight = "-1*A*1.e-9*rate";
	T->Draw(Form("hit.r>>h_main[%d]",ifile),Form("%s*(%s)",weight.Data(),ring5_Cut.Data()),"goff");
	T->Draw(Form("hit.r>>h_sm[%d]",ifile),Form("%s*(%s)",weight.Data(),sm_Cut.Data()),"goff");
	}

	TString title[]={"Main","sm","Main hist","sm hist"};
	int nc = sizeof(title)/sizeof(*title);
	TCanvas* c[4];
	TH1F* h_mainQ[nfile];
	TH1F* h_smQ[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	h_mainQ[ifile] = (TH1F*)h_main[ifile]->Clone(Form("h_mainQ[%d]",ifile));
	h_mainQ[ifile]->GetXaxis()->SetRangeUser(ring5_rmin,ring5_rmax);
	h_mainQ[ifile]->SetLineColor(colorQ[ifile]);
	h_smQ[ifile] = (TH1F*)h_sm[ifile]->Clone(Form("h_smQ[%d]",ifile));
	h_smQ[ifile]->GetXaxis()->SetRangeUser(sm_rmin,sm_rmax);
	h_smQ[ifile]->SetLineColor(colorQ[ifile]);
	}
	for(int ic=0;ic<nc;ic++){
	c[ic] = new TCanvas(Form("c[%d]",ic),Form("%s",title[ic].Data()));
	gPad->SetLogy();
	if(ic==0){
	h_main[0]->Draw();
	h_main[1]->Draw("same");
	h_main[2]->Draw("same");
	h_mainQ[0]->Draw("same");
	h_mainQ[1]->Draw("same");
	h_mainQ[2]->Draw("same");
	}
	if(ic==1){
	h_sm[0]->Draw();
	h_sm[1]->Draw("same");
	h_sm[2]->Draw("same");
	h_smQ[0]->Draw("same");
	h_smQ[1]->Draw("same");
	h_smQ[2]->Draw("same");
	}
	if(ic==2){
	h_main[0]->Draw("hist");
	h_main[1]->Draw("hist same");
	h_main[2]->Draw("hist same");
	h_mainQ[0]->Draw("hist same");
	h_mainQ[1]->Draw("hist same");
	h_mainQ[2]->Draw("hist same");
	}
	if(ic==3){
	h_sm[0]->Draw("hist");
	h_sm[1]->Draw("hist same");
	h_sm[2]->Draw("hist same");
	h_smQ[0]->Draw("hist same");
	h_smQ[1]->Draw("hist same");
	h_smQ[2]->Draw("hist same");
	}
	
	TLegend* leg = new TLegend(0.73,0.65,0.9,0.9);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	TLegendEntry* leg1[6];
	if(title[ic].Contains("Main")){
	leg1[0] = leg->AddEntry(h_main[0],"ee","le");
	leg1[1] = leg->AddEntry(h_main[1],"ep_el","le");
	leg1[2] = leg->AddEntry(h_main[2],"ep_inel","le");
	leg1[3] = leg->AddEntry(h_mainQ[0],"eeQ","le");
	leg1[4] = leg->AddEntry(h_mainQ[1],"ep_elQ","le");
	leg1[5] = leg->AddEntry(h_mainQ[2],"ep_inelQ","le");
	}
	if(title[ic].Contains("sm")){
	leg1[0] = leg->AddEntry(h_sm[0],"ee","le");
	leg1[1] = leg->AddEntry(h_sm[1],"ep_el","le");
	leg1[2] = leg->AddEntry(h_sm[2],"ep_inel","le");
	leg1[3] = leg->AddEntry(h_smQ[0],"eeQ","le");
	leg1[4] = leg->AddEntry(h_smQ[1],"ep_elQ","le");
	leg1[5] = leg->AddEntry(h_smQ[2],"ep_inelQ","le");
	}
	leg->SetTextSize(0.05);
	for(int ifile;ifile<nfile;ifile++){
	leg1[ifile]->SetTextColor(color[ifile]);
	leg1[ifile+3]->SetTextColor(colorQ[ifile]);
	}
	leg->Draw();
	TLine* line[2];
	TArrow* hline;
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.SetTextColor(kRed-2);
	c[ic]->Modified();
	c[ic]->Update();
	if(title[ic].Contains("Main")){
	latex.DrawLatex(0.43,0.2,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
	line[0] = new TLine(ring5_rmin,0.0,ring5_rmin,h_main[0]->GetMaximum());
	line[1] = new TLine(ring5_rmax,0.0,ring5_rmax,h_main[0]->GetMaximum());
	hline = new TArrow(ring5_rmin,0.2,ring5_rmax,0.2,0.02,"<|>");
	}
	if(title[ic].Contains("sm")){
	latex.DrawLatex(0.53,0.2,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
	line[0] = new TLine(sm_rmin,0.0,sm_rmin,h_sm[0]->GetMaximum());
	line[1] = new TLine(sm_rmax,0.0,sm_rmax,h_sm[0]->GetMaximum());
	hline = new TArrow(sm_rmin,0.2,sm_rmax,0.2,0.02,"<|>");
	}
	line[0]->SetLineWidth(2);
	line[1]->SetLineWidth(2);
	line[0]->SetLineColor(kRed-2);
	line[1]->SetLineColor(kRed-2);
	line[0]->Draw();
	line[1]->Draw();
	hline->SetLineWidth(2);
	hline->SetLineColor(kRed-2);
	hline->Draw();
	double ee_int, ee_accept_int, ep_el_accept_int, eff, bkg;
	if(title[ic].Contains("Main")){
	ee_int = h_main[0]->Integral();
	ee_accept_int = h_main[0]->Integral(h_main[0]->FindBin(ring5_rmin),h_main[0]->FindBin(ring5_rmax));
	ep_el_accept_int = h_main[1]->Integral(h_main[1]->FindBin(ring5_rmin),h_main[1]->FindBin(ring5_rmax));
	}
	if(title[ic].Contains("sm")){
	ee_int = h_sm[0]->Integral();
	ee_accept_int = h_sm[0]->Integral(h_sm[0]->FindBin(sm_rmin),h_sm[0]->FindBin(sm_rmax));
	ep_el_accept_int = h_sm[1]->Integral(h_sm[1]->FindBin(sm_rmin),h_sm[1]->FindBin(sm_rmax));
	}
	eff = ee_accept_int/ee_int*100.;
	bkg = ep_el_accept_int/ee_accept_int*100.;
	latex.SetTextColor(1);
	latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff,"%"));
	latex.SetTextColor(2);
	latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg,"%"));
	c[ic]->SaveAs(Form("../temp/plot%d_offset%.1f.pdf",ic,offset));
	}
	gSystem->Exec(Form("pdfunite ../temp/plot*_offset%.1f.pdf ../plots/ring5_sm_radial_offset%.1f.pdf",offset,offset));
	gSystem->Exec(Form("rm -rf ../temp/plot*_offset%.1f.pdf",offset));
}
