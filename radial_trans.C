//This script is used for testing the ep background in the detector acceptance
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void radial_trans(){
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
	int colorQ[] = {3,6,7};
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/";
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH1F* h_main[nfile];
	TH1F* h_sm[nfile];
	TH2F* h_main_transverse[nfile];
	TH2F* h_sm_transverse[nfile];

	for(int ifile=0;ifile<nfile;ifile++){
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 500;
	int x_min = 500;
	int x_max = 1500;
	double bin_width = (x_max-x_min)/nbin;

	h_main[ifile] = new TH1F(Form("h_main[%d]",ifile),Form("hit distribution on main det (rate*Asym weighted);hit_r (mm);rate*Asym (GHz*ppb)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_sm[ifile] = new TH1F(Form("h_sm[%d]",ifile),Form("hit distribution on sm det (rate*Asym weighted);hit_r (mm);rate*Asym (GHz*ppb)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main_transverse[ifile] = new TH2F(Form("h_main_transverse[%d]",ifile), "transverse on main det;-x (mm);y (mm)", 420, -1400, 1400,420,-1400,1400);
	h_sm_transverse[ifile] = new TH2F(Form("h_sm_transverse[%d]",ifile), "transverse on sm det;- x (mm);y (mm)", 420, -1400, 1400,420,-1400,1400);
	h_main[ifile]->SetLineColor(color[ifile]);
	h_sm[ifile]->SetLineColor(color[ifile]);
	h_main_transverse[ifile]->SetMarkerColor(color[ifile]);
	h_sm_transverse[ifile]->SetMarkerColor(color[ifile]);
	h_main[ifile]->Sumw2();
	h_sm[ifile]->Sumw2();

	if(ifile==0){	
	T->Draw(Form("hit.r>>h_main[%d]",ifile),"-1*A*1.e-9*rate/2.0*(hit.det==28&&hit.e>1000&&hit.r>500)","goff");
	T->Draw(Form("hit.r>>h_sm[%d]",ifile),"-1*A*1.e-9*rate/2.0*(hit.det==53&&hit.e>1000&&hit.r>500)","goff");
	}else{
	T->Draw(Form("hit.r>>h_main[%d]",ifile),"-1*A*1.e-9*rate*(hit.det==28&&hit.e>1000&&hit.r>500)","goff");
	T->Draw(Form("hit.r>>h_sm[%d]",ifile),"-1*A*1.e-9*rate*(hit.det==53&&hit.e>1000&&hit.r>500)","goff");
	}
	T->Draw(Form("hit.y:hit.x>>h_main_transverse[%d]",ifile),"hit.det==28&&hit.e>1000&&hit.r>500","goff");
	T->Draw(Form("hit.y:hit.x>>h_sm_transverse[%d]",ifile),"hit.det==53&&hit.e>1000&&hit.r>500","goff");

	}

	TCanvas* c1 = new TCanvas("c1","Main");
	gPad->SetLogy();
	h_main[0]->Draw();
	h_main[1]->Draw("same");
	h_main[2]->Draw("same");
	TH1F* h_mainQ[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	h_mainQ[ifile] = (TH1F*)h_main[ifile]->Clone(Form("h_mainQ[%d]",ifile));
	h_mainQ[ifile]->GetXaxis()->SetRangeUser(ring5_rmin,ring5_rmax);
	h_mainQ[ifile]->SetLineColor(colorQ[ifile]);
	h_mainQ[ifile]->Draw("same");
	}
	TLegend* legM = new TLegend(0.73,0.65,0.9,0.9);
	legM->SetBorderSize(0);
	legM->SetFillColor(0);
	legM->SetFillStyle(0);
	TLegendEntry* leg1[6];
	leg1[0] = legM->AddEntry(h_main[0],"ee","le");
	leg1[1] = legM->AddEntry(h_main[1],"ep_el","le");
	leg1[2] = legM->AddEntry(h_main[2],"ep_inel","le");
	leg1[3] = legM->AddEntry(h_mainQ[0],"eeQ","le");
	leg1[4] = legM->AddEntry(h_mainQ[1],"ep_elQ","le");
	leg1[5] = legM->AddEntry(h_mainQ[2],"ep_inelQ","le");
	legM->SetTextSize(0.05);
	for(int ifile;ifile<nfile;ifile++){
	leg1[ifile]->SetTextColor(color[ifile]);
	leg1[ifile+3]->SetTextColor(colorQ[ifile]);
	}
	legM->Draw();
	TLine* line[2];
	line[0] = new TLine(ring5_rmin,0.0,ring5_rmin,h_main[0]->GetMaximum());
	line[1] = new TLine(ring5_rmax,0.0,ring5_rmax,h_main[0]->GetMaximum());	
	line[0]->SetLineWidth(2);
	line[1]->SetLineWidth(2);
	line[0]->SetLineColor(kRed-2);
	line[1]->SetLineColor(kRed-2);
	line[0]->Draw();
	line[1]->Draw();
	TArrow* hlineM = new TArrow(ring5_rmin,0.2,ring5_rmax,0.2,0.02,"<|>");
	hlineM->SetLineWidth(2);
	hlineM->SetLineColor(kRed-2);
	hlineM->Draw();
	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.04);
	latex.SetTextColor(kRed-2);
	latex.DrawLatex(0.43,0.2,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
	double ring5_ee_int = h_main[0]->Integral();
	double ring5_ee_accept_int = h_main[0]->Integral(h_main[0]->FindBin(ring5_rmin),h_main[0]->FindBin(ring5_rmax));
	double ring5_ep_el_accept_int = h_main[1]->Integral(h_main[1]->FindBin(ring5_rmin),h_main[1]->FindBin(ring5_rmax));
	double ring5_eff = ring5_ee_accept_int/ring5_ee_int*100.;
	double ring5_bkg = ring5_ep_el_accept_int/ring5_ee_accept_int*100.;
	latex.SetTextColor(1);
	latex.DrawLatex(0.75,0.60,Form("Eff = %.1f %s",ring5_eff,"%"));
	latex.SetTextColor(2);
	latex.DrawLatex(0.75,0.55,Form("Bkg = %.1f %s",ring5_bkg,"%"));
	cout<<"ring5 efficiency: "<<ring5_eff<<"%; ring5 bkg: "<<ring5_bkg<<"%"<<endl;
	c1->SaveAs(Form("../temp/plot1.pdf"));

	TCanvas* c2 = new TCanvas("c2","sm");
	gPad->SetLogy();
	h_sm[0]->Draw();
	h_sm[1]->Draw("same");
	h_sm[2]->Draw("same");
	TH1F* h_smQ[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	h_smQ[ifile] = (TH1F*)h_sm[ifile]->Clone(Form("h_smQ[%d]",ifile));
	h_smQ[ifile]->GetXaxis()->SetRangeUser(sm_rmin,sm_rmax);
	h_smQ[ifile]->SetLineColor(colorQ[ifile]);
	h_smQ[ifile]->Draw("same");
	}
	TLegend* legS = new TLegend(0.73,0.65,0.9,0.9);
	legS->SetBorderSize(0);
	legS->SetFillColor(0);
	legS->SetFillStyle(0);
	TLegendEntry* leg2[nfile];;
	leg2[0] = legS->AddEntry(h_sm[0],"ee","le");
	leg2[1] = legS->AddEntry(h_sm[1],"ep_el","le");
	leg2[2] = legS->AddEntry(h_sm[2],"ep_inel","le");
	leg2[3] = legS->AddEntry(h_smQ[0],"eeQ","le");
	leg2[4] = legS->AddEntry(h_smQ[1],"ep_elQ","le");
	leg2[5] = legS->AddEntry(h_smQ[2],"ep_inelQ","le");
	legS->SetTextSize(0.05);
	for(int ifile=0;ifile<nfile;ifile++){
	leg2[ifile]->SetTextColor(color[ifile]);
	leg2[ifile+3]->SetTextColor(colorQ[ifile]);
	}
	legS->Draw();

	line[0] = new TLine(sm_rmin,0.0,sm_rmin,h_sm[0]->GetMaximum());
	line[1] = new TLine(sm_rmax,0.0,sm_rmax,h_sm[0]->GetMaximum());	
	line[0]->SetLineWidth(2);
	line[1]->SetLineWidth(2);
	line[0]->SetLineColor(kRed-2);
	line[1]->SetLineColor(kRed-2);
	line[0]->Draw();
	line[1]->Draw();
	TArrow* hlineS = new TArrow(sm_rmin,0.2,sm_rmax,0.2,0.02,"<|>");
	hlineS->SetLineWidth(2);
	hlineS->SetLineColor(kRed-2);
	hlineS->Draw();
	latex.SetTextColor(kRed-2);
	latex.DrawLatex(0.52,0.20,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
	double sm_ee_int = h_sm[0]->Integral();
	double sm_ee_accept_int = h_sm[0]->Integral(h_sm[0]->FindBin(sm_rmin),h_sm[0]->FindBin(sm_rmax));
	double sm_ep_el_accept_int = h_sm[1]->Integral(h_sm[1]->FindBin(sm_rmin),h_sm[1]->FindBin(sm_rmax));
	double sm_eff = sm_ee_accept_int/sm_ee_int*100.;
	double sm_bkg = sm_ep_el_accept_int/sm_ee_accept_int*100.;
	latex.SetTextColor(1);
	latex.DrawLatex(0.75,0.60,Form("Eff = %.1f %s",sm_eff,"%"));
	latex.SetTextColor(2);
	latex.DrawLatex(0.75,0.55,Form("Bkg = %.1f %s",sm_bkg,"%"));
	cout<<"sm efficiency: "<<sm_eff<<"%; sm bkg: "<<sm_bkg<<"%"<<endl;
	c2->SaveAs(Form("../temp/plot2.pdf"));

	TCanvas* c3 = new TCanvas("c3","main transverse",600,600);
	h_main_transverse[0]->Draw();
	h_main_transverse[1]->Draw("same");
	h_main_transverse[2]->Draw("same");
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
	c3->SaveAs(Form("../temp/plot3.pdf"));

	TCanvas* c4 = new TCanvas("c4","sm transverse",600,600);
	h_sm_transverse[0]->Draw();
	h_sm_transverse[1]->Draw("same");
	h_sm_transverse[2]->Draw("same");
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
	c4->SaveAs(Form("./temp/plot4.pdf"));
	gSystem->Exec(Form("pdfunite ./temp/* ./plots/ring5_showermax_hit_study.pdf"));
	gSystem->Exec(Form("rm -rf ./temp/*"));
}
