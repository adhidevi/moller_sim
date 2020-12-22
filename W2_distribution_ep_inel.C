//This script is used to plot the W2 distributions for ep-inel.
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void W2_distribution_ep_inel(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	TGaxis::SetMaxDigits(3);
	
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";//directory where the rootfiles exist
	TString file ="main_sm_EP_inelastic/main_sm_EP_inelastic*.root";
	
	TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph)
	TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sectors
	TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sectors
	TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sectors

	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file.Data()));

	int nbin = 500;
	double x_min = 0;
	double x_max = 20;
	double bin_width = (x_max-x_min)/nbin;
	TString weight[] = {"1e-9*rate","1e-9*rate*ev.Q2/GeV^2"};
	TString title[] = {"rate weighted","rate*Q^{2} weighted"};
	TString ytitle[] = {"rate (GHz)","rate*Q^{2} (GHz*(GeV/c)^{2})"};
	int nweight = sizeof(weight)/sizeof(*weight);
	TH1F* h_W2[nweight];//histogram for total th_W2 distribution
	TH1F* h_W2_c[nweight];//histogram for closed sector th_W2 distribution
	TH1F* h_W2_o[nweight];//histogram for open sector th_W2 distribution
	TH1F* h_W2_t[nweight];//histogram for transition sector th_W2 distribution
	for(int iweight=0;iweight<nweight;iweight++){
	h_W2[iweight] = new TH1F(Form("h_W2[%d]",iweight),Form("W^{2} distributions (%s);W^{2} (GeV/c)^{2};%s/%.2f(GeV/c)^{2}",title[iweight].Data(),ytitle[iweight].Data(),bin_width),nbin,x_min,x_max);
	h_W2_o[iweight] = new TH1F(Form("h_W2_o[%d]",iweight),Form("W^{2} distributions (%s);W^{2} (GeV/c)^{2};%s/%.2f(GeV/c)^{2}",title[iweight].Data(),ytitle[iweight].Data(),bin_width),nbin,x_min,x_max);
	h_W2_c[iweight] = new TH1F(Form("h_W2_c[%d]",iweight),Form("W^{2} distributions (%s);W^{2} (GeV/c)^{2};%s/%.2f(GeV/c)^{2}",title[iweight].Data(),ytitle[iweight].Data(),bin_width),nbin,x_min,x_max);
	h_W2_t[iweight] = new TH1F(Form("h_W2_t[%d]",iweight),Form("W^{2} distributions (%s);W^{2} (GeV/c)^{2};%s/%.2f(GeV/c)^{2}",title[iweight].Data(),ytitle[iweight].Data(),bin_width),nbin,x_min,x_max);
	h_W2[iweight]->SetLineColor(1);
	h_W2_o[iweight]->SetLineColor(2);
	h_W2_c[iweight]->SetLineColor(4);
	h_W2_t[iweight]->SetLineColor(3);
	T->Draw(Form("ev.W2/GeV^2>>h_W2[%d]",iweight),Form("%s*(%s)",weight[iweight].Data(),ring5Cut.Data()),"goff");
	T->Draw(Form("ev.W2/GeV^2>>h_W2_o[%d]",iweight),Form("%s*(%s&&%s)",weight[iweight].Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("ev.W2/GeV^2>>h_W2_c[%d]",iweight),Form("%s*(%s&&%s)",weight[iweight].Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("ev.W2/GeV^2>>h_W2_t[%d]",iweight),Form("%s*(%s&&%s)",weight[iweight].Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	}

	TLegend* leg[nweight];
	TLegendEntry* leg0[4];
	TCanvas* c1[nweight];
	TCanvas* c2[nweight];
	TCanvas* c3[nweight];
	TCanvas* c4[nweight];
	for(int iweight=0;iweight<nweight;iweight++){
	leg[iweight] = new TLegend(0.70,0.65,0.9,0.9);
	leg[iweight]->SetBorderSize(0);
	leg[iweight]->SetFillStyle(0);
	leg[iweight]->SetFillColor(0);
	c1[iweight] = new TCanvas(Form("c1[%d]",iweight),"ep-inel");
	h_W2[iweight]->Draw();
	h_W2_o[iweight]->Draw("same");
	h_W2_c[iweight]->Draw("same");
	h_W2_t[iweight]->Draw("same");

	leg0[0] = leg[iweight]->AddEntry(h_W2[iweight],"total ep-inel","le");
	leg0[1] = leg[iweight]->AddEntry(h_W2_o[iweight],"open","le");
	leg0[2] = leg[iweight]->AddEntry(h_W2_c[iweight],"close","le");
	leg0[3] = leg[iweight]->AddEntry(h_W2_t[iweight],"transition","le");
	leg[iweight]->SetTextSize(0.05);
	leg0[0]->SetTextColor(1);
	leg0[1]->SetTextColor(2);
	leg0[2]->SetTextColor(4);
	leg0[3]->SetTextColor(3);
	leg[iweight]->Draw();
	c1[iweight]->SaveAs(Form("../temp/plot%d_W2.pdf",iweight));

	c2[iweight] = new TCanvas(Form("c2[%d]",iweight),"ep-inel logy");
	gPad->SetLogy();
	h_W2[iweight]->Draw();
	h_W2_o[iweight]->Draw("same");
	h_W2_c[iweight]->Draw("same");
	h_W2_t[iweight]->Draw("same");
	leg[iweight]->Draw();
	c2[iweight]->SaveAs(Form("../temp/plot%d_W2_logy.pdf",iweight));

	c3[iweight] = new TCanvas(Form("c3[%d]",iweight),"ep-inel hist");
	h_W2[iweight]->Draw("hist");
	h_W2_o[iweight]->Draw("hist same");
	h_W2_c[iweight]->Draw("hist same");
	h_W2_t[iweight]->Draw("hist same");
	leg[iweight]->Draw();
	c3[iweight]->SaveAs(Form("../temp/plot%d_W2_hist.pdf",iweight));

	c4[iweight] = new TCanvas(Form("c4[%d]",iweight),"ep-inel logy hist");
	gPad->SetLogy();
	h_W2[iweight]->Draw("hist");
	h_W2_o[iweight]->Draw("hist same");
	h_W2_c[iweight]->Draw("hist same");
	h_W2_t[iweight]->Draw("hist same");
	leg[iweight]->Draw();
	c4[iweight]->SaveAs(Form("../temp/plot%d_W2_logy_hist.pdf",iweight));
	}
	gSystem->Exec(Form("pdfunite ../temp/plot*_W2*.pdf ../plots/W2_ep_inel.pdf"));
	gSystem->Exec(Form("rm -rf ../temp/plot*_W2*.pdf"));
}
