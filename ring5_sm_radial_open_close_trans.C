//This script is used for testing the ee radial distribution on det plane(rate weighted), open, close, and transition region.
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ring5_sm_radial_open_close_trans(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	TGaxis::SetMaxDigits(3);

	int color[] = {1,2,3};//colors for whole ring histograms
	int coloro[] = {41,51,97};//colors for open sector histograms
	int colorc[] = {4,6,7};//colors for closed sector histograms
	int colort[] = {8,9,30};//colors for transition sector histograms
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/";
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH1F* h_main[nfile];//histogram for whole ring5
	TH1F* h_main_o[nfile];//histogram for open phi sector ring5
	TH1F* h_main_c[nfile];//histogram for closed phi sector ring5
	TH1F* h_main_t[nfile];//histogram for transition phi sector ring5
	TH1F* h_sm[nfile];//histogram whole showermax
	TH1F* h_sm_o[nfile];//histogram for open phi sector showermax
	TH1F* h_sm_c[nfile];//histogram for closed phi sector showermax
	TH1F* h_sm_t[nfile];//histogram for transition phi sector showermax

        TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
        TString smCut = "(hit.det==53&&hit.e>1000&&hit.r>500)";//default cut for det 53 (sm)
	
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph).
        TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sector
        TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sector
        TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sector
	TString weight;//define conditional weight inside the loop below

	for(int ifile=0;ifile<nfile;ifile++){//loop over ee, ep-el, and ep-inel
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 500;
	int x_min = 500;
	int x_max = 1500;
	double bin_width = (x_max-x_min)/nbin;

	h_main[ifile] = new TH1F(Form("h_main[%d]",ifile),Form("hit distribution on main det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main_o[ifile] = new TH1F(Form("h_main_o[%d]",ifile),Form("hit distribution on main det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main_c[ifile] = new TH1F(Form("h_main_c[%d]",ifile),Form("hit distribution on main det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main_t[ifile] = new TH1F(Form("h_main_t[%d]",ifile),Form("hit distribution on main det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_sm[ifile] = new TH1F(Form("h_sm[%d]",ifile),Form("hit distribution on sm det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_sm_o[ifile] = new TH1F(Form("h_sm_o[%d]",ifile),Form("hit distribution on sm det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_sm_c[ifile] = new TH1F(Form("h_sm_c[%d]",ifile),Form("hit distribution on sm det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_sm_t[ifile] = new TH1F(Form("h_sm_t[%d]",ifile),Form("hit distribution on sm det (rate weighted);hit_r (mm);rate (GHz)/%.1fmm",bin_width),nbin,x_min,x_max);
	h_main[ifile]->SetLineColor(color[ifile]);
	h_main_o[ifile]->SetLineColor(coloro[ifile]);
	h_main_c[ifile]->SetLineColor(colorc[ifile]);
	h_main_t[ifile]->SetLineColor(colort[ifile]);
	h_sm[ifile]->SetLineColor(color[ifile]);
	h_sm_o[ifile]->SetLineColor(coloro[ifile]);
	h_sm_c[ifile]->SetLineColor(colorc[ifile]);
	h_sm_t[ifile]->SetLineColor(colort[ifile]);

	if(ifile==0)
	weight = "1e-9*rate/2.0";//moller rate needs to be divided by 2.0
	else
	weight = "1e-9*rate";

	T->Draw(Form("hit.r>>h_main[%d]",ifile),Form("%s*(%s)",weight.Data(),ring5Cut.Data()),"goff");
	T->Draw(Form("hit.r>>h_main_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_main_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_main_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_sm[%d]",ifile),Form("%s*(%s)",weight.Data(),smCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_sm_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_sm_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("hit.r>>h_sm_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_tCut.Data()),"goff");
	}

	TString sim[] = {"ee","ep-el","ep-inel"};

	TCanvas* cM[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cM[ifile] = new TCanvas(Form("cM[%d]",ifile),Form("Main %s",sim[ifile].Data()));
	h_main[ifile]->Draw("hist");
	h_main_o[ifile]->Draw("hist same");
	h_main_c[ifile]->Draw("hist same");
	h_main_t[ifile]->Draw("hist same");

	TLegend* legM0 = new TLegend(0.70,0.65,0.9,0.9);
	legM0->SetBorderSize(0);
	legM0->SetFillColor(0);
	legM0->SetFillStyle(0);
	TLegendEntry* leg0[6];
	leg0[0] = legM0->AddEntry(h_main[ifile],Form("total %s",sim[ifile].Data()),"le");
	leg0[1] = legM0->AddEntry(h_main_o[ifile],"open","le");
	leg0[2] = legM0->AddEntry(h_main_c[ifile],"close","le");
	leg0[3] = legM0->AddEntry(h_main_t[ifile],"transition","le");
	legM0->SetTextSize(0.05);
	leg0[0]->SetTextColor(color[ifile]);
	leg0[1]->SetTextColor(coloro[ifile]);
	leg0[2]->SetTextColor(colorc[ifile]);
	leg0[3]->SetTextColor(colort[ifile]);
	legM0->Draw();

	cM[ifile]->SaveAs(Form("../temp/plot%d.pdf",ifile+1));
	}

	TCanvas* cS[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cS[ifile] = new TCanvas(Form("cS[%d]",ifile),Form("SM %s",sim[ifile].Data()));
	h_sm[ifile]->Draw("hist");
	h_sm_o[ifile]->Draw("hist same");
	h_sm_c[ifile]->Draw("hist same");
	h_sm_t[ifile]->Draw("hist same");

	TLegend* legM3 = new TLegend(0.70,0.65,0.9,0.9);
	legM3->SetBorderSize(0);
	legM3->SetFillColor(0);
	legM3->SetFillStyle(0);
	TLegendEntry* leg3[6];
	leg3[0] = legM3->AddEntry(h_sm[ifile],Form("total %s",sim[ifile].Data()),"le");
	leg3[1] = legM3->AddEntry(h_sm_o[ifile],"open","le");
	leg3[2] = legM3->AddEntry(h_sm_c[ifile],"close","le");
	leg3[3] = legM3->AddEntry(h_sm_t[ifile],"transition","le");
	legM3->SetTextSize(0.05);
	leg3[0]->SetTextColor(color[ifile]);
	leg3[1]->SetTextColor(coloro[ifile]);
	leg3[2]->SetTextColor(colorc[ifile]);
	leg3[3]->SetTextColor(colort[ifile]);
	legM3->Draw();

	cS[ifile]->SaveAs(Form("../temp/plot%d.pdf",ifile+4));
	}

	for(int i=0;i<3;i++){
	cout<<h_main[i]->GetEntries()<<"\t"<<h_main[i]->GetMean()<<"\t"<<h_main[i]->GetRMS()<<endl;
	cout<<h_main_o[i]->GetEntries()<<"\t"<<h_main_o[i]->GetMean()<<"\t"<<h_main_o[i]->GetRMS()<<endl;
	cout<<h_main_c[i]->GetEntries()<<"\t"<<h_main_c[i]->GetMean()<<"\t"<<h_main_c[i]->GetRMS()<<endl;
	cout<<h_main_t[i]->GetEntries()<<"\t"<<h_main_t[i]->GetMean()<<"\t"<<h_main_t[i]->GetRMS()<<endl;
	cout<<h_sm[i]->GetEntries()<<"\t"<<h_sm[i]->GetMean()<<"\t"<<h_sm[i]->GetRMS()<<endl;
	cout<<h_sm_o[i]->GetEntries()<<"\t"<<h_sm_o[i]->GetMean()<<"\t"<<h_sm_o[i]->GetRMS()<<endl;
	cout<<h_sm_c[i]->GetEntries()<<"\t"<<h_sm_c[i]->GetMean()<<"\t"<<h_sm_c[i]->GetRMS()<<endl;
	cout<<h_sm_t[i]->GetEntries()<<"\t"<<h_sm_t[i]->GetMean()<<"\t"<<h_sm_t[i]->GetRMS()<<endl;
	}
}
