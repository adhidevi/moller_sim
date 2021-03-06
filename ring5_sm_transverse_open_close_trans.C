//This script is used to plot the transverse (y vs x) distribution on ring5 and showermax with proper cuts in closed, open and transition sectors.
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void ring5_sm_transverse_open_close_trans(){
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

	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";//directory where the rootfiles exist
	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};//list of root files.
	int nfile = sizeof(file)/sizeof(*file);
	
	TH2F* h_main_tr_c[nfile];//histogram for ring5 transverse distribution
	TH2F* h_main_tr_o[nfile];//histogram for ring5 transverse distribution
	TH2F* h_main_tr_t[nfile];//histogram for ring5 transverse distribution
	TH2F* h_sm_tr_c[nfile];//histogram for showermax transverse distribution
	TH2F* h_sm_tr_o[nfile];//histogram for showermax transverse distribution
	TH2F* h_sm_tr_t[nfile];//histogram for showermax transverse distribution

	int coloro[] = {2,2,2};
	int colorc[] = {4,4,4};
	int colort[] = {3,3,3};

	TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
	TString smCut = "(hit.det==53&&hit.e>1000&&hit.r>500)";//default fut for det 53 (sm)
	
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph).
	TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sector
	TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sector
	TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sector
	TString weight;//define weight inside loop

	TString sim[] = {"ee","ep-el","ep-inel"};//for histogram title

	for(int ifile=0;ifile<nfile;ifile++){//loop over ee, ep-el, and ep-inel
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbinx = 420;
	int nbiny = 420;
	int x_min = -1400;
	int x_max = 1400;
	int y_min = -1400;
	int y_max = 1400;

	h_main_tr_c[ifile] = new TH2F(Form("h_main_tr_c[%d]",ifile),Form("%s transverse on main det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_main_tr_o[ifile] = new TH2F(Form("h_main_tr_o[%d]",ifile),Form("%s transverse on main det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_main_tr_t[ifile] = new TH2F(Form("h_main_tr_t[%d]",ifile),Form("%s transverse on main det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_sm_tr_c[ifile] = new TH2F(Form("h_sm_tr_c[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_sm_tr_o[ifile] = new TH2F(Form("h_sm_tr_o[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
	h_sm_tr_t[ifile] = new TH2F(Form("h_sm_tr_t[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);

	h_main_tr_c[ifile]->SetMarkerColor(colorc[ifile]);
	h_main_tr_o[ifile]->SetMarkerColor(coloro[ifile]);
	h_main_tr_t[ifile]->SetMarkerColor(colort[ifile]);
	h_sm_tr_c[ifile]->SetMarkerColor(colorc[ifile]);
	h_sm_tr_o[ifile]->SetMarkerColor(coloro[ifile]);
	h_sm_tr_t[ifile]->SetMarkerColor(colort[ifile]);

	if(ifile==0)
	weight = "1e-9*rate/2.0";
	else
	weight = "1e-9*rate";
	
	T->Draw(Form("hit.y:-1*hit.x>>h_main_tr_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("hit.y:-1*hit.x>>h_main_tr_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("hit.y:-1*hit.x>>h_main_tr_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	T->Draw(Form("hit.y:-1*hit.x>>h_sm_tr_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("hit.y:-1*hit.x>>h_sm_tr_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("hit.y:-1*hit.x>>h_sm_tr_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),smCut.Data(),phi_tCut.Data()),"goff");
	}
	
	TString label[] = {"open","close","transition"};
	TArc* mainIn = new TArc(0,0,ring5_rmin,0,360);
	TArc* mainOut = new TArc(0,0,ring5_rmax,0,360);
	TArc* smIn = new TArc(0,0,sm_rmin,0,360);
	TArc* smOut = new TArc(0,0,sm_rmax,0,360);
	mainIn->SetLineColor(6);
	mainIn->SetFillColor(0);
	mainIn->SetFillStyle(0);
	mainIn->SetLineWidth(2);
	mainIn->SetLineStyle(7);
	mainOut->SetLineColor(6);
	mainOut->SetFillColor(0);
	mainOut->SetFillStyle(0);
	mainOut->SetLineWidth(2);
	mainOut->SetLineStyle(7);
	smIn->SetLineColor(6);
	smIn->SetFillColor(0);
	smIn->SetFillStyle(0);
	smIn->SetLineWidth(2);
	smIn->SetLineStyle(7);
	smOut->SetLineColor(6);
	smOut->SetFillColor(0);
	smOut->SetFillStyle(0);
	smOut->SetLineWidth(2);
	smOut->SetLineStyle(7);
	TCanvas* cM[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cM[ifile]=new TCanvas(Form("cM[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
	h_main_tr_o[ifile]->Draw();
	h_main_tr_c[ifile]->Draw("same");
	h_main_tr_t[ifile]->Draw("same");
	TLegend* legM = new TLegend(0.70,0.75,0.9,0.9);
	legM->SetBorderSize(0);
	legM->SetFillColor(0);
	legM->SetFillStyle(0);
	TLegendEntry* leg1[3];
	leg1[0] = legM->AddEntry(h_main_tr_o[ifile],"open","p");
	leg1[1] = legM->AddEntry(h_main_tr_c[ifile],"close","p");
	leg1[2] = legM->AddEntry(h_main_tr_t[ifile],"transition","p");
	legM->SetTextSize(0.05);
	leg1[0]->SetTextColor(coloro[ifile]);
	leg1[1]->SetTextColor(colorc[ifile]);
	leg1[2]->SetTextColor(colort[ifile]);
	legM->Draw();
	mainIn->Draw();
	mainOut->Draw();
	cM[ifile]->SaveAs(Form("../temp/ring5_transverse_open_close_transition_%s.pdf",sim[ifile].Data()));
	}
	TCanvas* cMcolz[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cMcolz[ifile]=new TCanvas(Form("cMcolz[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
	h_main_tr_o[ifile]->Draw("colz");
	h_main_tr_c[ifile]->Draw("colz same");
	h_main_tr_t[ifile]->Draw("colz same");
	mainIn->Draw();
	mainOut->Draw();
	cMcolz[ifile]->SaveAs(Form("../temp/ring5_transverse_colz_%s.pdf",sim[ifile].Data()));
	}
	TCanvas* cS[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cS[ifile]=new TCanvas(Form("cS[%d]",ifile),Form("SM %s",sim[ifile].Data()),600,600);
	h_sm_tr_o[ifile]->Draw();
	h_sm_tr_c[ifile]->Draw("same");
	h_sm_tr_t[ifile]->Draw("same");
	TLegend* legS = new TLegend(0.70,0.75,0.9,0.9);
	legS->SetBorderSize(0);
	legS->SetFillColor(0);
	legS->SetFillStyle(0);
	TLegendEntry* leg2[3];
	leg2[0] = legS->AddEntry(h_sm_tr_o[ifile],"open","p");
	leg2[1] = legS->AddEntry(h_sm_tr_c[ifile],"close","p");
	leg2[2] = legS->AddEntry(h_sm_tr_t[ifile],"transition","p");
	legS->SetTextSize(0.05);
	leg2[0]->SetTextColor(coloro[ifile]);
	leg2[1]->SetTextColor(colorc[ifile]);
	leg2[2]->SetTextColor(colort[ifile]);
	legS->Draw();
	smIn->Draw();
	smOut->Draw();
	cS[ifile]->SaveAs(Form("../temp/sm_transverse_open_close_transition_%s.pdf",sim[ifile].Data()));
	}
	TCanvas* cScolz[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	cScolz[ifile]=new TCanvas(Form("cScolz[%d]",ifile),Form("SM %s",sim[ifile].Data()),600,600);
	h_sm_tr_o[ifile]->Draw("colz");
	h_sm_tr_c[ifile]->Draw("colz same");
	h_sm_tr_t[ifile]->Draw("colz same");
	smIn->Draw();
	smOut->Draw();
	cScolz[ifile]->SaveAs(Form("../temp/sm_transverse_colz_%s.pdf",sim[ifile].Data()));
	}
	gSystem->Exec(Form("pdfunite ../temp/ring5_transverse* ../temp/sm_transverse* ../plots/ring5_sm_transverse_open_close_trans_ee_ep-el_ep-inel.pdf"));
	gSystem->Exec(Form("rm -rf ../temp/ring5_transverse* ../temp/sm_transverse*"));
}
