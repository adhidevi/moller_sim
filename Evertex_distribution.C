//This script is used to plot the ee, ep-el, and ep-inel E_vertex distribution (rate weighted) in open, close, and transition region.
//It draws histograms directry using T-Draw.
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void Evertex_distribution(){
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetTitleYOffset(1.3);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	TGaxis::SetMaxDigits(3);

	int color[] = {1,1,1};//colors for whole ring histograms
	int coloro[] = {2,2,2};//colors for open sector histograms
	int colorc[] = {4,4,4};//colors for closed sector histograms
	int colort[] = {3,3,3};//colors for transition sector histograms
	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";
//	TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/";
	TString outdir = "/w/halla-scifs17exp/parity/disk1/devi/MOLLER/myScripts";

	TString file[] = {"main_sm_moller/main_sm_moller*.root","main_sm_EP_elastic/main_sm_EP_elastic*.root","main_sm_EP_inelastic/main_sm_EP_inelastic*.root"};
	int nfile = sizeof(file)/sizeof(*file);

	TH1F* h_Evertex[nfile];//histogram for whole ring5
	TH1F* h_Evertex_o[nfile];//histogram for open phi sector ring5
	TH1F* h_Evertex_c[nfile];//histogram for closed phi sector ring5
	TH1F* h_Evertex_t[nfile];//histogram for transition phi sector ring5

        TString ring5Cut = "(hit.det==28&&hit.e>1000&&hit.r>500)";//default cut for det 28 (ring5, main)
	
	TString modphi = "fmod(fabs(hit.ph),2*3.14159/7.)";
	double phiperdet = 3.14159/28.;//this is actually 2*pi/28. but I am looking for fabs(hit.ph).
        TString phi_cCut = Form("(%s<%f||%s>7.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining close sector
        TString phi_oCut = Form("(%s>3.*%f&&%s<5.*%f)",modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining open sector
        TString phi_tCut = Form("((%s>%f&&%s<3.*%f)||(%s>5.*%f&&%s<7.*%f))",modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet,modphi.Data(),phiperdet);//phi cut defining transition sector
	TString weight;//define conditional weight inside the loop below
	TString sim[] = {"ee","ep-el","ep-inel"};

	for(int ifile=0;ifile<nfile;ifile++){//loop over ee, ep-el, and ep-inel
	TChain* T = new TChain("T");
	T->Add(rootfile_dir+Form("%s",file[ifile].Data()));

	int nbin = 110;
	double x_min = 0;
	double x_max = 11;
	double bin_width = (x_max-x_min)/nbin;

	h_Evertex[ifile] = new TH1F(Form("h_Evertex[%d]",ifile),Form("%s E'_{vertex} distribution (rate weighted);E'_{vertex} (GeV/c);rate (GHz)/%.2e(GeV/c)",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_Evertex_o[ifile] = new TH1F(Form("h_Evertex_o[%d]",ifile),Form("%s E'_{vertex} distribution (rate weighted);E'_{vertex} (GeV/c);rate (GHz)/%.2e(GeV/c)",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_Evertex_c[ifile] = new TH1F(Form("h_Evertex_c[%d]",ifile),Form("%s E'_{vertex} distribution (rate weighted);E'_{vertex} (GeV/c);rate (GHz)/%.2e(GeV/c)",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_Evertex_t[ifile] = new TH1F(Form("h_Evertex_t[%d]",ifile),Form("%s E'_{vertex} distribution (rate weighted);E'_{vertex} (GeV/c);rate (GHz)/%.2e(GeV/c)",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
	h_Evertex[ifile]->SetLineColor(color[ifile]);
	h_Evertex_o[ifile]->SetLineColor(coloro[ifile]);
	h_Evertex_c[ifile]->SetLineColor(colorc[ifile]);
	h_Evertex_t[ifile]->SetLineColor(colort[ifile]);

	if(ifile==0)
	weight = "1e-9*rate/2.0";//moller rate needs to be divided by 2.0
	else
	weight = "1e-9*rate";

	T->Draw(Form("hit.e/GeV>>h_Evertex[%d]",ifile),Form("%s*(%s)",weight.Data(),ring5Cut.Data()),"goff");
	T->Draw(Form("hit.e/GeV>>h_Evertex_o[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_oCut.Data()),"goff");
	T->Draw(Form("hit.e/GeV>>h_Evertex_c[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_cCut.Data()),"goff");
	T->Draw(Form("hit.e/GeV>>h_Evertex_t[%d]",ifile),Form("%s*(%s&&%s)",weight.Data(),ring5Cut.Data(),phi_tCut.Data()),"goff");
	cout<<"done with drawing "<<sim[ifile].Data()<<endl;
	}

	TCanvas* c1[nfile];
	TCanvas* c2[nfile];
	TCanvas* c3[nfile];
	TCanvas* c4[nfile];
	for(int ifile=0;ifile<nfile;ifile++){
	c1[ifile] = new TCanvas(Form("c1[%d]",ifile),Form("Main %s",sim[ifile].Data()));
	h_Evertex[ifile]->Draw("hist");
	h_Evertex_o[ifile]->Draw("hist same");
	h_Evertex_c[ifile]->Draw("hist same");
	h_Evertex_t[ifile]->Draw("hist same");

	TLegend* legM0 = new TLegend(0.70,0.65,0.9,0.9);
	legM0->SetBorderSize(0);
	legM0->SetFillColor(0);
	legM0->SetFillStyle(0);
	TLegendEntry* leg0[6];
	leg0[0] = legM0->AddEntry(h_Evertex[ifile],Form("total %s",sim[ifile].Data()),"le");
	leg0[1] = legM0->AddEntry(h_Evertex_o[ifile],"open","le");
	leg0[2] = legM0->AddEntry(h_Evertex_c[ifile],"close","le");
	leg0[3] = legM0->AddEntry(h_Evertex_t[ifile],"transition","le");
	legM0->SetTextSize(0.05);
	leg0[0]->SetTextColor(color[ifile]);
	leg0[1]->SetTextColor(coloro[ifile]);
	leg0[2]->SetTextColor(colorc[ifile]);
	leg0[3]->SetTextColor(colort[ifile]);
	legM0->Draw();
	c1[ifile]->SaveAs(Form("%s/temp/plot%d_Evertex_hist.pdf",outdir.Data(),ifile+1));

	c2[ifile] = new TCanvas(Form("c2[%d]",ifile),Form("Main %s",sim[ifile].Data()));
	gPad->SetLogy();
	h_Evertex[ifile]->Draw("hist");
	h_Evertex_o[ifile]->Draw("hist same");
	h_Evertex_c[ifile]->Draw("hist same");
	h_Evertex_t[ifile]->Draw("hist same");
	legM0->Draw();
	c2[ifile]->SaveAs(Form("%s/temp/plot%d_Evertex_hist_logy.pdf",outdir.Data(),ifile+1));

	c3[ifile] = new TCanvas(Form("c3[%d]",ifile),Form("Main %s",sim[ifile].Data()));
	h_Evertex[ifile]->Draw();
	h_Evertex_o[ifile]->Draw("same");
	h_Evertex_c[ifile]->Draw("same");
	h_Evertex_t[ifile]->Draw("same");
	legM0->Draw();
	c3[ifile]->SaveAs(Form("%s/temp/plot%d_Evertex.pdf",outdir.Data(),ifile+1));

	c4[ifile] = new TCanvas(Form("c4[%d]",ifile),Form("Main %s",sim[ifile].Data()));
	gPad->SetLogy();
	h_Evertex[ifile]->Draw();
	h_Evertex_o[ifile]->Draw("same");
	h_Evertex_c[ifile]->Draw("same");
	h_Evertex_t[ifile]->Draw("same");
	legM0->Draw();
	c4[ifile]->SaveAs(Form("%s/temp/plot%d_Evertex_logy.pdf",outdir.Data(),ifile+1));
	}

	gSystem->Exec(Form("pdfunite %s/temp/plot*_Evertex.pdf %s/plots/Evertex_open_close_trans.pdf",outdir.Data(),outdir.Data()));
//	gSystem->Exec(Form("rm -rf %s/temp/plot*_Evertex.pdf",outdir.Data()));
/*	for(int i=0;i<3;i++){
	cout<<h_Evertex[i]->GetEntries()<<"\t"<<h_Evertex[i]->GetMean()<<"\t"<<h_Evertex[i]->GetRMS()<<endl;
	cout<<h_Evertex_o[i]->GetEntries()<<"\t"<<h_Evertex_o[i]->GetMean()<<"\t"<<h_Evertex_o[i]->GetRMS()<<endl;
	cout<<h_Evertex_c[i]->GetEntries()<<"\t"<<h_Evertex_c[i]->GetMean()<<"\t"<<h_Evertex_c[i]->GetRMS()<<endl;
	cout<<h_Evertex_t[i]->GetEntries()<<"\t"<<h_Evertex_t[i]->GetMean()<<"\t"<<h_Evertex_t[i]->GetRMS()<<endl;
	}
*/
}
