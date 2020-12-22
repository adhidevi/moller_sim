//plot radial Q2 and Energy distribution
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void qsq_E_distribution(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   gStyle->SetPadGridY(1);
   TGaxis::SetMaxDigits(3);

   int color[] = {1,1,1};//for total events in sector cuts
   int coloro[] = {2,2,2};//for open sectors
   int colorc[] = {4,4,4};//for close sectors
   int colort[] = {3,3,3};//for transition sectors

   TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";
   TString file[] = {"main_sm_moller/main_sm_moller_*.root","main_sm_EP_elastic/main_sm_EP_elastic_*.root","main_sm_EP_inelastic/main_sm_EP_inelastic_*.root"};
        int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"Moller (ee)", "Elastic (ep)", "Inelastic (ep)"};
   
   int nbinQ2 = 500;
   double x_minQ2 = 0;
   double x_maxQ2 = 0.06;
   int nbinE = 110;
   double x_minE = 0;
   double x_maxE = 11;

   double bin_widthQ2 = (x_maxQ2-x_minQ2)/nbinQ2;
   double bin_widthE = (x_maxE-x_minE)/nbinE;

   TH1F* h_Q2[nfile];//Q-square total
   TH1F* h_E[nfile];//W-square total
   TH1F* h_Q2_o[nfile];//Q-square open sectors
   TH1F* h_E_o[nfile];//W-square open sectors
   TH1F* h_Q2_c[nfile];//Q-square close sectors
   TH1F* h_E_c[nfile];//W-square open sectors
   TH1F* h_Q2_t[nfile];//Q-square transition sectors
   TH1F* h_E_t[nfile];//W-square transition sectors

   for(int ifile=0;ifile<nfile;ifile++){
      h_Q2[ifile] = new TH1F(Form("h_Q2[%d]",ifile),Form("%s Q^{2} distribution (rate weighted);Q^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_widthQ2),nbinQ2,x_minQ2,x_maxQ2);
      h_Q2_o[ifile] = new TH1F(Form("h_Q2o[%d]",ifile),Form("%s Q^{2} distribution (rate weighted);Q^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_widthQ2),nbinQ2,x_minQ2,x_maxQ2);
      h_Q2_c[ifile] = new TH1F(Form("h_Q2c[%d]",ifile),Form("%s Q^{2} distribution (rate weighted);Q^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_widthQ2),nbinQ2,x_minQ2,x_maxQ2);
      h_Q2_t[ifile] = new TH1F(Form("h_Q2t[%d]",ifile),Form("%s Q^{2} distribution (rate weighted);Q^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_widthQ2),nbinQ2,x_minQ2,x_maxQ2);
      h_E[ifile] = new TH1F(Form("h_E[%d]",ifile),Form("%s Energy distribution (rate weighted);Energy (GeV);rate (GHz)/%.2e(GeV)",sim[ifile].Data(),bin_widthE),nbinE,x_minE,x_maxE);
      h_E_o[ifile] = new TH1F(Form("h_E_o[%d]",ifile),Form("%s Energy distribution (rate weighted);Energy (GeV);rate (GHz)/%.2e(GeV)",sim[ifile].Data(),bin_widthE),nbinE,x_minE,x_maxE);
      h_E_c[ifile] = new TH1F(Form("h_E_c[%d]",ifile),Form("%s Energy distribution (rate weighted);Energy (GeV);rate (GHz)/%.2e(GeV)",sim[ifile].Data(),bin_widthE),nbinE,x_minE,x_maxE);
      h_E_t[ifile] = new TH1F(Form("h_E_t[%d]",ifile),Form("%s Energy distribution (rate weighted);Energy (GeV);rate (GHz)/%.2e(GeV)",sim[ifile].Data(),bin_widthE),nbinE,x_minE,x_maxE);
      h_Q2[ifile]->SetLineColor(color[ifile]);
      h_Q2_o[ifile]->SetLineColor(coloro[ifile]);
      h_Q2_c[ifile]->SetLineColor(colorc[ifile]);
      h_Q2_t[ifile]->SetLineColor(colort[ifile]);
      h_E[ifile]->SetLineColor(color[ifile]);
      h_E_o[ifile]->SetLineColor(coloro[ifile]);
      h_E_c[ifile]->SetLineColor(colorc[ifile]);
      h_E_t[ifile]->SetLineColor(colort[ifile]);
      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("ev", &fEv);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), asym(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(-1.e-12), Qsq(-1.e-12);
      Int_t pid(0);
   
      for(int ientry=0;ientry<nentry;ientry++){
         if(ientry%(nentry/10)==0)
           cout<<"analyzed "<<ientry<<" events!!"<<endl;
         T->GetEntry(ientry);
         for(size_t pk=0;pk<fHit->size();pk++){
            pid = (Int_t)TMath::Abs(fHit->at(pk).pid);
            detector = fHit->at(pk).det;
            energy = fHit->at(pk).e;
            hitr = fHit->at(pk).r;
            asym = -1*fEv->A;
            Qsq = fEv->Q2;
//            if(ifile==0)
//            rate = fRate/2./1.e9;
//            else
            rate = fRate/1.e9;
            phi = fHit->at(pk).ph;
            if(phi<0) phi +=2.0*3.14159;
            modphi = fmod(phi,2.0*3.14159/7.);
            if(detector==28 && energy>1000 && hitr>500){
              h_Q2[ifile]->Fill(Qsq/1.e6,rate);
              h_E[ifile]->Fill(energy/1.e3,rate);
              if(modphi<3.14159/28.){
                h_Q2_c[ifile]->Fill(Qsq/1.e6,rate);
                h_E_c[ifile]->Fill(energy/1.e3,rate);
              }else if(modphi<3.0*3.14159/28.){
                h_Q2_t[ifile]->Fill(Qsq/1.e6,rate);
                h_E_t[ifile]->Fill(energy/1.e3,rate);
              }else if(modphi<5.0*3.14159/28.){
                h_Q2_o[ifile]->Fill(Qsq/1.e6,rate);
                h_E_o[ifile]->Fill(energy/1.e3,rate);
              }else if(modphi<7.0*3.14159/28.){
                h_Q2_t[ifile]->Fill(Qsq/1.e6,rate);
                h_E_t[ifile]->Fill(energy/1.e3,rate);
              }else{
                h_Q2_c[ifile]->Fill(Qsq/1.e6,rate);
                h_E_c[ifile]->Fill(energy/1.e3,rate);
              }
            }
         }
      }
   }
//Now we plot the W-square rate weighted with open, close, and transition sectors cuts
   TCanvas* cq2_1[nfile];
   TCanvas* cq2_2[nfile];
   TCanvas* cq2_3[nfile];
   TCanvas* cq2_4[nfile];
   TCanvas* ce_1[nfile];
   TCanvas* ce_2[nfile];
   TCanvas* ce_3[nfile];
   TCanvas* ce_4[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
   cq2_1[ifile] = new TCanvas(Form("cq2_1[%d]",ifile),"cq2");
   h_E[ifile]->Draw();
   h_E_o[ifile]->Draw("same");
   h_E_c[ifile]->Draw("same");
   h_E_t[ifile]->Draw("same");
   TLegend* leg1 = new TLegend(0.70,0.65,0.9,0.9);
   leg1->SetBorderSize(0);
   leg1->SetFillColor(0);
   leg1->SetFillStyle(0);
   leg1->SetTextSize(0.05);
   TLegendEntry* leg11[4];
   leg11[0] = leg1->AddEntry(h_E[ifile],"total","le");
   leg11[1] = leg1->AddEntry(h_E_o[ifile],"open","le");
   leg11[2] = leg1->AddEntry(h_E_c[ifile],"close","le");
   leg11[3] = leg1->AddEntry(h_E_t[ifile],"transition","le");
   leg11[0]->SetLineColor(1);
   leg11[1]->SetLineColor(2);
   leg11[2]->SetLineColor(4);
   leg11[3]->SetLineColor(3);
   leg1->Draw();
   cq2_1[ifile]->SaveAs(Form("../temp/plot1_qsq_%d.pdf",ifile));

   cq2_2[ifile] = new TCanvas(Form("cq2_2[%d]",ifile),"cq2");
   gPad->SetLogy();
   h_E[ifile]->Draw();
   h_E_o[ifile]->Draw("same");
   h_E_c[ifile]->Draw("same");
   h_E_t[ifile]->Draw("same");
   leg1->Draw();
   cq2_2[ifile]->SaveAs(Form("../temp/plot2_qsq_%d.pdf",ifile));

   cq2_3[ifile] = new TCanvas(Form("cq2_3[%d]",ifile),"cq2");
   h_E[ifile]->Draw("hist");
   h_E_o[ifile]->Draw("hist same");
   h_E_c[ifile]->Draw("hist same");
   h_E_t[ifile]->Draw("hist same");
   leg1->Draw();
   cq2_3[ifile]->SaveAs(Form("../temp/plot3_qsq_%d.pdf",ifile));

   cq2_4[ifile] = new TCanvas(Form("cq2_4[%d]",ifile),"cq2");
   gPad->SetLogy();
   h_E[ifile]->Draw("hist");
   h_E_o[ifile]->Draw("hist same");
   h_E_c[ifile]->Draw("hist same");
   h_E_t[ifile]->Draw("hist same");
   leg1->Draw();
   cq2_4[ifile]->SaveAs(Form("../temp/plot4_qsq_%d.pdf",ifile));

//Now we plot the W-square rate*Q2 weighted with open, close, and transition sectors cuts
   ce_1[ifile] = new TCanvas(Form("ce_1[%d]",ifile),"ce");
   h_Q2[ifile]->Draw();
   h_Q2_o[ifile]->Draw("same");
   h_Q2_c[ifile]->Draw("same");
   h_Q2_t[ifile]->Draw("same");
   TLegend* leg2 = new TLegend(0.70,0.65,0.9,0.9);
   leg2->SetBorderSize(0);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   leg2->SetTextSize(0.05);
   TLegendEntry* leg21[4];
   leg21[0] = leg2->AddEntry(h_Q2[ifile],"total","le");
   leg21[1] = leg2->AddEntry(h_Q2_o[ifile],"open","le");
   leg21[2] = leg2->AddEntry(h_Q2_c[ifile],"close","le");
   leg21[3] = leg2->AddEntry(h_Q2_t[ifile],"transition","le");
   leg21[0]->SetLineColor(1);
   leg21[1]->SetLineColor(2);
   leg21[2]->SetLineColor(4);
   leg21[3]->SetLineColor(3);
   leg2->Draw();
   ce_1[ifile]->SaveAs(Form("../temp/plot1_energy_%d.pdf",ifile));

   ce_2[ifile] = new TCanvas(Form("ce_2[%d]",ifile),"ce");
   gPad->SetLogy();
   h_Q2[ifile]->Draw();
   h_Q2_o[ifile]->Draw("same");
   h_Q2_c[ifile]->Draw("same");
   h_Q2_t[ifile]->Draw("same");
   leg2->Draw();
   ce_2[ifile]->SaveAs(Form("../temp/plot2_energy_%d.pdf",ifile));

   ce_3[ifile] = new TCanvas(Form("ce_3[%d]",ifile),"ce");
   h_Q2[ifile]->Draw("hist");
   h_Q2_o[ifile]->Draw("hist same");
   h_Q2_c[ifile]->Draw("hist same");
   h_Q2_t[ifile]->Draw("hist same");
   leg2->Draw();
   ce_3[ifile]->SaveAs(Form("../temp/plot3_energy_%d.pdf",ifile));

   ce_4[ifile] = new TCanvas(Form("ce_4[%d]",ifile),"ce");
   gPad->SetLogy();
   h_Q2[ifile]->Draw("hist");
   h_Q2_o[ifile]->Draw("hist same");
   h_Q2_c[ifile]->Draw("hist same");
   h_Q2_t[ifile]->Draw("hist same");
   leg2->Draw();
   ce_4[ifile]->SaveAs(Form("../temp/plot4_energy_%d.pdf",ifile));
   }
//Now combine all pdf files saved in ../temp/ directory and save a single pdf file in ../plots/ directory
   gSystem->Exec(Form("pdfunite ../temp/plot*_qsq_0.pdf ../temp/plot*_qsq_1.pdf ../temp/plot*_qsq_2.pdf ../temp/plot*_energy_0.pdf ../temp/plot*_energy_1.pdf ../temp/plot*_energy_2.pdf ../plots/Q2_E_distribution.pdf"));
   gSystem->Exec(Form("rm -rf ../temp/plot*.pdf"));
}
