//plot radial W2 distribution
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void wsq_distribution(){
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
   
   int nbin = 500;
   double x_min = 0;
   double x_max = 20;

   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_W2q2[nfile];//Q-square total
   TH1F* h_W2[nfile];//W-square total
   TH1F* h_W2q2_o[nfile];//Q-square open sectors
   TH1F* h_W2_o[nfile];//W-square open sectors
   TH1F* h_W2q2_c[nfile];//Q-square close sectors
   TH1F* h_W2_c[nfile];//W-square open sectors
   TH1F* h_W2q2_t[nfile];//Q-square transition sectors
   TH1F* h_W2_t[nfile];//W-square transition sectors

   for(int ifile=0;ifile<nfile;ifile++){
      h_W2q2[ifile] = new TH1F(Form("h_W2q2[%d]",ifile),Form("%s W^{2} distribution (rate*Q2 weighted);W^{2} (GeV/c)^{2};rate*Q2 (GHz*(GeV/c)^{2})/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2q2_o[ifile] = new TH1F(Form("h_W2q2o[%d]",ifile),Form("%s W^{2} distribution (rate*Q2 weighted);W^{2} (GeV/c)^{2};rate*Q2 (GHz*(GeV/c)^{2})/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2q2_c[ifile] = new TH1F(Form("h_W2q2c[%d]",ifile),Form("%s W^{2} distribution (rate*Q2 weighted);W^{2} (GeV/c)^{2};rate*Q2 (GHz*(GeV/c)^{2})/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2q2_t[ifile] = new TH1F(Form("h_W2q2t[%d]",ifile),Form("%s W^{2} distribution (rate*Q2 weighted);W^{2} (GeV/c)^{2};rate*Q2 (GHz*(GeV/c)^{2})/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2[ifile] = new TH1F(Form("h_W2[%d]",ifile),Form("%s W^{2} distribution (rate weighted);W^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2_o[ifile] = new TH1F(Form("h_W2_o[%d]",ifile),Form("%s W^{2} distribution (rate weighted);W^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2_c[ifile] = new TH1F(Form("h_W2_c[%d]",ifile),Form("%s W^{2} distribution (rate weighted);W^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2_t[ifile] = new TH1F(Form("h_W2_t[%d]",ifile),Form("%s W^{2} distribution (rate weighted);W^{2} (GeV/c)^{2};rate (GHz)/%.2f(GeV/c)^{2}",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_W2q2[ifile]->SetLineColor(color[ifile]);
      h_W2q2_o[ifile]->SetLineColor(coloro[ifile]);
      h_W2q2_c[ifile]->SetLineColor(colorc[ifile]);
      h_W2q2_t[ifile]->SetLineColor(colort[ifile]);
      h_W2[ifile]->SetLineColor(color[ifile]);
      h_W2_o[ifile]->SetLineColor(coloro[ifile]);
      h_W2_c[ifile]->SetLineColor(colorc[ifile]);
      h_W2_t[ifile]->SetLineColor(colort[ifile]);
      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("ev", &fEv);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), asym(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(-1.e-12), Wsq(-1.e-12), Qsq(-1.e-12);
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
            Wsq = fEv->W2;
            Qsq = fEv->Q2;
//            if(ifile==0)
//            rate = fRate/2./1.e9;
//            else
            rate = fRate/1.e9;
            phi = fHit->at(pk).ph;
            if(phi<0) phi +=2.0*3.14159;
            modphi = fmod(phi,2.0*3.14159/7.);
            if(detector==28 && energy>1000 && hitr>500){
              h_W2q2[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
              h_W2[ifile]->Fill(Wsq/1.e6,rate);
              if(modphi<3.14159/28.){
                h_W2q2_c[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
                h_W2_c[ifile]->Fill(Wsq/1.e6,rate);
              }else if(modphi<3.0*3.14159/28.){
                h_W2q2_t[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
                h_W2_t[ifile]->Fill(Wsq/1.e6,rate);
              }else if(modphi<5.0*3.14159/28.){
                h_W2q2_o[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
                h_W2_o[ifile]->Fill(Wsq/1.e6,rate);
              }else if(modphi<7.0*3.14159/28.){
                h_W2q2_t[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
                h_W2_t[ifile]->Fill(Wsq/1.e6,rate);
              }else{
                h_W2q2_c[ifile]->Fill(Wsq/1.e6,rate*Qsq/1.e6);
                h_W2_c[ifile]->Fill(Wsq/1.e6,rate);
              }
            }
         }
      }
   }
//Now we plot the W-square rate weighted with open, close, and transition sectors cuts
   TCanvas* cw2_1 = new TCanvas("cw2_1","cw2");
   h_W2[2]->Draw();
   h_W2_o[2]->Draw("same");
   h_W2_c[2]->Draw("same");
   h_W2_t[2]->Draw("same");
   TLegend* leg1 = new TLegend(0.70,0.65,0.9,0.9);
   leg1->SetBorderSize(0);
   leg1->SetFillColor(0);
   leg1->SetFillStyle(0);
   leg1->SetTextSize(0.05);
   TLegendEntry* leg11[4];
   leg11[0] = leg1->AddEntry(h_W2[2],"total","le");
   leg11[1] = leg1->AddEntry(h_W2_o[2],"open","le");
   leg11[2] = leg1->AddEntry(h_W2_c[2],"close","le");
   leg11[3] = leg1->AddEntry(h_W2_t[2],"transition","le");
   leg11[0]->SetLineColor(1);
   leg11[1]->SetLineColor(2);
   leg11[2]->SetLineColor(4);
   leg11[3]->SetLineColor(3);
   leg1->Draw();
   cw2_1->SaveAs(Form("../temp/plot1_wsq.pdf"));

   TCanvas* cw2_2 = new TCanvas("cw2_2","cw2");
   gPad->SetLogy();
   h_W2[2]->Draw();
   h_W2_o[2]->Draw("same");
   h_W2_c[2]->Draw("same");
   h_W2_t[2]->Draw("same");
   leg1->Draw();
   cw2_2->SaveAs(Form("../temp/plot2_wsq.pdf"));

   TCanvas* cw2_3 = new TCanvas("cw2_3","cw2");
   h_W2[2]->Draw("hist");
   h_W2_o[2]->Draw("hist same");
   h_W2_c[2]->Draw("hist same");
   h_W2_t[2]->Draw("hist same");
   leg1->Draw();
   cw2_3->SaveAs(Form("../temp/plot3_wsq.pdf"));

   TCanvas* cw2_4 = new TCanvas("cw2_4","cw2");
   gPad->SetLogy();
   h_W2[2]->Draw("hist");
   h_W2_o[2]->Draw("hist same");
   h_W2_c[2]->Draw("hist same");
   h_W2_t[2]->Draw("hist same");
   leg1->Draw();
   cw2_4->SaveAs(Form("../temp/plot4_wsq.pdf"));

//Now we plot the W-square rate*Q2 weighted with open, close, and transition sectors cuts
   TCanvas* cw2q2_1 = new TCanvas("cw2q2_1","cw2q2");
   h_W2q2[2]->Draw();
   h_W2q2_o[2]->Draw("same");
   h_W2q2_c[2]->Draw("same");
   h_W2q2_t[2]->Draw("same");
   TLegend* leg2 = new TLegend(0.70,0.65,0.9,0.9);
   leg2->SetBorderSize(0);
   leg2->SetFillColor(0);
   leg2->SetFillStyle(0);
   leg2->SetTextSize(0.05);
   TLegendEntry* leg21[4];
   leg21[0] = leg2->AddEntry(h_W2q2[2],"total","le");
   leg21[1] = leg2->AddEntry(h_W2q2_o[2],"open","le");
   leg21[2] = leg2->AddEntry(h_W2q2_c[2],"close","le");
   leg21[3] = leg2->AddEntry(h_W2q2_t[2],"transition","le");
   leg21[0]->SetLineColor(1);
   leg21[1]->SetLineColor(2);
   leg21[2]->SetLineColor(4);
   leg21[3]->SetLineColor(3);
   leg2->Draw();
   cw2q2_1->SaveAs(Form("../temp/plot1_wsq_q2.pdf"));

   TCanvas* cw2q2_2 = new TCanvas("cw2q2_2","cw2q2");
   gPad->SetLogy();
   h_W2q2[2]->Draw();
   h_W2q2_o[2]->Draw("same");
   h_W2q2_c[2]->Draw("same");
   h_W2q2_t[2]->Draw("same");
   leg2->Draw();
   cw2q2_2->SaveAs(Form("../temp/plot2_wsq_q2.pdf"));

   TCanvas* cw2q2_3 = new TCanvas("cw2q2_3","cw2q2");
   h_W2q2[2]->Draw("hist");
   h_W2q2_o[2]->Draw("hist same");
   h_W2q2_c[2]->Draw("hist same");
   h_W2q2_t[2]->Draw("hist same");
   leg2->Draw();
   cw2q2_3->SaveAs(Form("../temp/plot3_wsq_q2.pdf"));

   TCanvas* cw2q2_4 = new TCanvas("cw2q2_4","cw2q2");
   gPad->SetLogy();
   h_W2q2[2]->Draw("hist");
   h_W2q2_o[2]->Draw("hist same");
   h_W2q2_c[2]->Draw("hist same");
   h_W2q2_t[2]->Draw("hist same");
   leg2->Draw();
   cw2q2_4->SaveAs(Form("../temp/plot4_wsq_q2.pdf"));

//Now combine all pdf files saved in ../temp/ directory and save a single pdf file in ../plots/ directory
   gSystem->Exec(Form("pdfunite ../temp/plot*_wsq.pdf ../temp/plot*_wsq_*.pdf ../plots/W2_distribution.pdf"));
   gSystem->Exec(Form("rm -rf ../temp/plot*.pdf"));
}
