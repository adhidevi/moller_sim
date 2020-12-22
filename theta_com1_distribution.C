//This script is used to plot the theta com (center of mass) distributions for ee (moller).
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void theta_com1_distribution(){
   gROOT->Reset();
   gStyle->SetOptStat(0);
   gStyle->SetTitleYOffset(1.3);
   gStyle->SetPadGridX(1);
   TGaxis::SetMaxDigits(3);
   
   TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";//directory where the rootfiles exist
   TString file = "main_sm_moller/main_sm_moller_*.root";
   int nfile = sizeof(file)/sizeof(*file);
   
   TH1F* h_thcom;//histograms for total th_thcom distributions
   TH1F* h_thcom_c;//histograms for closed sector th_thcom distributions
   TH1F* h_thcom_o;//histograms for open sector th_thcom distribution
   TH1F* h_thcom_t;//histograms for transition sector th_thcom distributions
   int nbin = 500;
   double x_min = 0;
   double x_max = 180;
   double bin_width = (x_max-x_min)/nbin;
   h_thcom = new TH1F("h_thcom",Form("#theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",bin_width),nbin,x_min,x_max);
   h_thcom_o = new TH1F("h_thcom_o",Form("#theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",bin_width),nbin,x_min,x_max);
   h_thcom_c = new TH1F("h_thcom_c",Form("#theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",bin_width),nbin,x_min,x_max);
   h_thcom_t = new TH1F("h_thcom_t",Form("#theta_{com} distributions (rate weighted);#theta_{com} (deg);rate (GHz)/%.1fdeg",bin_width),nbin,x_min,x_max);
   h_thcom->SetLineColor(1);
   h_thcom_o->SetLineColor(2);
   h_thcom_c->SetLineColor(4);
   h_thcom_t->SetLineColor(3);
   
   TChain* T = new TChain("T");
   T->Add(rootfile_dir+Form("%s",file.Data()));
   Long64_t nentry = T->GetEntries();
   std::vector<remollGenericDetectorHit_t> *fHit =0;
   remollEvent_t *fEv =0;
   Double_t fRate=0.;
   T->SetBranchAddress("hit", &fHit);
   T->SetBranchAddress("ev", &fEv);
   T->SetBranchAddress("rate", &fRate);

   Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), asym(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(-1.e-12), th_com(-1.e-12);
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
         th_com = 57.295779513*fEv->thcom;//convert to deg
         rate = fRate/1.e9;//Convert to GHz
         phi = fHit->at(pk).ph;
        
         if(phi<0) phi +=2.0*3.14159;
         modphi = fmod(phi,2.0*3.14159/7.);//to account for 7 fold collimator acceptance
         if(detector==28 && energy>1000 && hitr>500){
           h_thcom->Fill(th_com,rate);
           if(modphi<3.14159/28.)
           h_thcom_c->Fill(th_com,rate);
           else if(modphi<3.0*3.14159/28.)
           h_thcom_t->Fill(th_com,rate);
           else if(modphi<5.0*3.14159/28.)
           h_thcom_o->Fill(th_com,rate);
           else if(modphi<7.0*3.14159/28.)
           h_thcom_t->Fill(th_com,rate);
           else
           h_thcom_c->Fill(th_com,rate);
         }
      }
   }
   TCanvas* c1 = new TCanvas("c1","thcom");
   h_thcom->Draw();
   h_thcom_o->Draw("same");
   h_thcom_c->Draw("same");
   h_thcom_t->Draw("same");

   TLegend* leg = new TLegend(0.70,0.65,0.9,0.9);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetFillColor(0);
   TLegendEntry* leg1[4];
   leg1[0] = leg->AddEntry(h_thcom,"total ee","le");
   leg1[1] = leg->AddEntry(h_thcom_o,"open","le");
   leg1[2] = leg->AddEntry(h_thcom_c,"close","le");
   leg1[3] = leg->AddEntry(h_thcom_t,"transition","le");
   leg->SetTextSize(0.05);
   leg1[0]->SetTextColor(1);
   leg1[1]->SetTextColor(2);
   leg1[2]->SetTextColor(4);
   leg1[3]->SetTextColor(3);
   leg->Draw();
   c1->SaveAs("../temp/plot1_thcom.pdf");
   
   TCanvas* c2 = new TCanvas("c2","thcom");
   gPad->SetLogy();
   h_thcom->Draw();
   h_thcom_o->Draw("same");
   h_thcom_c->Draw("same");
   h_thcom_t->Draw("same");
   leg->Draw();
   c2->SaveAs("../temp/plot2_thcom_logy.pdf");
   
   TCanvas* c3 = new TCanvas("c3","thcom");
   h_thcom->Draw("hist");
   h_thcom_o->Draw("hist same");
   h_thcom_c->Draw("hist same");
   h_thcom_t->Draw("hist same");
   leg->Draw();
   c3->SaveAs("../temp/plot3_thcom_hist.pdf");
   
   TCanvas* c4 = new TCanvas("c4","thcom");
   gPad->SetLogy();
   h_thcom->Draw("hist");
   h_thcom_o->Draw("hist same");
   h_thcom_c->Draw("hist same");
   h_thcom_t->Draw("hist same");
   leg->Draw();
   c4->SaveAs("../temp/plot4_thcom_logy_hist.pdf");
   gSystem->Exec(Form("pdfunite ../temp/plot*_thcom*.pdf ../plots/thcom_moller.pdf"));
   gSystem->Exec(Form("rm -rf ../temp/plot*_thcom*.pdf"));
}
