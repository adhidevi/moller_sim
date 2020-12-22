//plot radial hit distribution on ring5 and showermax plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void radial_distribution(){
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
   int colorsim[] = {1,2,4};//for total events of ee, ep-el, and ep-inel
   int colorsimQ[] = {3,6,7};//for quartz accepted ee, ep-el, and ep-inel
   int color[] = {1,1,1};//for total events in sector cuts
   int coloro[] = {2,2,2};//for open sectors
   int colorc[] = {4,4,4};//for close sectors
   int colort[] = {3,3,3};//for transition sectors

   TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";
   TString file[] = {"main_sm_moller/main_sm_moller_*.root","main_sm_EP_elastic/main_sm_EP_elastic_*.root","main_sm_EP_inelastic/main_sm_EP_inelastic_*.root"};
        int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"Moller (ee)", "Elastic (ep)", "Inelastic (ep)"};
   
   int nbin = 500;
   double x_min = 500;
   double x_max = 1500;
   double bin_width = (x_max-x_min)/nbin;
   TH1F* h_ring5[nfile];
   TH1F* h_sm[nfile];
   TH1F* h_ring5o[nfile];
   TH1F* h_smo[nfile];
   TH1F* h_ring5c[nfile];
   TH1F* h_smc[nfile];
   TH1F* h_ring5t[nfile];
   TH1F* h_smt[nfile];

   for(int ifile=0;ifile<nfile;ifile++){
      h_ring5[ifile] = new TH1F(Form("h_fing5[%d]",ifile),Form("%s distribution on ring 5 (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_ring5o[ifile] = new TH1F(Form("h_fing5o[%d]",ifile),Form("%s distribution on ring 5 (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_ring5c[ifile] = new TH1F(Form("h_fing5c[%d]",ifile),Form("%s distribution on ring 5 (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_ring5t[ifile] = new TH1F(Form("h_fing5t[%d]",ifile),Form("%s distribution on ring 5 (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_sm[ifile] = new TH1F(Form("h_sm[%d]",ifile),Form("%s distribution on sm det (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_smo[ifile] = new TH1F(Form("h_smo[%d]",ifile),Form("%s distribution on sm det (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_smc[ifile] = new TH1F(Form("h_smc[%d]",ifile),Form("%s distribution on sm det (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_smt[ifile] = new TH1F(Form("h_smt[%d]",ifile),Form("%s distribution on sm det (rate*Asym weighted);hit.r (mm);rate*Asym (GHz*ppb)/%.1fmm",sim[ifile].Data(),bin_width),nbin,x_min,x_max);
      h_ring5[ifile]->SetLineColor(colorsim[ifile]);
      h_ring5o[ifile]->SetLineColor(coloro[ifile]);
      h_ring5c[ifile]->SetLineColor(colorc[ifile]);
      h_ring5t[ifile]->SetLineColor(colort[ifile]);
      h_sm[ifile]->SetLineColor(colorsim[ifile]);
      h_smo[ifile]->SetLineColor(coloro[ifile]);
      h_smc[ifile]->SetLineColor(colorc[ifile]);
      h_smt[ifile]->SetLineColor(colort[ifile]);
      h_ring5[ifile]->Sumw2();
      h_ring5o[ifile]->Sumw2();
      h_ring5c[ifile]->Sumw2();
      h_ring5t[ifile]->Sumw2();
      h_sm[ifile]->Sumw2();
      h_smo[ifile]->Sumw2();
      h_smc[ifile]->Sumw2();
      h_smt[ifile]->Sumw2();
      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("ev", &fEv);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), asym(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(-1.e-12);
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
//            if(ifile==0)
//            rate = fRate/2./1.e9;
//            else
            rate = fRate/1.e9;
            phi = fHit->at(pk).ph;
            if(phi<0) phi +=2.0*3.14159;
            modphi = fmod(phi,2.0*3.14159/7.);
            if(detector==28 && energy>1000 && hitr>500){
              h_ring5[ifile]->Fill(hitr,asym*rate);
              if(modphi<3.14159/28.)
              h_ring5c[ifile]->Fill(hitr,asym*rate);
              else if(modphi<3.0*3.14159/28.)
              h_ring5t[ifile]->Fill(hitr,asym*rate);
              else if(modphi<5.0*3.14159/28.)
              h_ring5o[ifile]->Fill(hitr,asym*rate);
              else if(modphi<7.0*3.14159/28.)
              h_ring5t[ifile]->Fill(hitr,asym*rate);
              else
              h_ring5c[ifile]->Fill(hitr,asym*rate);
            }
            if(detector==53 && energy>1000 && hitr>500){
              h_sm[ifile]->Fill(hitr,asym*rate);
              if(modphi<3.14159/28.)
              h_smc[ifile]->Fill(hitr,asym*rate);
              else if(modphi<3.0*3.14159/28.)
              h_smt[ifile]->Fill(hitr,asym*rate);
              else if(modphi<5.0*3.14159/28.)
              h_smo[ifile]->Fill(hitr,asym*rate);
              else if(modphi<7.0*3.14159/28.)
              h_smt[ifile]->Fill(hitr,asym*rate);
              else
              h_smc[ifile]->Fill(hitr,asym*rate);
            }
         }
      }
   }
   TH1F* h_ring5Q[nfile];
   TH1F* h_smQ[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
      h_ring5Q[ifile] = (TH1F*)h_ring5[ifile]->Clone(Form("h_ring5Q[%d]",ifile));
      h_ring5Q[ifile]->GetXaxis()->SetRangeUser(ring5_rmin,ring5_rmax);
      h_ring5Q[ifile]->SetLineColor(colorsimQ[ifile]);
      h_smQ[ifile] = (TH1F*)h_sm[ifile]->Clone(Form("h_smQ[%d]",ifile));
      h_smQ[ifile]->GetXaxis()->SetRangeUser(sm_rmin,sm_rmax);
      h_smQ[ifile]->SetLineColor(colorsimQ[ifile]);
    }
//Let's plot the ring5 radial distribution with explicit error bars and log scale in y
   TCanvas* ring5_1 = new TCanvas("ring5_1");
   gPad->SetLogy(1);
   h_ring5[1]->Draw();
   h_ring5[0]->Draw("same");
   h_ring5[2]->Draw("same");
   h_ring5Q[0]->Draw("same");
   h_ring5Q[1]->Draw("same");
   h_ring5Q[2]->Draw("same");

   TLegend* legM = new TLegend(0.73,0.65,0.9,0.9);
   legM->SetBorderSize(0);
   legM->SetFillColor(0);
   legM->SetFillStyle(0);
   TLegendEntry* legM1[6];
   legM1[0] = legM->AddEntry(h_ring5[0],"ee","le");
   legM1[1] = legM->AddEntry(h_ring5[1],"ep_el","le");
   legM1[2] = legM->AddEntry(h_ring5[2],"ep_inel","le");
   legM1[3] = legM->AddEntry(h_ring5Q[0],"eeQ","le");
   legM1[4] = legM->AddEntry(h_ring5Q[1],"ep_elQ","le");
   legM1[5] = legM->AddEntry(h_ring5Q[2],"ep_inelQ","le");
   legM->SetTextSize(0.05);
   for(int ifile=0;ifile<nfile;ifile++){
   legM1[ifile]->SetTextColor(colorsim[ifile]);
   legM1[ifile+3]->SetTextColor(colorsimQ[ifile]);
   }
   legM->Draw();
   TLine* line[2];
   TArrow* hline;
   TLatex latex;
   latex.SetNDC(1);
   latex.SetTextSize(0.04);
   latex.SetTextColor(kRed-2);

   latex.DrawLatex(0.43,0.80,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
   line[0] = new TLine(ring5_rmin,0.0,ring5_rmin,h_ring5[1]->GetMaximum());
   line[1] = new TLine(ring5_rmax,0.0,ring5_rmax,h_ring5[1]->GetMaximum());
   hline = new TArrow(ring5_rmin,0.2,ring5_rmax,0.2,0.02,"<|>");
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();
   hline->SetLineWidth(2);
   hline->SetLineColor(kRed-2);
   hline->Draw();
   double ee_int, ee_accept_int, ep_el_accept_int, ep_inel_accept_int, eff, bkg_el, bkg_inel;
   ee_int = h_ring5[0]->Integral();
   ee_accept_int = h_ring5Q[0]->Integral();
   ep_el_accept_int = h_ring5Q[1]->Integral();
   ep_inel_accept_int = h_ring5Q[2]->Integral();
   eff = ee_accept_int/ee_int*100.;
   bkg_el = ep_el_accept_int/ee_accept_int*100.;
   bkg_inel = ep_inel_accept_int/ee_accept_int*100.;
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel,"%"));
   ring5_1->SaveAs(Form("../temp/plot1_ring5.pdf"));

//Let's plot the ring5 radial distribution without explicit error bars and log scale in y
   TCanvas* ring5_2 = new TCanvas("ring5_2");
   gPad->SetLogy(1);
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("hist same");
   h_ring5[2]->Draw("hist same");
   h_ring5Q[0]->Draw("hist same");
   h_ring5Q[1]->Draw("hist same");
   h_ring5Q[2]->Draw("hist same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.43,0.80,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel,"%"));
   ring5_2->SaveAs(Form("../temp/plot2_ring5.pdf"));

//Let's plot the ring5 radial distribution with explicit error bars and linear scale
   TCanvas* ring5_3 = new TCanvas("ring5_3");
   h_ring5[1]->Draw();
   h_ring5[0]->Draw("same");
   h_ring5[2]->Draw("same");
   h_ring5Q[0]->Draw("same");
   h_ring5Q[1]->Draw("same");
   h_ring5Q[2]->Draw("same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.43,0.2,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel,"%"));
   ring5_3->SaveAs(Form("../temp/plot3_ring5.pdf"));

//Let's plot the ring5 radial distribution without explicit error bars and linear scale

   TCanvas* ring5_4 = new TCanvas("ring5_4");
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("hist same");
   h_ring5[2]->Draw("hist same");
   h_ring5Q[0]->Draw("hist same");
   h_ring5Q[1]->Draw("hist same");
   h_ring5Q[2]->Draw("hist same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.43,0.2,Form("#splitline{%.0f mm}{quartz}",ring5_rmax-ring5_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel,"%"));
   ring5_4->SaveAs(Form("../temp/plot4_ring5.pdf"));

//Let's plot radial distribution showing various rings locations
   Double_t ring_border[] = {640,680,720,770,885,1045,1145};
   TString ring[] = {"ring 1","ring 2","ring 3","ring 4","ring 5","ring 6"};
   int ringN = sizeof(ring)/sizeof(*ring);

   TLegend* legR = new TLegend(0.73,0.65,0.9,0.9);
   legR->SetBorderSize(0);
   legR->SetFillColor(0);
   legR->SetFillStyle(0);
   TLegendEntry* legR1[nfile];
   legR1[0] = legR->AddEntry(h_ring5[0],"ee","le");
   legR1[1] = legR->AddEntry(h_ring5[1],"ep_el","le");
   legR1[2] = legR->AddEntry(h_ring5[2],"ep_inel","le");

   legR->SetTextSize(0.05);
   for(int ifile;ifile<nfile;ifile++){
      legR1[ifile]->SetTextColor(colorsim[ifile]);
   }
   TLine* lineR[ringN];
   for(int iring=0;iring<ringN+1;iring++){
      lineR[iring] = new TLine(ring_border[iring],0.0,ring_border[iring],h_ring5[1]->GetMaximum());
      lineR[iring]->SetLineWidth(2);
      lineR[iring]->SetLineColor(kRed-2);
   }
   latex.SetTextColor(kRed-2);
   latex.SetTextAngle(90);
   latex.SetNDC(0);

   TCanvas* cR1 = new TCanvas("cR1","Main");
   gPad->SetLogy(0);
   h_ring5[1]->Draw();
   h_ring5[0]->Draw("same");
   h_ring5[2]->Draw("same");
   legR->Draw();
   for(int iring=0;iring<ringN+1;iring++)
      lineR[iring]->Draw();
   for(int iring=0;iring<ringN;iring++)
      latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.95*h_ring5[1]->GetMaximum(),Form("%s",ring[iring].Data()));
   cR1->SaveAs(Form("../temp/plot0_rings.pdf"));

   TCanvas* cR2 = new TCanvas("cR2","Main logy");
   gPad->SetLogy();
   h_ring5[1]->Draw();
   h_ring5[0]->Draw("same");
   h_ring5[2]->Draw("same");
   legR->Draw();
   for(int iring=0;iring<ringN+1;iring++)
   lineR[iring]->Draw();
   for(int iring=0;iring<ringN;iring++)
   latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.5*h_ring5[1]->GetMaximum(),Form("%s",ring[iring].Data()));
   cR2->SaveAs(Form("../temp/plot1_rings.pdf"));

   TCanvas* cR3 = new TCanvas("cR3","Main hist");
   gPad->SetLogy(0);
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("hist same");
   h_ring5[2]->Draw("hist same");
   legR->Draw();
   for(int iring=0;iring<ringN+1;iring++)
   lineR[iring]->Draw();
   for(int iring=0;iring<ringN;iring++)
   latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.90*h_ring5[1]->GetMaximum(),Form("%s",ring[iring].Data()));
   cR3->SaveAs(Form("../temp/plot2_rings.pdf"));

   TCanvas* cR4 = new TCanvas("cR4","Main hist logy");
   gPad->SetLogy(1);
   h_ring5[1]->Draw("hist");
   h_ring5[0]->Draw("hist same");
   h_ring5[2]->Draw("hist same");
   legR->Draw();
   for(int iring=0;iring<ringN+1;iring++)
   lineR[iring]->Draw();
   for(int iring=0;iring<ringN;iring++)
   latex.DrawLatex(ring_border[iring]+(ring_border[iring+1]-ring_border[iring])/2.0,0.5*h_ring5[1]->GetMaximum(),Form("%s",ring[iring].Data()));
   cR4->SaveAs(Form("../temp/plot3_rings.pdf"));

//Let's plot the showermax radial distribution with explicit error bars and log scale in y
   TCanvas* sm_1 = new TCanvas("sm_1");
   gPad->SetLogy(1);
   h_sm[1]->Draw();
   h_sm[0]->Draw("same");
   h_sm[2]->Draw("same");
   h_smQ[0]->Draw("same");
   h_smQ[1]->Draw("same");
   h_smQ[2]->Draw("same");

   TLegend* legS = new TLegend(0.73,0.65,0.9,0.9);
   legS->SetBorderSize(0);
   legS->SetFillColor(0);
   legS->SetFillStyle(0);
   legS->SetTextSize(0.05);
   TLegendEntry* legS1[6];
   legS1[0] = legS->AddEntry(h_sm[0],"ee","le");
   legS1[1] = legS->AddEntry(h_sm[1],"ep_el","le");
   legS1[2] = legS->AddEntry(h_sm[2],"ep_inel","le");
   legS1[3] = legS->AddEntry(h_smQ[0],"eeQ","le");
   legS1[4] = legS->AddEntry(h_smQ[1],"ep_elQ","le");
   legS1[5] = legS->AddEntry(h_smQ[2],"ep_inelQ","le");
   for(int ifile=0;ifile<nfile;ifile++){
   legS1[ifile]->SetTextColor(colorsim[ifile]);
   legS1[ifile+3]->SetTextColor(colorsimQ[ifile]);
   }
   legS->Draw();

   latex.SetTextColor(kRed-2);
   latex.SetTextAngle(0);
   latex.SetNDC(1);
   latex.DrawLatex(0.52,0.80,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
   line[0] = new TLine(sm_rmin,0.0,sm_rmin,h_sm[1]->GetMaximum());
   line[1] = new TLine(sm_rmax,0.0,sm_rmax,h_sm[1]->GetMaximum());
   hline = new TArrow(sm_rmin,0.2,sm_rmax,0.2,0.02,"<|>");
   line[0]->SetLineWidth(2);
   line[1]->SetLineWidth(2);
   line[0]->SetLineColor(kRed-2);
   line[1]->SetLineColor(kRed-2);
   line[0]->Draw();
   line[1]->Draw();
   hline->SetLineWidth(2);
   hline->SetLineColor(kRed-2);
   hline->Draw();
   double ee_int_sm, ee_accept_int_sm, ep_el_accept_int_sm, ep_inel_accept_int_sm, eff_sm, bkg_el_sm, bkg_inel_sm;
   ee_int_sm = h_sm[0]->Integral();
   ee_accept_int_sm = h_smQ[0]->Integral();
   ep_el_accept_int_sm = h_smQ[1]->Integral();
   ep_inel_accept_int_sm = h_smQ[2]->Integral();
   eff_sm = ee_accept_int_sm/ee_int_sm*100.;
   bkg_el_sm = ep_el_accept_int_sm/ee_accept_int_sm*100.;
   bkg_inel_sm = ep_inel_accept_int_sm/ee_accept_int_sm*100.;
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff_sm,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el_sm,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel_sm,"%"));
   sm_1->SaveAs(Form("../temp/plot1_sm.pdf"));

//Let's plot the showermax radial distribution without explicit error bars and log scale in y
   TCanvas* sm_2 = new TCanvas("sm_2");
   gPad->SetLogy(1);
   h_sm[1]->Draw("hist");
   h_sm[0]->Draw("hist same");
   h_sm[2]->Draw("hist same");
   h_smQ[0]->Draw("hist same");
   h_smQ[1]->Draw("hist same");
   h_smQ[2]->Draw("hist same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.52,0.80,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff_sm,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el_sm,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel_sm,"%"));
   sm_2->SaveAs(Form("../temp/plot2_sm.pdf"));

//Let's plot the showermax radial distribution with explicit error bars and linear scale
   TCanvas* sm_3 = new TCanvas("sm_3");
   h_sm[1]->Draw();
   h_sm[0]->Draw("same");
   h_sm[2]->Draw("same");
   h_smQ[0]->Draw("same");
   h_smQ[1]->Draw("same");
   h_smQ[2]->Draw("same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.52,0.2,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff_sm,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el_sm,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel_sm,"%"));
   sm_3->SaveAs(Form("../temp/plot3_sm.pdf"));

//Let's plot the showermax radial distribution without explicit error bars and linear scale
   TCanvas* sm_4 = new TCanvas("sm_4");
   h_sm[1]->Draw("hist");
   h_sm[0]->Draw("hist same");
   h_sm[2]->Draw("hist same");
   h_smQ[0]->Draw("hist same");
   h_smQ[1]->Draw("hist same");
   h_smQ[2]->Draw("hist same");
   legM->Draw();

   latex.SetTextColor(kRed-2);
   latex.DrawLatex(0.52,0.2,Form("#splitline{%.0f mm}{quartz}",sm_rmax-sm_rmin));
   line[0]->Draw();
   line[1]->Draw();
   hline->Draw();
   latex.SetTextColor(3);
   latex.DrawLatex(0.75,0.60,Form("Eff = %.2f %s",eff_sm,"%"));
   latex.SetTextColor(6);
   latex.DrawLatex(0.75,0.55,Form("Bkg = %.2f %s",bkg_el_sm,"%"));
   latex.SetTextColor(7);
   latex.DrawLatex(0.75,0.50,Form("Bkg = %.2f %s",bkg_inel_sm,"%"));
   sm_4->SaveAs(Form("../temp/plot4_sm.pdf"));

//Now we plot the radial distributions on ring5 and showermax planes with open, close, and transition sectors cuts
   TH1F* h_ring5Cpy[nfile];
   TH1F* h_smCpy[nfile];
   TCanvas* cphM1[nfile];
   TCanvas* cphS1[nfile];
   TCanvas* cphM2[nfile];
   TCanvas* cphS2[nfile];
   TCanvas* cphM3[nfile];
   TCanvas* cphS3[nfile];
   TCanvas* cphM4[nfile];
   TCanvas* cphS4[nfile];
   for(int ii=0;ii<nfile;ii++){
      h_ring5Cpy[ii] = (TH1F*)h_ring5[ii]->Clone(Form("h_ring5Cpy[%d]",ii));
      h_smCpy[ii] = (TH1F*)h_sm[ii]->Clone(Form("h_smCpy[%d]",ii));
      h_ring5Cpy[ii]->SetLineColor(color[ii]);
      h_smCpy[ii]->SetLineColor(color[ii]);

      cphM1[ii] = new TCanvas(Form("cphM1[%d]",ii),"ring5");
      h_ring5Cpy[ii]->Draw();
      h_ring5c[ii]->Draw("same");
      h_ring5o[ii]->Draw("same");
      h_ring5t[ii]->Draw("same");
      TLegend* legMph = new TLegend(0.73,0.65,0.9,0.9);
      legMph->SetBorderSize(0);
      legMph->SetFillColor(0);
      legMph->SetFillStyle(0);
      legMph->SetTextSize(0.05);
      TLegendEntry* legMph1[4];
      legMph1[0] = legMph->AddEntry(h_ring5Cpy[ii],"total","le");
      legMph1[1] = legMph->AddEntry(h_ring5o[ii],"open","le");
      legMph1[2] = legMph->AddEntry(h_ring5c[ii],"close","le");
      legMph1[3] = legMph->AddEntry(h_ring5t[ii],"transition","le");
      legMph1[0]->SetTextColor(color[ii]);
      legMph1[1]->SetTextColor(coloro[ii]);
      legMph1[2]->SetTextColor(colorc[ii]);
      legMph1[3]->SetTextColor(colort[ii]);
      legMph->Draw();
      cphM1[ii]->SaveAs(Form("../temp/plot1_ring5_%d.pdf",ii));

      cphM2[ii] = new TCanvas(Form("cphM2[%d]",ii),"ring5");
      gPad->SetLogy();
      h_ring5Cpy[ii]->Draw();
      h_ring5c[ii]->Draw("same");
      h_ring5o[ii]->Draw("same");
      h_ring5t[ii]->Draw("same");
      legMph->Draw();
      cphM2[ii]->SaveAs(Form("../temp/plot2_ring5_%d.pdf",ii));

      cphM3[ii] = new TCanvas(Form("cphM3[%d]",ii),"ring5");
      h_ring5Cpy[ii]->Draw("hist");
      h_ring5c[ii]->Draw("hist same");
      h_ring5o[ii]->Draw("hist same");
      h_ring5t[ii]->Draw("hist same");
      legMph->Draw();
      cphM3[ii]->SaveAs(Form("../temp/plot3_ring5_%d.pdf",ii));

      cphM4[ii] = new TCanvas(Form("cphM4[%d]",ii),"ring5");
      gPad->SetLogy();
      h_ring5Cpy[ii]->Draw("hist");
      h_ring5c[ii]->Draw("hist same");
      h_ring5o[ii]->Draw("hist same");
      h_ring5t[ii]->Draw("hist same");
      legMph->Draw();
      cphM4[ii]->SaveAs(Form("../temp/plot4_ring5_%d.pdf",ii));

      cphS1[ii] = new TCanvas(Form("cphS1[%d]",ii),"sm");
      h_smCpy[ii]->Draw();
      h_smc[ii]->Draw("same");
      h_smo[ii]->Draw("same");
      h_smt[ii]->Draw("same");
      TLegend* legSph = new TLegend(0.73,0.65,0.9,0.9);
      legSph->SetBorderSize(0);
      legSph->SetFillColor(0);
      legSph->SetFillStyle(0);
      legSph->SetTextSize(0.05);
      TLegendEntry* legSph1[4];
      legSph1[0] = legSph->AddEntry(h_smCpy[ii],"total","le");
      legSph1[1] = legSph->AddEntry(h_smo[ii],"open","le");
      legSph1[2] = legSph->AddEntry(h_smc[ii],"close","le");
      legSph1[3] = legSph->AddEntry(h_smt[ii],"transition","le");
      legSph1[0]->SetTextColor(color[ii]);
      legSph1[1]->SetTextColor(coloro[ii]);
      legSph1[2]->SetTextColor(colorc[ii]);
      legSph1[3]->SetTextColor(colort[ii]);
      legSph->Draw();
      cphS1[ii]->SaveAs(Form("../temp/plot1_sm_%d.pdf",ii));

      cphS2[ii] = new TCanvas(Form("cphS2[%d]",ii),"sm");
      gPad->SetLogy();
      h_smCpy[ii]->Draw();
      h_smc[ii]->Draw("same");
      h_smo[ii]->Draw("same");
      h_smt[ii]->Draw("same");
      legMph->Draw();
      cphS2[ii]->SaveAs(Form("../temp/plot2_sm_%d.pdf",ii));

      cphS3[ii] = new TCanvas(Form("cphS3[%d]",ii),"sm");
      h_smCpy[ii]->Draw("hist");
      h_smc[ii]->Draw("hist same");
      h_smo[ii]->Draw("hist same");
      h_smt[ii]->Draw("hist same");
      legMph->Draw();
      cphS3[ii]->SaveAs(Form("../temp/plot3_sm_%d.pdf",ii));

      cphS4[ii] = new TCanvas(Form("cphS4[%d]",ii),"sm");
      gPad->SetLogy();
      h_smCpy[ii]->Draw("hist");
      h_smc[ii]->Draw("hist same");
      h_smo[ii]->Draw("hist same");
      h_smt[ii]->Draw("hist same");
      legMph->Draw();
      cphS4[ii]->SaveAs(Form("../temp/plot4_sm_%d.pdf",ii));
   }
//Now combine all pdf files saved in ../temp/ directory and save a single pdf file in ../plots/ directory
   gSystem->Exec(Form("pdfunite ../temp/plot*_ring5.pdf ../temp/plot*_sm.pdf ../temp/plot*_rings.pdf ../temp/plot*_ring5_*.pdf ../temp/plot*_sm_*.pdf ../plots/radial_distribution.pdf"));
   gSystem->Exec(Form("rm -rf ../temp/plot*.pdf"));
}
