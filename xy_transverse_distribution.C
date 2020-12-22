//plot y vs x (transverse) looking downstream from the scattering chamber hit distribution on ring5 and showermax plane
//
#include "remolltypes.hh"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void xy_transverse_distribution(){
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
   int coloro[] = {2,2,2};//for open sectors
   int colorc[] = {4,4,4};//for close sectors
   int colort[] = {3,3,3};//for transition sectors

   TString rootfile_dir = "/lustre/expphy/volatile/halla/parity/adhidevi/remoll_rootfiles/fieldmap_v2_50M/";
   TString file[] = {"main_sm_moller/main_sm_moller_*.root","main_sm_EP_elastic/main_sm_EP_elastic_*.root","main_sm_EP_inelastic/main_sm_EP_inelastic_*.root"};
        int nfile = sizeof(file)/sizeof(*file);
   TString sim[] = {"Moller (ee)", "Elastic (ep)", "Inelastic (ep)"};
   
   int nbinx = 420;
   int nbiny = 420;
   double x_min = -1400;
   double x_max = 1400;
   double y_min = -1400;
   double y_max = 1400;
   TH2F* h_ring5_tr[nfile];//histogram for ring5 transverse distribution
   TH2F* h_ring5_tr_c[nfile];//histogram for ring5 transverse distribution
   TH2F* h_ring5_tr_o[nfile];//histogram for ring5 transverse distribution
   TH2F* h_ring5_tr_t[nfile];//histogram for ring5 transverse distribution
   TH2F* h_sm_tr[nfile];//histogram for showermax transverse distribution
   TH2F* h_sm_tr_c[nfile];//histogram for showermax transverse distribution
   TH2F* h_sm_tr_o[nfile];//histogram for showermax transverse distribution
   TH2F* h_sm_tr_t[nfile];//histogram for showermax transverse distribution

   for(int ifile=0;ifile<nfile;ifile++){
      h_ring5_tr[ifile] = new TH2F(Form("h_ring5_tr[%d]",ifile),Form("%s transverse on ring5 det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_ring5_tr_c[ifile] = new TH2F(Form("h_ring5_tr_c[%d]",ifile),Form("%s transverse on ring5 det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_ring5_tr_o[ifile] = new TH2F(Form("h_ring5_tr_o[%d]",ifile),Form("%s transverse on ring5 det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_ring5_tr_t[ifile] = new TH2F(Form("h_ring5_tr_t[%d]",ifile),Form("%s transverse on ring5 det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_sm_tr[ifile] = new TH2F(Form("h_sm_tr[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_sm_tr_c[ifile] = new TH2F(Form("h_sm_tr_c[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_sm_tr_o[ifile] = new TH2F(Form("h_sm_tr_o[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
      h_sm_tr_t[ifile] = new TH2F(Form("h_sm_tr_t[%d]",ifile),Form("%s transverse on sm det (rate weighted); -x (mm);y (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);

      h_ring5_tr[ifile]->SetMarkerColor(colorsim[ifile]);
      h_ring5_tr_c[ifile]->SetMarkerColor(colorc[ifile]);
      h_ring5_tr_o[ifile]->SetMarkerColor(coloro[ifile]);
      h_ring5_tr_t[ifile]->SetMarkerColor(colort[ifile]);
      h_sm_tr[ifile]->SetMarkerColor(colorsim[ifile]);
      h_sm_tr_c[ifile]->SetMarkerColor(colorc[ifile]);
      h_sm_tr_o[ifile]->SetMarkerColor(coloro[ifile]);
      h_sm_tr_t[ifile]->SetMarkerColor(colort[ifile]);

      TChain* T = new TChain("T");
      T->Add(rootfile_dir+Form("%s",file[ifile].Data()));
      Long64_t nentry = T->GetEntries();
      std::vector<remollGenericDetectorHit_t> *fHit =0;
      remollEvent_t *fEv =0;
      Double_t fRate=0.;
      T->SetBranchAddress("hit", &fHit);
      T->SetBranchAddress("rate", &fRate);
   
      Float_t energy(-1.e-12), hitr(-1.e-12), detector(-1.e-12), hitx(-1.e-12), hity(-1.e-12), phi(-1.e-12), modphi(-1.e-12), rate(1.e-12);
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
            hitx = fHit->at(pk).x;
            hity = fHit->at(pk).y;
//            if(ifile==0)
//            rate = fRate/2./1.e9;
//            else
            rate = fRate/1.e9;
            phi = fHit->at(pk).ph;
            if(phi<0) phi +=2.0*3.14159;
            modphi = fmod(phi,2.0*3.14159/7.);
            if(detector==28 && energy>1000 && hitr>500){
              h_ring5_tr[ifile]->Fill(-hitx,hity,rate);
              if(modphi<3.14159/28.)
              h_ring5_tr_c[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<3.0*3.14159/28.)
              h_ring5_tr_t[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<5.0*3.14159/28.)
              h_ring5_tr_o[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<7.0*3.14159/28.)
              h_ring5_tr_t[ifile]->Fill(-hitx,hity,rate);
              else
              h_ring5_tr_c[ifile]->Fill(-hitx,hity,rate);
            }
            if(detector==53 && energy>1000 && hitr>500){
              h_sm_tr[ifile]->Fill(-hitx,hity,rate);
              if(modphi<3.14159/28.)
              h_sm_tr_c[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<3.0*3.14159/28.)
              h_sm_tr_t[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<5.0*3.14159/28.)
              h_sm_tr_o[ifile]->Fill(-hitx,hity,rate);
              else if(modphi<7.0*3.14159/28.)
              h_sm_tr_t[ifile]->Fill(-hitx,hity,rate);
              else
              h_sm_tr_c[ifile]->Fill(-hitx,hity,rate);
            }
         }
      }
   }
    TCanvas* c1 = new TCanvas("c1","ring5 transverse",600,600);
    h_ring5_tr[2]->Draw();
    h_ring5_tr[1]->Draw("same");
    h_ring5_tr[0]->Draw("same");
    TLegend* legMT = new TLegend(0.75,0.75,0.9,0.9);
    legMT->SetBorderSize(0);
    legMT->SetFillColor(0);
    legMT->SetFillStyle(0);
    TLegendEntry* legM1[nfile];;
    legM1[0] = legMT->AddEntry(h_ring5_tr[0],"ee","p");
    legM1[1] = legMT->AddEntry(h_ring5_tr[1],"ep_el","p");
    legM1[2] = legMT->AddEntry(h_ring5_tr[2],"ep_inel","p");
    legMT->SetTextSize(0.05);
    for(int ifile=0;ifile<nfile;ifile++){
    legM1[ifile]->SetTextColor(colorsim[ifile]);
    }
    legMT->Draw();
    TArc* ring5In = new TArc(0,0,ring5_rmin,0,360);
    TArc* ring5Out = new TArc(0,0,ring5_rmax,0,360);
    ring5In->SetLineColor(6);
    ring5Out->SetLineColor(6);
    ring5In->SetLineWidth(2);
    ring5In->SetFillColor(0);
    ring5In->SetFillStyle(0);
    ring5In->SetLineStyle(7);
    ring5Out->SetLineWidth(2);
    ring5Out->SetFillColor(0);
    ring5Out->SetFillStyle(0);
    ring5Out->SetLineStyle(7);
    ring5In->Draw();
    ring5Out->Draw();
    c1->SaveAs(Form("../temp/plot1_transverse.pdf"));

    TCanvas* c2 = new TCanvas("c2","sm transverse",600,600);
    h_sm_tr[2]->Draw();
    h_sm_tr[1]->Draw("same");
    h_sm_tr[0]->Draw("same");
    TLegend* legST = new TLegend(0.75,0.75,0.9,0.9);
    legST->SetBorderSize(0);
    legST->SetFillColor(0);
    legST->SetFillStyle(0);
    TLegendEntry* legS1[nfile];;
    legS1[0] = legST->AddEntry(h_sm_tr[0],"ee","p");
    legS1[1] = legST->AddEntry(h_sm_tr[1],"ep_el","p");
    legS1[2] = legST->AddEntry(h_sm_tr[2],"ep_inel","p");
    legST->SetTextSize(0.05);
    for(int ifile=0;ifile<nfile;ifile++){
    legS1[ifile]->SetTextColor(colorsim[ifile]);
    }
    legST->Draw();
    TArc* smIn = new TArc(0,0,sm_rmin,0,360);
    TArc* smOut = new TArc(0,0,sm_rmax,0,360);
    smIn->SetLineColor(6);
    smOut->SetLineColor(6);
    smIn->SetLineWidth(2);
    smIn->SetFillColor(0);
    smIn->SetFillStyle(0);
    smIn->SetLineStyle(7);
    smOut->SetLineWidth(2);
    smOut->SetFillColor(0);
    smOut->SetFillStyle(0);
    smOut->SetLineStyle(7);
    smIn->Draw();
    smOut->Draw();
    c2->SaveAs(Form("../temp/plot2_transverse.pdf"));

    TCanvas* cM[nfile];
    for(int ifile=0;ifile<nfile;ifile++){
    cM[ifile]=new TCanvas(Form("cM[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
    h_ring5_tr_o[ifile]->Draw();
    h_ring5_tr_c[ifile]->Draw("same");
    h_ring5_tr_t[ifile]->Draw("same");
    TLegend* legM = new TLegend(0.70,0.75,0.9,0.9);
    legM->SetBorderSize(0);
    legM->SetFillColor(0);
    legM->SetFillStyle(0);
    TLegendEntry* leg1[3];
    leg1[0] = legM->AddEntry(h_ring5_tr_o[ifile],"open","p");
    leg1[1] = legM->AddEntry(h_ring5_tr_c[ifile],"close","p");
    leg1[2] = legM->AddEntry(h_ring5_tr_t[ifile],"transition","p");
    legM->SetTextSize(0.05);
    leg1[0]->SetTextColor(coloro[ifile]);
    leg1[1]->SetTextColor(colorc[ifile]);
    leg1[2]->SetTextColor(colort[ifile]);
    legM->Draw();
    ring5In->Draw();
    ring5Out->Draw();
    cM[ifile]->SaveAs(Form("../temp/ring5_transverse_%s.pdf",sim[ifile].Data()));
    }

    TCanvas* cMcolz[nfile];
    for(int ifile=0;ifile<nfile;ifile++){
    cMcolz[ifile]=new TCanvas(Form("cMcolz[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
    h_ring5_tr_o[ifile]->Draw("colz");
    h_ring5_tr_c[ifile]->Draw("colz same");
    h_ring5_tr_t[ifile]->Draw("colz same");
    ring5In->Draw();
    ring5Out->Draw();
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
    cS[ifile]->SaveAs(Form("../temp/sm_transverse_%s.pdf",sim[ifile].Data()));
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
    gSystem->Exec(Form("pdfunite ../temp/ring5_transverse* ../temp/sm_transverse* ../plots/transverse_distributions.pdf"));
    gSystem->Exec(Form("rm -rf ../temp/ring5_transverse* ../temp/sm_transverse* ../temp/plot*"));
}
