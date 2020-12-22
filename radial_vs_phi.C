//This script is used to plot the radial vs phi distribution on ring5 and showermax with proper cuts for open, close and transition sectors.
//
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

void radial_vs_phi(){
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
   
   TH2F* h_ring5[nfile];//histogram for ring5 radial vs phi distribution
   TH2F* h_ring5_c[nfile];//histogram for ring5 close sector
   TH2F* h_ring5_o[nfile];//histogram for ring5 open sector
   TH2F* h_ring5_t[nfile];//histogram for ring5 transition sector
   TH2F* h_sm[nfile];//histogram for showermax radial vs phi distribution
   TH2F* h_sm_c[nfile];//histogram for showermax close sector
   TH2F* h_sm_o[nfile];//histogram for showermax open sector
   TH2F* h_sm_t[nfile];//histogram for showermax transition sector

   int color[] = {1,1,1};
   int coloro[] = {2,2,2};
   int colorc[] = {4,4,4};
   int colort[] = {3,3,3};

   TString sim[] = {"ee","ep-el","ep-inel"};//for histogram title

   int nbinx = 500;
   int nbiny = 500;
   int x_min = 0;
   int x_max = 360;
   int y_min = 500;
   int y_max = 1500;

   for(int ifile=0;ifile<nfile;ifile++){//loop over ee, ep-el, and ep-inel
   h_ring5[ifile] = new TH2F(Form("h_ring5[%d]",ifile),Form("%s radial vs phi on main det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_ring5_c[ifile] = new TH2F(Form("h_ring5_c[%d]",ifile),Form("%s radial vs phi on main det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_ring5_o[ifile] = new TH2F(Form("h_ring5_o[%d]",ifile),Form("%s radial vs phi on main det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_ring5_t[ifile] = new TH2F(Form("h_ring5_t[%d]",ifile),Form("%s radial vs phi on main det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_sm[ifile] = new TH2F(Form("h_sm[%d]",ifile),Form("%s radial vs phi on sm det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_sm_c[ifile] = new TH2F(Form("h_sm_c[%d]",ifile),Form("%s radial vs phi on sm det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_sm_o[ifile] = new TH2F(Form("h_sm_o[%d]",ifile),Form("%s radial vs phi on sm det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);
   h_sm_t[ifile] = new TH2F(Form("h_sm_t[%d]",ifile),Form("%s radial vs phi on sm det (rate weighted); #phi (deg);r (mm)",sim[ifile].Data()),nbinx,x_min,x_max,nbiny,y_min,y_max);

   h_ring5[ifile]->SetMarkerColor(color[ifile]);
   h_ring5_c[ifile]->SetMarkerColor(colorc[ifile]);
   h_ring5_o[ifile]->SetMarkerColor(coloro[ifile]);
   h_ring5_t[ifile]->SetMarkerColor(colort[ifile]);
   h_sm[ifile]->SetMarkerColor(color[ifile]);
   h_sm_c[ifile]->SetMarkerColor(colorc[ifile]);
   h_sm_o[ifile]->SetMarkerColor(coloro[ifile]);
   h_sm_t[ifile]->SetMarkerColor(colort[ifile]);

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
//         if(ifile==0)
//         rate = fRate/2./1.e9;
//         else
         rate = fRate/1.e9;
         phi = fHit->at(pk).ph;
         if(phi<0) phi +=2.0*3.14159;
         modphi = fmod(phi,2.0*3.14159/7.);
         phi = 57.295779513*phi;
         if(detector==28 && energy>1000 && hitr>500){
           h_ring5[ifile]->Fill(phi,hitr);
           if(modphi<3.14159/28. && hitr>ring5_rmin && hitr<ring5_rmax)
           h_ring5_c[ifile]->Fill(phi,hitr);
           else if(modphi<3.0*3.14159/28. && hitr>ring5_rmin && hitr<ring5_rmax)
           h_ring5_t[ifile]->Fill(phi,hitr);
           else if(modphi<5.0*3.14159/28. && hitr>ring5_rmin && hitr<ring5_rmax)
           h_ring5_o[ifile]->Fill(phi,hitr);
           else if(modphi<7.0*3.14159/28. && hitr>ring5_rmin && hitr<ring5_rmax)
           h_ring5_t[ifile]->Fill(phi,hitr);
           else if(hitr>ring5_rmin && hitr<ring5_rmax)
           h_ring5_c[ifile]->Fill(phi,hitr);
         }
         if(detector==53 && energy>1000 && hitr>500){
           h_sm[ifile]->Fill(phi,hitr);
           if(modphi<3.14159/28. && hitr>sm_rmin && hitr<sm_rmax)
           h_sm_c[ifile]->Fill(phi,hitr);
           else if(modphi<3.0*3.14159/28. && hitr>sm_rmin && hitr<sm_rmax)
           h_sm_t[ifile]->Fill(phi,hitr);
           else if(modphi<5.0*3.14159/28. && hitr>sm_rmin && hitr<sm_rmax)
           h_sm_o[ifile]->Fill(phi,hitr);
           else if(modphi<7.0*3.14159/28. && hitr>sm_rmin && hitr<sm_rmax)
           h_sm_t[ifile]->Fill(phi,hitr);
           else if(hitr>sm_rmin && hitr<sm_rmax)
           h_sm_c[ifile]->Fill(phi,hitr);
         }
      }
   }
   }
  
   TString label[] = {"open","close","transition"};
   TLine* mainIn = new TLine(0,ring5_rmin,360,ring5_rmin);
   TLine* mainOut = new TLine(0,ring5_rmax,360,ring5_rmax);
   TLine* smIn = new TLine(0,sm_rmin,360,sm_rmin);
   TLine* smOut = new TLine(0,sm_rmax,360,sm_rmax);
   mainIn->SetLineColor(6);
   mainIn->SetLineWidth(2);
   mainIn->SetLineStyle(7);
   mainOut->SetLineColor(6);
   mainOut->SetLineWidth(2);
   mainOut->SetLineStyle(7);
   smIn->SetLineColor(6);
   smIn->SetLineWidth(2);
   smIn->SetLineStyle(7);
   smOut->SetLineColor(6);
   smOut->SetLineWidth(2);
   smOut->SetLineStyle(7);
   TCanvas* cM[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
   cM[ifile]=new TCanvas(Form("cM[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
   h_ring5[ifile]->Draw();
   h_ring5_o[ifile]->Draw("same");
   h_ring5_c[ifile]->Draw("same");
   h_ring5_t[ifile]->Draw("same");
   TLegend* legM = new TLegend(0.70,0.75,0.9,0.9);
   legM->SetBorderSize(0);
   legM->SetFillColor(0);
   legM->SetFillStyle(0);
   TLegendEntry* leg1[4];
   leg1[0] = legM->AddEntry(h_ring5_o[ifile],"total","p");
   leg1[1] = legM->AddEntry(h_ring5_o[ifile],"open","p");
   leg1[2] = legM->AddEntry(h_ring5_c[ifile],"close","p");
   leg1[3] = legM->AddEntry(h_ring5_t[ifile],"transition","p");
   legM->SetTextSize(0.05);
   leg1[0]->SetTextColor(color[ifile]);
   leg1[1]->SetTextColor(coloro[ifile]);
   leg1[2]->SetTextColor(colorc[ifile]);
   leg1[3]->SetTextColor(colort[ifile]);
   legM->Draw();
   mainIn->Draw();
   mainOut->Draw();
   cM[ifile]->SaveAs(Form("../temp/ring5_open_close_transition_%s.pdf",sim[ifile].Data()));
   }
   TCanvas* cMcolz[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
   cMcolz[ifile]=new TCanvas(Form("cMcolz[%d]",ifile),Form("Main %s",sim[ifile].Data()),600,600);
   h_ring5[ifile]->Draw("colz");
   h_ring5_o[ifile]->Draw("colz same");
   h_ring5_c[ifile]->Draw("colz same");
   h_ring5_t[ifile]->Draw("colz same");
   mainIn->Draw();
   mainOut->Draw();
   cMcolz[ifile]->SaveAs(Form("../temp/ring5_colz_%s.pdf",sim[ifile].Data()));
   }
   TCanvas* cS[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
   cS[ifile]=new TCanvas(Form("cS[%d]",ifile),Form("SM %s",sim[ifile].Data()),600,600);
   h_sm[ifile]->Draw();
   h_sm_o[ifile]->Draw("same");
   h_sm_c[ifile]->Draw("same");
   h_sm_t[ifile]->Draw("same");
   TLegend* legS = new TLegend(0.70,0.75,0.9,0.9);
   legS->SetBorderSize(0);
   legS->SetFillColor(0);
   legS->SetFillStyle(0);
   TLegendEntry* leg2[4];
   leg2[0] = legS->AddEntry(h_sm_o[ifile],"total","p");
   leg2[1] = legS->AddEntry(h_sm_o[ifile],"open","p");
   leg2[2] = legS->AddEntry(h_sm_c[ifile],"close","p");
   leg2[3] = legS->AddEntry(h_sm_t[ifile],"transition","p");
   legS->SetTextSize(0.05);
   leg2[0]->SetTextColor(color[ifile]);
   leg2[1]->SetTextColor(coloro[ifile]);
   leg2[2]->SetTextColor(colorc[ifile]);
   leg2[3]->SetTextColor(colort[ifile]);
   legS->Draw();
   smIn->Draw();
   smOut->Draw();
   cS[ifile]->SaveAs(Form("../temp/sm_open_close_transition_%s.pdf",sim[ifile].Data()));
   }
   TCanvas* cScolz[nfile];
   for(int ifile=0;ifile<nfile;ifile++){
   cScolz[ifile]=new TCanvas(Form("cScolz[%d]",ifile),Form("SM %s",sim[ifile].Data()),600,600);
   h_sm[ifile]->Draw("colz");
   h_sm_o[ifile]->Draw("colz same");
   h_sm_c[ifile]->Draw("colz same");
   h_sm_t[ifile]->Draw("colz same");
   smIn->Draw();
   smOut->Draw();
   cScolz[ifile]->SaveAs(Form("../temp/sm_colz_%s.pdf",sim[ifile].Data()));
   }
   gSystem->Exec(Form("pdfunite ../temp/ring5_* ../temp/sm_* ../plots/radial_vs_phi.pdf"));
   gSystem->Exec(Form("rm -rf ../temp/ring5_* ../temp/sm_*"));
}
