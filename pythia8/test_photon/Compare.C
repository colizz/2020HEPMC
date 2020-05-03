{
#include "TAxis.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector.h>



	TFile * input1 = new TFile ("mg267p8mlm_hepmc.root");
        input1->cd("");
        TTree *Outtree1 = (TTree*)input1->Get("h");

 	TFile * input2 = new TFile ("za.root");
        input2->cd("");
        TTree *Outtree2 = (TTree*)input2->Get("h");

	TFile * input3 = new TFile ("mg267p8mlm_hepmc_nomuphoton.root");
        input3->cd("");
        TTree *Outtree3 = (TTree*)input3->Get("h");

       int nentries1 = (int)Outtree1->GetEntries();
       cout<<"nentries1 = "<<nentries1<<endl;
       int nentries2 = (int)Outtree2->GetEntries();
       cout<<"nentries2 = "<<nentries2<<endl;
       int nentries3 = (int)Outtree3->GetEntries();
       cout<<"nentries3 = "<<nentries3<<endl; 

       TH1D *PTA1   = new TH1D("a"," ",50,0,100);
       PTA1->Sumw2();
       TH1D *PTA2   = new TH1D("b"," ",50,0,100);
       PTA2->Sumw2();
       TH1D *PTA3   = new TH1D("c"," ",50,0,100);
       PTA3->Sumw2();
  
       double ptz,pta;
       double ptz2,pta2;
       double ptz3,pta3;

       Outtree1->SetBranchAddress("ptz" ,&ptz);
       Outtree1->SetBranchAddress("pta" ,&pta);

       Outtree2->SetBranchAddress("ptz" ,&ptz2);
       Outtree2->SetBranchAddress("pta" ,&pta2);

       Outtree3->SetBranchAddress("ptz" ,&ptz3);
       Outtree3->SetBranchAddress("pta" ,&pta3);


       for(int i=0; i<nentries1; i++)
     {
         Outtree1->GetEntry(i);
         PTA1->Fill(pta,1.0*1476.0/float(nentries1));
     } 

       for(int i=0; i<nentries2; i++)
     {
         Outtree2->GetEntry(i);
         PTA2->Fill(pta2,1.0*3.538/float(nentries2));
     } 
 
       for(int i=0; i<nentries3; i++)
     {
         Outtree3->GetEntry(i);
         PTA3->Fill(pta3,1.0*1476.0/float(nentries3));
     } 

      TCanvas *c01 = new TCanvas("c01","",700,500);
      c01->SetLogy();
      PTA1->SetTitle("ZJets vs Z+Photon at 13TeV LHC");
      PTA1->GetXaxis()->SetTitle("Photon PT [GeV] ");
      PTA1->GetYaxis()->SetTitle("a.u.");
      PTA1->GetXaxis()->CenterTitle();
      PTA1->GetYaxis()->CenterTitle();
      PTA1->SetStats(kFALSE);
      PTA1->SetLineColor(kBlack);
      PTA1->SetLineWidth(3);
      PTA1->SetMarkerStyle(20);
      PTA1->GetXaxis()->SetTitleOffset(1.4);
      PTA1->Draw("e HIST");

      PTA2->SetLineColor(kRed);
      PTA2->SetLineWidth(2);
      PTA2->Draw("e same");
 
      PTA3->SetLineColor(kBlue);
      PTA3->SetLineWidth(2);
      PTA3->Draw("e same");
 

     TLegend *l1 = new TLegend(0.62,0.65,0.86,0.85);
     l1->SetBorderSize(1);
     l1->SetFillColor(0);
     l1->AddEntry(PTA1,"ZJets");
     l1->AddEntry(PTA2,"Zgamma");
     l1->AddEntry(PTA3,"ZJets, gamma not from muon"); 
     l1->Draw();


     c01->SaveAs("compare-ptll2.png");

}
