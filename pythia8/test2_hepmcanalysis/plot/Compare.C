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

	TFile * input2 = new TFile ("mg267p8mlm.root");
        input2->cd("");
        TTree *Outtree2 = (TTree*)input2->Get("t");


       int nentries1 = (int)Outtree1->GetEntries();
       cout<<"nentries1 = "<<nentries1<<endl;
       int nentries2 = (int)Outtree2->GetEntries();
       cout<<"nentries2 = "<<nentries2<<endl;
    

       TH1D *PTA1   = new TH1D("a"," ",20,0,200);
       PTA1->Sumw2();
       TH1D *PTA2   = new TH1D("b"," ",20,0,200);
       PTA2->Sumw2();
 
 
       double ptz1,ptz2;

       Outtree1->SetBranchAddress("ptz" ,&ptz1);
       Outtree2->SetBranchAddress("ptz" ,&ptz2);
 

       for(int i=0; i<nentries1; i++)
     {
         Outtree1->GetEntry(i);
         PTA1->Fill(ptz1,1.0/float(nentries1));
     } 

       for(int i=0; i<nentries2; i++)
     {
         Outtree2->GetEntry(i);
         PTA2->Fill(ptz2,1.0/float(nentries2));
     } 
 
 

      TCanvas *c01 = new TCanvas("c01","",700,500);
      c01->SetLogy();
      PTA1->SetTitle("pp #rightarrow z #rightarrow #mu^{+}#mu^{-} 012Jets MLM at 13TeV LHC");
      PTA1->GetXaxis()->SetTitle("Z PT [GeV] ");
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
 

     TLegend *l1 = new TLegend(0.62,0.65,0.86,0.85);
     l1->SetBorderSize(1);
     l1->SetFillColor(0);
     l1->AddEntry(PTA1,"MG267_LHE_Pythia8 MLM");
     l1->AddEntry(PTA2,"MG627_HEP_Analysis");

     l1->Draw();


     c01->SaveAs("compare-ptll2.png");

}
