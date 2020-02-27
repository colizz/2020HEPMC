{
#include "TAxis.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector.h>



	TFile * input1 = new TFile ("official.root");
        input1->cd("");
        TTree *Outtree1 = (TTree*)input1->Get("MCPart");

	TFile * input2 = new TFile ("signal.root");
        input2->cd("");
        TTree *Outtree2 = (TTree*)input2->Get("MCPart");
 

       int nentries1 = (int)Outtree1->GetEntries();
       cout<<"nentries1 = "<<nentries1<<endl;
       int nentries2 = (int)Outtree2->GetEntries();
       cout<<"nentries2 = "<<nentries2<<endl;
 
       TH1F *PTA1   = new TH1F("a"," ",20,0,400);
       TH1F *PTA2   = new TH1F("b"," ",20,0,400);

       float RecoilMass1,RecoilMass2;

       Outtree1->SetBranchAddress("RecoilMass" ,&RecoilMass1);
       Outtree2->SetBranchAddress("RecoilMass" ,&RecoilMass2);

       for(int i=0; i<nentries1; i++)
     {
         Outtree1->GetEntry(i);
         PTA1->Fill(RecoilMass1,1.0/float(nentries1));
     } 

       for(int i=0; i<nentries2; i++)
     {
         Outtree2->GetEntry(i);
         PTA2->Fill(RecoilMass2,1.0/float(nentries2));
     } 
 


      TCanvas *c01 = new TCanvas("c01","",700,500);
      c01->SetLogy();
      PTA1->SetTitle("e^{+}e^{-} #rightarrow #mu^{+}#mu^{-} H at CEPC");
      PTA1->GetXaxis()->SetTitle("Recoil Mass [GeV] ");
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
 
      
 

     TLegend *l1 = new TLegend(0.52,0.6,0.76,0.8);
     l1->SetBorderSize(1);
     l1->SetFillColor(0);
     l1->AddEntry(PTA1,"Official Whizard+Pythia6");
     l1->AddEntry(PTA2,"MG+Pythia8");
     l1->Draw();


     c01->SaveAs("compare-ptll2.png");

}
