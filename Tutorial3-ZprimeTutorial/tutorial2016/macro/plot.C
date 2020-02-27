{
#include "TAxis.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector.h>


     TFile *file1 = new TFile("out.root","R"); 
     file1->cd("");

     TH1D* h1= (TH1D*) mll_DY->Clone();
     TH1D* h2= (TH1D*) mll_TT->Clone();
     TH1D* h3= (TH1D*) mll_ST->Clone();
     TH1D* h4= (TH1D*) mll_WW->Clone();
     h1->SetFillColor(kGreen);
     h2->SetFillColor(kYellow);
     h3->SetFillColor(kRed);
     h4->SetFillColor(kBlue);

     THStack * Mstack = new THStack("Mstack","");
     Mstack->Add(h3);
     Mstack->Add(h4);
     Mstack->Add(h2);
     Mstack->Add(h1);
     Mstack->SetMaximum(float(4.0)*Mstack->GetMaximum());

     TH1D* hz2000= (TH1D*) mll_ZP2000->Clone();
     TH1D* hz2500= (TH1D*) mll_ZP2500->Clone();
     TH1D* hz3000= (TH1D*) mll_ZP3000->Clone();
     TH1D* hz3500= (TH1D*) mll_ZP3500->Clone();


     TFile *file2 = new TFile("outdata.root","R"); 
     file2->cd("");
     TH1D* hdata= (TH1D*) mll_data->Clone();


     TCanvas *c01 = new TCanvas("c01","",700,500);
     c01->SetLogy();
     Mstack->Draw("HIST");
     Mstack->Draw("E same");
     Mstack->GetXaxis()->SetTitle("M_{#mu#mu} [GeV] ");
     Mstack->GetYaxis()->SetTitle("Events/bin");

     Mstack->GetXaxis()->CenterTitle();
     Mstack->GetYaxis()->CenterTitle();

     hdata->SetLineWidth(1);
     hdata->SetMarkerColor(4);
     hdata->SetMarkerStyle(20);
     hdata->Draw("EP, same");

     hz2000->SetLineWidth(1);
     hz2000->SetMarkerColor(6);
     hz2000->SetMarkerStyle(20);
     hz2000->Draw("EP, same");
     hz2500->SetLineWidth(1.);
     hz2500->SetMarkerColor(6);
     hz2500->SetMarkerStyle(22);
     hz2500->Draw("EP, same");
     hz3500->SetLineWidth(1);
     hz3500->SetMarkerColor(6);
     hz3500->SetMarkerStyle(25);
     hz3500->Draw("EP, same");



     TLegend *l1 = new TLegend(0.72,0.6,0.88,0.88);
     l1->SetBorderSize(1);
     l1->SetFillColor(0);
     l1->AddEntry(h1,"DY","f");
     l1->AddEntry(h2,"TT","f");
     l1->AddEntry(h3,"ST","f");
     l1->AddEntry(h4,"WW","f");
     l1->AddEntry(hdata,"data","EP");
     l1->AddEntry(hz2000,"Zprime 2000","EP");
     l1->AddEntry(hz2500,"Zprime 2500","EP");
     l1->AddEntry(hz3500,"Zprime 3500","EP");
     l1->Draw();


     TLatex *t = new TLatex();
     t->SetNDC();
     t->SetTextAlign(22);
     t->SetTextFont(64);
     t->SetTextSizePixels(17);
     t->DrawLatex(0.4,0.84,"ISTEP Preliminary 1. fb^{-1} (mu-channel) at #sqrt{s}=13TeV");
     t->Draw();
     c01->SaveAs("Mll-ISTEP.png");

}
