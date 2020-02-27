{
// cut couting based on TH1
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

     TH1D* hbkg= (TH1D*) mll_DY->Clone();
     hbkg->Add(h2);
     hbkg->Add(h3);
     hbkg->Add(h4);


     TH1D* hz2000= (TH1D*) mll_ZP2000->Clone();
     TH1D* hz2100= (TH1D*) mll_ZP2100->Clone();
     TH1D* hz2200= (TH1D*) mll_ZP2200->Clone();
     TH1D* hz2300= (TH1D*) mll_ZP2300->Clone();
     TH1D* hz2400= (TH1D*) mll_ZP2400->Clone();
     TH1D* hz2500= (TH1D*) mll_ZP2500->Clone();
     TH1D* hz2600= (TH1D*) mll_ZP2600->Clone();
     TH1D* hz2700= (TH1D*) mll_ZP2700->Clone();
     TH1D* hz2800= (TH1D*) mll_ZP2800->Clone();
     TH1D* hz2900= (TH1D*) mll_ZP2900->Clone();
     TH1D* hz3000= (TH1D*) mll_ZP3000->Clone();
     TH1D* hz3100= (TH1D*) mll_ZP3100->Clone();
     TH1D* hz3200= (TH1D*) mll_ZP3200->Clone();
     TH1D* hz3300= (TH1D*) mll_ZP3300->Clone();
     TH1D* hz3400= (TH1D*) mll_ZP3400->Clone();
     TH1D* hz3500= (TH1D*) mll_ZP3500->Clone();
     TH1D* hz3800= (TH1D*) mll_ZP3800->Clone();
     TH1D* hz4000= (TH1D*) mll_ZP4000->Clone();
     TH1D* hz4500= (TH1D*) mll_ZP4500->Clone();

     int inibin=18;
     double sigError;
     double sigCount = hz2000->IntegralAndError(inibin,50,sigError);
     double bkgError;
     double bkgCount = hbkg->IntegralAndError(inibin,50,bkgError);
     cout<<"mll>"<<0.0+(inibin-1)*100.0<<endl;
     cout<<"Signal="<<sigCount<<"   Bkg="<<bkgCount<<endl;
     double sig=TMath::Sqrt(2*((sigCount+bkgCount)*TMath::Log(1+(sigCount/(bkgCount)))-sigCount));
     cout<<"significance="<<sig<<endl; 
}
