{
#include "TAxis.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include <stdio.h>
#include <stdlib.h>
#include <vector.h>

        double lumi=1.0; 


        TFile *fout = new TFile("out.root", "RECREATE");
 
        vector<TString> files;
        files.push_back("DY");
        files.push_back("TT");
        files.push_back("ST");
        files.push_back("WW");
        files.push_back("ZP2000");
        files.push_back("ZP2100");
        files.push_back("ZP2200");
        files.push_back("ZP2300");
        files.push_back("ZP2400");
        files.push_back("ZP2500");
        files.push_back("ZP2600");
        files.push_back("ZP2700");
        files.push_back("ZP2800");
        files.push_back("ZP2900");
        files.push_back("ZP3000");
        files.push_back("ZP3100");
        files.push_back("ZP3200");
        files.push_back("ZP3300");
        files.push_back("ZP3400");
        files.push_back("ZP3500");
        files.push_back("ZP3800");
        files.push_back("ZP4000");
        files.push_back("ZP4500");
        files.push_back("ZPspecial");

        

        const int nfiles = files.size();
        cout<<"file size "<<nfiles<<endl;

 

  for(int i=0; i<nfiles; i++)
  {
        TString filename = "../ntuple/ntuple_"+files.at(i)+".root";
        TString hname = "mll_"+files.at(i);
	TFile * input1 = new TFile (filename);
        input1->cd("");
        TTree *Outtree1 = (TTree*)input1->Get("tree");

        int nentries1 = (int)Outtree1->GetEntries();
        cout<<"Sample = "<<filename<<" "<<"nentries1 = "<<nentries1<<endl;
     
        double invmll, weight;

        Outtree1->SetBranchAddress("mll",&invmll);
        Outtree1->SetBranchAddress("SF",&weight);

        TH1F *PTA1   = new TH1F(hname," ",50,0,5000);
        PTA1->Sumw2();

        for(int j=0; j<nentries1; j++)
        {
         Outtree1->GetEntry(j);
         PTA1->Fill(invmll,lumi*weight);
        } 

 
         fout->cd(); 
         PTA1->Write();   


     }


        fout->Write(); 
        fout->Close();
 
 

}
