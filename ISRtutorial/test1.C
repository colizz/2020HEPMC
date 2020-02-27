void test1(){

      gStyle->SetPadBorderMode(0);
      gStyle->SetOptStat(0);
      gStyle->SetPadGridX(1);
      gStyle->SetPadGridY(1);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
      gStyle->SetPadTickX(1);
      gStyle->SetPadTickY(1);
      gStyle->SetAxisColor(1, "XYZ");
      gStyle->SetStripDecimals(kTRUE);
      gStyle->SetTickLength(0.03, "XYZ");
      gStyle->SetNdivisions(510, "XYZ");

// Load shared library
 gSystem->Load("/home/qliphy/Desktop/MG5_aMC_v2_4_2/ExRootAnalysis/libExRootAnalysis.so");
 gSystem->Load("/home/qliphy/Desktop/common/root/lib/root/libPhysics.so");
// Create chain of root trees
 TChain chain("LHEF");
 TChain chain2("LHEF");

 TCanvas *c1= new TCanvas("c1","test graph",800,800);
 float inix= 210.0;
 float finx= 260.0;
 float nbin= 50.0;

 
TString IFile1 =   "1.root";
chain.Add(IFile1); 
// Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
 Long64_t numberOfEntries = treeReader->GetEntries();
// Get pointers to branches used in this analysis
 TClonesArray *branchEvent = treeReader->UseBranch("Event");
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
//  TClonesArray *branchJet = treeReader->UseBranch("Jet");
 TH1F *h1= new TH1F("1_Mjj","",nbin,inix,finx);
 h1->Sumw2();
 double mee;
for(int i=0; i<=numberOfEntries-1; i++)
{
////*****************************************************************
    treeReader->ReadEntry(i);
    TRootLHEFParticle *particle1=(TRootLHEFParticle*) branchParticle->At(0);
    TRootLHEFParticle *particle2=(TRootLHEFParticle*) branchParticle->At(1);
    mee=sqrt((particle2->E+particle1->E)*(particle2->E+particle1->E)-(particle2->Pz+particle1->Pz)*(particle2->Pz+particle1->Pz));
    h1->Fill(mee,1.0/1000.); 
////*****************************************************************
}

 
 
TString IFile2 =   "2.root";
chain2.Add(IFile2); 
// Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader2 = new ExRootTreeReader(&chain2);
 Long64_t numberOfEntries2 = treeReader2->GetEntries();
// Get pointers to branches used in this analysis
 TClonesArray *branchEvent2 = treeReader2->UseBranch("Event");
 TClonesArray *branchParticle2 = treeReader2->UseBranch("Particle");
//  TClonesArray *branchJet = treeReader->UseBranch("Jet");
 TH1F *h2= new TH1F("2_Mjj","test histogram",nbin,inix,finx);
 h2->Sumw2();
 double mee2;
for(int i=0; i<=numberOfEntries2-1; i++)
{
////*****************************************************************
    treeReader2->ReadEntry(i);
    TRootLHEFParticle *particle12=(TRootLHEFParticle*) branchParticle2->At(0);
    TRootLHEFParticle *particle22=(TRootLHEFParticle*) branchParticle2->At(1);
    mee2=sqrt((particle22->E+particle12->E)*(particle22->E+particle12->E)-(particle22->Pz+particle12->Pz)*(particle22->Pz+particle12->Pz));
    h2->Fill(mee2, 1.0/1000.); ////*****************************************************************
}
 

/*      TPad top_pad("top_pad", "top_pad",0,0.2, 1.0, 1.0);
      top_pad.Draw();
      top_pad.cd();
      top_pad->SetLogy();

      top_pad.SetBottomMargin(0);
*/     
//      h1->SetTitle("p p #rightarrow W #gamma ,  W #rightarrow l #nu at 13TeV LHC");
      c1->SetLogy();
      h1->GetXaxis()->SetTitle("Center-of-mass energy (GeV)");
      h1->GetYaxis()->SetTitle("a.u.");
      h1->GetXaxis()->CenterTitle();
      h1->GetYaxis()->CenterTitle();
      h1->SetStats(kFALSE);
      h1->SetLineColor(kBlue);
      h1->SetLineWidth(2);
      h1->SetMarkerStyle(24);
      h1->SetMarkerSize(0.8);
      h1->SetLineStyle(2);
      h1->GetXaxis()->SetTitleOffset(1.4);
      h1->GetYaxis()->SetTitleOffset(1.4);
      h1->SetMinimum(0.0001);
      h1->Draw("e HIST");

      h1->SetMaximum(0.8);
      h1->SetMinimum(0.0001);
 
      h2->SetLineColor(kRed);
      h2->SetLineWidth(2);
      h2->SetMarkerStyle(22);
      h2->SetMarkerColor(kRed);
      h2->Draw("e HIST same");

      TLegend *l1 = new TLegend(0.2,0.52,0.57,0.68);
      l1->SetBorderSize(1);
      l1->SetFillColor(0);
      l1->AddEntry(h1,"Whizard");
      l1->AddEntry(h2,"MG with ISR");
      l1->Draw();

      TLatex latex2;
      latex2.SetNDC();
      latex2.SetTextAngle(0);
      latex2.SetTextColor(kBlack);
      latex2.SetTextFont(52);
      latex2.SetTextSize(0.04);
      latex2.SetTextAlign(31);
      latex2.DrawLatex(0.66, 0.83, "e^{+}e^{-}#rightarrow W^{+}W^{-}#gamma @ 240GeV CEPC");

/* TH1F *h4=  (TH1F*) h1 -> Clone();
 h4->Sumw2(); 
 h4->Divide(h2);
 h4->Scale(100.0);

      TPad bottom_pad("bottom_pad", "bottom_pad", 0, 0., 1.0, 0.2);
      bottom_pad.Draw();
      bottom_pad.cd();
      bottom_pad.SetTopMargin(0);
      h4->SetStats(kFALSE);
      h4->SetTitle("(Full-EWK-QCD)/EWK (%)");
      h4->SetTitleSize(0.4);
      h4->GetXaxis()->SetLabelSize(0.1);
      h4->GetXaxis()->SetTitleSize(0.1);
      h4->Draw();
*/
      c1->SaveAs("COM.eps");
 
 
}
