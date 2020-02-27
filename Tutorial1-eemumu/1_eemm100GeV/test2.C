void test2(){

// Load shared library
 gSystem->Load("/home/qliphy/Desktop/MG5_aMC_v2_4_2/ExRootAnalysis/libExRootAnalysis.so");
 gSystem->Load("/home/qliphy/Desktop/common/root/lib/root/libPhysics.so");
// Create chain of root trees
 TChain chain("LHEF");
 TChain chain2("LHEF");

 TCanvas *c1= new TCanvas("c1","test graph",800,800);
 float inix= 0.0;
 float finx= 120.0;
 float nbin= 40.0;

 
TString IFile1 =   "1.root";
chain.Add(IFile1); 
// Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
 Long64_t numberOfEntries = treeReader->GetEntries();
// Get pointers to branches used in this analysis
 TClonesArray *branchEvent = treeReader->UseBranch("Event");
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
//  TClonesArray *branchJet = treeReader->UseBranch("Jet");
 TH1F *h1= new TH1F("1_Mjj","test histogram",nbin,inix,finx);
 h1->Sumw2();
 TLorentzVector J1, J2, JJ;
for(int i=0; i<=numberOfEntries-1; i++)
{
////*****************************************************************
    J1.SetPtEtaPhiE(0.,0.,0.,0.);
    J2.SetPtEtaPhiE(0.,0.,0.,0.);
    treeReader->ReadEntry(i);
    TRootLHEFEvent *event=(TRootLHEFEvent*) branchEvent->At(0);
    int np=event->Nparticles;
    for(int j=2; j<np; j++) {
     TRootLHEFParticle *particle1=(TRootLHEFParticle*) branchParticle->At(j);
      if(abs(particle1->PID)==13 ) {
          if(J1.E()>0){
          J2.SetPtEtaPhiE(particle1->PT,particle1->Eta,particle1->Phi,particle1->E);
             }
          else{  
          J1.SetPtEtaPhiE(particle1->PT,particle1->Eta,particle1->Phi,particle1->E);
             }
        }
     }
    JJ=J1+J2;
    h1->Fill(JJ.M(), 51.3/10000.); 
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
 TLorentzVector J12, J22, JJ2;
for(int i=0; i<=numberOfEntries2-1; i++)
{
////*****************************************************************
    J12.SetPtEtaPhiE(0.,0.,0.,0.);
    J22.SetPtEtaPhiE(0.,0.,0.,0.);
    treeReader2->ReadEntry(i);
    TRootLHEFEvent *event2=(TRootLHEFEvent*) branchEvent2->At(0);
    int np=event2->Nparticles;
    for(int j=2; j<np; j++) {
     TRootLHEFParticle *particle12=(TRootLHEFParticle*) branchParticle2->At(j);
      if(abs(particle12->PID)==13 ) {
          if(J12.E()>0){
          J22.SetPtEtaPhiE(particle12->PT,particle12->Eta,particle12->Phi,particle12->E);
             }
          else{  
          J12.SetPtEtaPhiE(particle12->PT,particle12->Eta,particle12->Phi,particle12->E);
             }
        }
     }
    JJ2=J12+J22;
    h2->Fill(JJ2.M(), 9.0/10000.); 
////*****************************************************************
} 
 

cout<< h2->GetBinContent(4)<<" "<<h2->GetBinError(4)<<endl;


      c1->SetLogy();
      h2->SetTitle("#mu#mu at 100GeV ILC");
      h2->GetXaxis()->SetTitle("M_{#mu#mu} [GeV] ");
      h2->GetYaxis()->SetTitle("a.u.");
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      h2->SetStats(kFALSE);
      h2->SetLineColor(kBlue);
      h2->SetLineWidth(3);
      h2->SetMarkerStyle(20);
      h2->GetXaxis()->SetTitleOffset(1.4);
      h2->GetYaxis()->SetTitleOffset(1.4);
      h2->Draw("HIST e");

      h1->SetLineColor(kRed);
      h1->SetLineWidth(3);
      h1->SetMarkerStyle(20);
      h1->Draw("HIST e same");
      



     TLegend *l1 = new TLegend(0.4,0.55,0.7,0.76);
     l1->SetBorderSize(1);
     l1->SetFillColor(0);
     l1->AddEntry(h1,"e e #rightarrow Z,#gamma  #rightarrow #mu#mu ");
     l1->AddEntry(h2,"e e #rightarrow #gamma #rightarrow #mu#mu ");
     l1->Draw();
     c1->SaveAs("Mmm.png");
 
 
}
