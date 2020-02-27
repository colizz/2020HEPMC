{
// Load shared library
 gSystem->Load("/home/qliphy/Desktop/MG5_aMC_v2_4_2/ExRootAnalysis/libExRootAnalysis.so");
 gSystem->Load("/home/qliphy/Desktop/common/root/lib/root/libPhysics.so");
// Create chain of root trees
 TChain chain("LHEF");

 chain.Add("y.root"); 

// Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
 Long64_t numberOfEntries = treeReader->GetEntries();
// Get pointers to branches used in this analysis
 TClonesArray *branchEvent = treeReader->UseBranch("Event");
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
//  TClonesArray *branchJet = treeReader->UseBranch("Jet");
 TCanvas *c1= new TCanvas("c1","test graph",800,800);

 float inix= -1.0;
 float finx= 1.;
 float nbin= 100.0; 
 TH1F *h1= new TH1F("h1"," ",nbin,inix,finx);
      h1->Sumw2();
  float mwp,mwm,phiwp,phiwm,phiww;
  float iwpd1,iwpd2;
  float iwmd1,iwmd2;
  float ptl,etal,phil,el;
  float ptv,etav,phiv,ev;
  float phi1,phi2,phitemp;
for(int i=0; i<=numberOfEntries-1; i++)
{
////*****************************************************************
    treeReader->ReadEntry(i);
    TRootLHEFEvent *event=(TRootLHEFEvent*) branchEvent->At(0);
    TRootLHEFParticle *particle1=(TRootLHEFParticle*) branchParticle->At(0);
    TRootLHEFParticle *particle2=(TRootLHEFParticle*) branchParticle->At(1);
 
//---MG does not output Mgraviton, not like JHU, to check--------------------

      int nne=event->Nparticles;
      int ifw,ifwb,pidwp,pidwm;
      int firstdwp=1;
      int firstdwm=1;

	TLorentzVector w;
	TLorentzVector lep;
	TLorentzVector neu;

//-------------------
      for(int j=3; j<=nne; j++){
      TRootLHEFParticle *ppw=(TRootLHEFParticle*) branchParticle->At(j-1);

      ifw=(fabs(ppw->PID)-13)*(fabs(ppw->PID)-11);
      if(ifw==0) {ptl=ppw->PT; etal=ppw->Eta; phil=ppw->Phi; el=ppw->E;}
      ifwb=(fabs(ppw->PID)-12)*(fabs(ppw->PID)-14);
      if(ifwb==0) {ptv=ppw->PT; etav=ppw->Eta; phiv=ppw->Phi;ev=ppw->E;}
       } 

       lep.SetPtEtaPhiE(ptl,etal,phil,el);
       neu.SetPtEtaPhiE(ptv,etav,phiv,ev);
       w=lep+neu;
       TVector3 v1 = w.Vect();
  	TLorentzVector lep_in_w(lep);
	TVector3 wboost = -(w.BoostVector());
	lep_in_w.Boost(wboost);
	TVector3 v2 = lep_in_w.Vect();
 	double theta = v1.Angle(v2);

//-------------------
      
     h1->Fill(cos(theta),nbin/(finx-inix)/float(numberOfEntries)); 


////*****************************************************************
}


//(((((((((((((((((((((((((((((((((((((((((((((((((((((((
 TChain chain2("LHEF");
 chain2.Add("h.root"); 

// Create object of class ExRootTreeReader
 ExRootTreeReader *treeReader2 = new ExRootTreeReader(&chain2);
 Long64_t numberOfEntries2 = treeReader2->GetEntries();
// Get pointers to branches used in this analysis
 TClonesArray *branchEvent2 = treeReader2->UseBranch("Event");
 TClonesArray *branchParticle2 = treeReader2->UseBranch("Particle");
 

  TH1F *h2= new TH1F("h2"," ",nbin,inix,finx);
  h2->Sumw2();
  float mwp2,mwm2,phiwp2,phiwm2,phiww2;
  float iwpd12,iwpd22;
  float iwmd12,iwmd22;
  float ptl2,etal2, phil2,el2;
  float ptv2,etav2,phiv2,ev2;
  float phi12,phi22,phitemp2;
for(int i2=0; i2<=numberOfEntries2-1; i2++)
{
////*****************************************************************
    treeReader2->ReadEntry(i2);
    TRootLHEFEvent *event2=(TRootLHEFEvent*) branchEvent2->At(0);
    TRootLHEFParticle *particle12=(TRootLHEFParticle*) branchParticle2->At(0);
    TRootLHEFParticle *particle22=(TRootLHEFParticle*) branchParticle2->At(1);
 
//---MG does not output Mgraviton, not like JHU, to check--------------------

      int nne2=event2->Nparticles;
      int ifw2,ifwb2,pidwp2,pidwm2;
      int firstdwp2=1;
      int firstdwm2=1;


	TLorentzVector w2;
	TLorentzVector lep2;
	TLorentzVector neu2;
//-------------------
      for(int j2=3; j2<=nne2; j2++){
      TRootLHEFParticle *ppw2=(TRootLHEFParticle*) branchParticle2->At(j2-1);

      ifw2=(fabs(ppw2->PID)-13)*(fabs(ppw2->PID)-11);
      if(ifw2==0) {ptl2=ppw2->PT; etal2=ppw2->Eta; phil2=ppw2->Phi; el2=ppw2->E;}
      ifwb2=(fabs(ppw2->PID)-12)*(fabs(ppw2->PID)-14);
      if(ifwb2==0) {ptv2=ppw2->PT; etav2=ppw2->Eta; phiv2=ppw2->Phi;ev2=ppw2->E;}
       } 

       lep2.SetPtEtaPhiE(ptl2,etal2,phil2,el2);
       neu2.SetPtEtaPhiE(ptv2,etav2,phiv2,ev2);
       w2=lep2+neu2;
       TVector3 v12 = w2.Vect();
  	TLorentzVector lep_in_w2(lep2);
	TVector3 wboost2 = -(w2.BoostVector());
	lep_in_w2.Boost(wboost2);
	TVector3 v22 = lep_in_w2.Vect();
 	double theta2 = v12.Angle(v22);
 
//-------------------
      
     h2->Fill(cos(theta2),nbin/(finx-inix)/float(numberOfEntries2)); 

////*****************************************************************
}
//)))))))))))))))))))))))))))))))))))))))))))))))))))))))



h1->GetXaxis()->SetTitle("cos#theta_{l}  ");
h1->GetYaxis()->SetTitle("dcos#sigma/d#theta_{l} (Events)");
h1->GetXaxis()->CenterTitle();
h1->GetYaxis()->CenterTitle();

h1->SetLineColor(kBlack);
h2->SetLineColor(kRed);
//h3->SetLineColor(kGreen);
//h4->SetLineColor(kBlue);

 h1->Draw("HIST");
 h2->Draw("HIST same");
TLegend *my =new TLegend(0.55,0.6,0.9,0.8,"");
my->SetTextSize(0.02);
my->AddEntry(h1,"MG","l");
my->AddEntry(h2,"JHU","l");
my->Draw();
 c1->Update();

 

}
