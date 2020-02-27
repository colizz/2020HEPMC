/*
This macro is used to produce the ntuple.root,
storing useful informations we need
//root -b -q ntuple_chain.C\(\"WW\"\)  
*/

//------------------------------------------------------------------------------

double deltaPhi(const double& phi1, const double& phi2)
{
        double deltaphi = fabs(phi1 - phi2);
        if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
        return deltaphi;
}

double deltaEta(const double& eta1, const double& eta2)
{
        double deltaeta = fabs(eta1 - eta2);
        return deltaeta;
}

double deltaR(const double& eta1, const double& phi1,
                const double& eta2, const double& phi2)
{
        double deltaphi = deltaPhi(phi1, phi2);
        double deltaeta = deltaEta(eta1, eta2);
        double deltar = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);
        return deltar;
}

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TString filename)
{
  // create a root file to store the variables
  TFile file(filename, "recreate");
  TTree *tree = new TTree("tree", "keep events after cuts");

  double SF=1.0;
  if(filename.Contains("DY")==1) {SF=7.48/80000.;}
  if(filename.Contains("TT")==1) {SF=0.76/40000.;}
  if(filename.Contains("ST")==1) {SF=0.089/10000.;}
  if(filename.Contains("WW")==1) {SF=0.338/40000.;}
  if(filename.Contains("ZP2000")==1) {SF=15.58/10000.;}
  if(filename.Contains("ZP2100")==1) {SF=12.06/10000.;}
  if(filename.Contains("ZP2200")==1) {SF=9.38/10000.;}
  if(filename.Contains("ZP2300")==1) {SF=7.33/10000.;}
  if(filename.Contains("ZP2400")==1) {SF=5.58/10000.;}
  if(filename.Contains("ZP2500")==1) {SF=4.55/10000.;}
  if(filename.Contains("ZP2600")==1) {SF=3.59/10000.;}
  if(filename.Contains("ZP2700")==1) {SF=2.86/10000.;}
  if(filename.Contains("ZP2800")==1) {SF=2.29/10000.;}
  if(filename.Contains("ZP2900")==1) {SF=1.83/10000.;}
  if(filename.Contains("ZP3000")==1) {SF=1.47/10000.;}
  if(filename.Contains("ZP3100")==1) {SF=1.18/10000.;}
  if(filename.Contains("ZP3200")==1) {SF=0.957/10000.;}
  if(filename.Contains("ZP3300")==1) {SF=0.772/10000.;}
  if(filename.Contains("ZP3400")==1) {SF=0.626/10000.;}
  if(filename.Contains("ZP3500")==1) {SF=0.512/10000.;}
  if(filename.Contains("ZP3800")==1) {SF=0.280/10000.;}
  if(filename.Contains("ZP4000")==1) {SF=0.189/10000.;}
  if(filename.Contains("ZP4500")==1) {SF=0.0752/10000.;}
  if(filename.Contains("ZPspecial")==1) {SF=1.263/10000.;}


  int inputNum=0;
  int event;
  int numJettot, numJet, numbJet;
  double ptJet[99],etaJet[99],phiJet[99],massJet[99];
  double ptBJet[99],etaBJet[99],phiBJet[99],massBJet[99];
  double met,mll,ptll;
  int numEl, numMu, numLep;
  double ptMu[99],etaMu[99],phiMu[99],chargeMu[99];
   
  tree->Branch("event", &event, "event/I");
  tree->Branch("SF", &SF, "SF/D");
  tree->Branch("inputNum", &inputNum, "inputNum/I");
  tree->Branch("numJettot", &numJettot, "numJettot/I");
  tree->Branch("numJet", &numJet, "numJet/I");
  tree->Branch("numbJet", &numbJet, "numbJet/I");
  tree->Branch("numEl", &numEl, "numEl/I");
  tree->Branch("numMu", &numMu, "numMu/I");
  tree->Branch("numLep", &numLep, "numLep/I");
  tree->Branch("ptJet", &ptJet, "ptJet[numJet]/D");
  tree->Branch("etaJet", &etaJet, "etaJet[numJet]/D");
  tree->Branch("phiJet", &phiJet, "phiJet[numJet]/D");
  tree->Branch("massJet", &massJet, "massJet[numJet]/D");
  tree->Branch("ptBJet", &ptBJet, "ptBJet[numbJet]/D");
  tree->Branch("etaBJet", &etaBJet, "etaBJet[numbJet]/D");
  tree->Branch("phiBJet", &phiBJet, "phiBJet[numbJet]/D");
  tree->Branch("massBJet", &massBJet, "massBJet[numbJet]/D");
  tree->Branch("met", &met, "met/D");
  tree->Branch("mll", &mll, "mll/D");
  tree->Branch("ptll", &ptll, "ptll/D");
  tree->Branch("ptMu", &ptMu, "ptMu[numMu]/D");
  tree->Branch("etaMu", &etaMu, "etaMu[numMu]/D");
  tree->Branch("phiMu", &phiMu, "phiMu[numMu]/D");
  tree->Branch("chargeMu", &chargeMu, "chargeMu[numMu]/I");

  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  Long64_t numberOfEntries = treeReader->GetEntries();
  inputNum = numberOfEntries;

  cout << "** Chain contains " << numberOfEntries << " events" << endl;

  for(int count=0; count < numberOfEntries; count++)
  {
      //initiate branch members
      event=-1; numJettot=-1; numJet=0; numbJet=0;
      met=-99.0; mll=-99.0; ptll=-99.0; 
      numEl=-1; numMu=-1; 
      for(int i=0; i<99; i++)
      {
          ptJet[i]=-99;etaJet[i]=-99;phiJet[i]=-99;massJet[i]=-99;
          ptBJet[i]=-99;etaBJet[i]=-99;phiBJet[i]=-99;massBJet[i]=-99;
          ptMu[i]=-99;etaMu[i]=-99;phiMu[i]=-99;chargeMu[i]=-99;
      }

      treeReader->ReadEntry(count) ;
      //********************************************************************


      Event* ev = (Event*)branchEvent->At(0);
      event=ev->Number;
      numJettot = branchJet->GetEntries();
      MissingET* Met = (MissingET *) branchMET->At(0);
      met=Met->MET;
      numEl=branchElectron->GetEntries();
      numMu=branchMuon->GetEntries();

// Preselection Cuts1:  veto additional soft lepton
      if(numMu!=2)continue;
      if(numEl>0)continue;

      Jet* pJet[1200]; // access to Jet Collection
      double btagJet[1200]; // btag value of jet in Jet Collection
      for(int i=0, j=0, k=0; i<numJettot; i++)
      {
          if(branchJet->At(i))
       	  {
	      pJet[i]=(Jet*)branchJet->At(i);
              btagJet[i]=pJet[i]->BTag;
              if( pJet[i]->PT>30 && fabs(pJet[i]->Eta)<4.7 )
	      {
		  ptJet[j]=pJet[i]->PT;
                  etaJet[j]=pJet[i]->Eta;
                  phiJet[j]=pJet[i]->Phi;
                  massJet[j]=pJet[i]->Mass;
		  numJet++;j++;
	      }
              if( (btagJet[i] & (1 << 0)) && pJet[i]->PT>30 && fabs(pJet[i]->Eta)<2.5 )
	      {
		  ptBJet[k]=pJet[i]->PT;
                  etaBJet[k]=pJet[i]->Eta;
                  phiBJet[k]=pJet[i]->Phi;
                  massBJet[k]=pJet[i]->Mass;
		  numbJet++;k++;
	      }
	  }
      }  

      //********************************************************************
 
      Muon* M[99] ;
      int mark2=0;
 
      for(int j=0; j<numMu; j++)
      {
         M[j]=(Muon*)branchMuon->At(j);
         if(M[j]->PT>40 && fabs(M[j]->Eta)<2.4)
         {
             ptMu[mark2]=M[j]->PT;
             etaMu[mark2]=M[j]->Eta;
             phiMu[mark2]=M[j]->Phi;
             chargeMu[mark2]=M[j]->Charge;
             mark2++;
   	 }
      }

// Preselection Cuts1:  2 and only 2 good muons
      if(mark2!=2)continue;

      float MMass = 0.105000 ;
      TLorentzVector l1,l2 ;
      l1.SetPtEtaPhiM(ptMu[0],etaMu[0],phiMu[0],MMass);
      l2.SetPtEtaPhiM(ptMu[1],etaMu[1],phiMu[1],MMass);
      mll=(l1+l2).M();
      ptll=(l1+l2).Pt();


      tree->Fill();
  } // end of chain loop
  //tree->Print();
  file.Write();
}

//------------------------------------------------------------------------------

void ntuple_chain(char* dataFileName)
{
  gSystem->Load("/home/qliphy/Desktop/MG5_aMC_v2_3_0/Delphes/libDelphes");
  //THStack *stack;
  TChain *chain = new TChain("Delphes");

  TString buffer = dataFileName;
  cout<<buffer<<endl;
  TString IFile =  "./data/"+buffer + ".root";
  TString OFile =  "./ntuple/ntuple_" + buffer+".root";

  chain->Add(IFile);
  ExRootTreeReader *treeReader1 = new ExRootTreeReader(chain);
  ExRootResult *result1 = new ExRootResult();
  AnalyseEvents(treeReader1, OFile);

  cout << "** Exiting..." << endl;

  delete result1;
  delete treeReader1;
  delete chain;
}

//------------------------------------------------------------------------------
