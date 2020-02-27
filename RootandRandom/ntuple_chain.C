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

  int inputNum=0;
  int event;
  double mll,ptll;
  int numEl, numMu;
  double ptEl[99],etaEl[99],phiEl[99],chargeEl[99];
  double genptEl[99],genetaEl[99],genphiEl[99],genchargeEl[99];
  double reso[99];
 
  tree->Branch("event", &event, "event/I");
  tree->Branch("inputNum", &inputNum, "inputNum/I");
  tree->Branch("numEl", &numEl, "numEl/I");
  tree->Branch("numMu", &numMu, "numMu/I");
  tree->Branch("mll", &mll, "mll/D");
  tree->Branch("ptll", &ptll, "ptll/D");
  tree->Branch("ptEl", &ptEl, "ptEl[numEl]/D");
  tree->Branch("etaEl", &etaEl, "etaEl[numEl]/D");
  tree->Branch("phiEl", &phiEl, "phiEl[numEl]/D");
  tree->Branch("chargeEl", &chargeEl, "chargeEl[numEl]/D");
  tree->Branch("genptEl", &genptEl, "genptEl[numEl]/D");
  tree->Branch("genetaEl", &genetaEl, "genetaEl[numEl]/D");
  tree->Branch("genphiEl", &genphiEl, "genphiEl[numEl]/D");
  tree->Branch("genchargeEl", &genchargeEl, "genchargeEl[numEl]/D");
  tree->Branch("reso", &reso, "reso[numEl]/D");

  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  //TClonesArray *branchJet = treeReader->UseBranch("Jet");
  //TClonesArray *branchMET = treeReader->UseBranch("MissingET");

  Long64_t numberOfEntries = treeReader->GetEntries();
  inputNum = numberOfEntries;

  cout << "** Chain contains " << numberOfEntries << " events" << endl;

  for(int count=0; count < numberOfEntries; count++)
  {
      //initiate branch members
      event=-1; mll=-99.0; ptll=-99.0; 
      numEl=-1; numMu=-1; 
      for(int i=0; i<99; i++)
      {
          ptEl[i]=-99;etaEl[i]=-99;phiEl[i]=-99;chargeEl[i]=-99;
          genptEl[i]=-99;genetaEl[i]=-99;genphiEl[i]=-99;genchargeEl[i]=-99;
          reso[i]=-1.;
      }

      treeReader->ReadEntry(count) ;
      //********************************************************************

      Event* ev = (Event*)branchEvent->At(0);
      event=ev->Number;
      numEl=branchElectron->GetEntries();
      numMu=branchMuon->GetEntries();

      GenParticle *particle;

// Preselection Cuts:
      if(numEl<2)continue;

      Electron* M[99] ;
      int mark2=0;
 
      for(int j=0; j<numEl; j++)
      {
         M[j]=(Electron*)branchElectron->At(j);
         if(M[j]->PT>1. && fabs(M[j]->Eta)<4)
         {
             ptEl[mark2]=M[j]->PT;
             etaEl[mark2]=M[j]->Eta;
             phiEl[mark2]=M[j]->Phi;
             chargeEl[mark2]=M[j]->Charge;

             particle = (GenParticle*) M[j]->Particle.GetObject();
             genptEl[mark2]=particle->PT;
             genetaEl[mark2]=particle->Eta;
             genphiEl[mark2]=particle->Phi;
             genchargeEl[mark2]=particle->Charge;

             reso[mark2]=(ptEl[mark2]-genptEl[mark2])/genptEl[mark2];

             mark2++;

   	 }
      }

      float MMass = 0.0 ;
      TLorentzVector l1,l2 ;
      l1.SetPtEtaPhiM(ptEl[0],etaEl[0],phiEl[0],MMass);
      l2.SetPtEtaPhiM(ptEl[1],etaEl[1],phiEl[1],MMass);
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
  gSystem->Load("/home/qliphy/Desktop/CEPC-Delphes/delphes-master/libDelphes");
  //THStack *stack;
  TChain *chain = new TChain("Delphes");

  TString buffer = dataFileName;
  cout<<buffer<<endl;
  TString IFile =  buffer + ".root";
  TString OFile =  "ntuple_" + buffer+".root";

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
