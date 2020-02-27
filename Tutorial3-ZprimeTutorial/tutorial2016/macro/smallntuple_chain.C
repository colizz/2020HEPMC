void smallntuple_chain(char* dataFileName)
{
  TString buffer = dataFileName;
  cout<<buffer<<endl;

  TString IFile =  "../ntuple/ntuple_" + buffer+".root";
  TFile *ifile = new TFile(IFile);
  TTree *tree = (TTree*)ifile->Get("tree");
  int nentries = tree->GetEntries();

  int numJet, numbJet, numMu;
  double met,mll,ptll,SF;
  double ptMu[2],etaMu[2],phiMu[2];
 

  tree->SetBranchAddress("numJet", &numJet);
  tree->SetBranchAddress("numbJet", &numbJet);
  tree->SetBranchAddress("numMu", &numMu);
  tree->SetBranchAddress("met", &met);
  tree->SetBranchAddress("mll", &mll);
  tree->SetBranchAddress("ptll", &ptll);
  tree->SetBranchAddress("ptMu", ptMu);
  tree->SetBranchAddress("etaMu", etaMu);
  tree->SetBranchAddress("phiMu", phiMu);
  tree->SetBranchAddress("SF",&SF);



  TString OFile =  "small"+buffer+".root";
  TFile  *ofile=new TFile(OFile, "recreate");
  TTree *otree = new TTree(buffer, buffer);
  
  double ptMu1,ptMu2,etaMu1,etaMu2,phiMu1,phiMu2;
  double detall, dphill;

  otree->Branch("numJet", &numJet, "numJet/I");
  otree->Branch("numbJet", &numbJet, "numbJet/I");
  otree->Branch("mll", &mll, "mll/D");
  otree->Branch("ptll", &ptll, "ptll/D");
  otree->Branch("SF", &SF, "SF/D");
  otree->Branch("ptMu1", &ptMu1, "ptMu1/D");
  otree->Branch("ptMu2", &ptMu2, "ptMu2/D");
  otree->Branch("etaMu1", &etaMu1, "etaMu1/D");
  otree->Branch("etaMu2", &etaMu2, "etaMu2/D");
  otree->Branch("phiMu1", &phiMu1, "phiMu1/D");
  otree->Branch("phiMu2", &phiMu2, "phiMu2/D");
  otree->Branch("detall", &detall, "detall/D");
  otree->Branch("dphill", &dphill, "dphill/D");

    for(int j=0; j<nentries; j++)
    {
      tree->GetEntry(j);
      //std::cout<<numMu<<std::endl;
      if(numMu!=2) continue;
      
      int ii=0, jj=1;
      if(ptMu[0]<ptMu[1]) {ii=1,jj=0;}
      ptMu1=ptMu[ii];
      etaMu1=etaMu[ii];
      phiMu1=phiMu[ii];
      ptMu2=ptMu[jj];
      etaMu2=etaMu[jj];
      phiMu2=phiMu[jj];
      detall=fabs(etaMu1-etaMu2);
      dphill=deltaPhi(phiMu1,phiMu2);
      otree->Fill();
    }

   ofile->Write();

}

double deltaPhi(const double& phi1, const double& phi2)
{
        double deltaphi = fabs(phi1 - phi2);
        if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
        return deltaphi;
}
