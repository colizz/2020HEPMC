void ana(){
 	gSystem->Load("/home/jing/jing/root5/root/lib/libPhysics.so");

	TCanvas *c1= new TCanvas("c1","test graph",800,800);
	float inix= 0.0;
	float finx=100.0;
	float nbin= 50.0;
	//input file
   	TFile *f = new TFile("../20190110/workspace/h13_10k_mc.root");
	TTree *t1 = (TTree*)f->Get("h");
   	double pt, m, px, py, pz, gm, theta, eta, phi, e ;
   	int event, pid, status, scale;
        t1->SetBranchAddress("scale",&scale);
   	t1->SetBranchAddress("px",&px);
   	t1->SetBranchAddress("py",&py);
   	t1->SetBranchAddress("pz",&pz);
   	t1->SetBranchAddress("m",&m);
        t1->SetBranchAddress("pt",&pt);
        t1->SetBranchAddress("gm",&gm);
        t1->SetBranchAddress("theta",&theta);
        t1->SetBranchAddress("phi",&phi);
   	t1->SetBranchAddress("e",&e);
        t1->SetBranchAddress("status",&status);
        t1->SetBranchAddress("pid",&pid);
        t1->SetBranchAddress("event",&event);

	//define histogram
	TH1F *h1= new TH1F("1","test histogram",nbin,inix,finx);
        h1->Sumw2();
	//define lorentz vector
	TLorentzVector J1, J2, JJ;
	J1.SetPtEtaPhiE(0.,0.,0.,0.);
	J2.SetPtEtaPhiE(0.,0.,0.,0.);
   	//read all entries and fill the histograms
   	Long64_t nentries = t1->GetEntries();
	int N=0, h, j;
	Long64_t i=0;
	cout<<"ok"<<endl;
   	for (Long64_t n=0;n<10000;n++) {
		cout<<"event= "<<n<<endl;
        	J1.SetPtEtaPhiE(0.,0.,0.,0.);
        	J2.SetPtEtaPhiE(0.,0.,0.,0.);
		h = 0;
		j=0;
		cout<<N<<" "<<event<<endl;
     		while(N==event){
	                t1->GetEntry(i);
			//cout<<"pid= "<<pid<<endl;
                        if(abs(pid)==25) {h = h + 1; i++;  continue;}
			if(scale!=250202402) {i++;continue;}
			if(abs(pid)==15||abs(pid)==13||abs(pid)==11){
          			if(J1.E()>0){
          				J2.SetPtEtaPhiE(pt,eta,phi,e);
					j++;
					//cout<<pid<<" "<<j<<endl;
             			}
          			else{
          				J1.SetPtEtaPhiE(pt,eta,phi,e);
					j++;
                                        //cout<<pid<<" "<<j<<endl;

				}
			}

			i++;

		}
		N++;		
		i = i - 1;
		JJ = J1 +J2;
		//cout<<h<<endl;
		if(j==2) {
			//cout<<"j "<<J1.E()/1000<<endl; 
			h1->Fill(J1.E()/1000);
		}
  	}

	h1->Draw("HIST");

}


