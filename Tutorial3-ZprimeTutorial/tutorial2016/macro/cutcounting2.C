void cutcounting2(){
// cut couting based on Ntuple
   Double_t lumi = 1.;

   TString Sig_name = "ZP2000";
   Double_t zpmass = 2000.;
   TString Bkg_name[4] = {"DY","TT","WW","ST"};
   TString massCut[4] = {"Mllpm300","Mllpm200","Mllpm100","Mllpm50"};
   Double_t cut[4] = {300.,200.,100,50};
   double SF,mll;

   TFile * input1 = new TFile ("zprimetraining.root");
   input1->cd("");


  double sum_sig, sum_bkg;

   for(int i_massCut=0;i_massCut<4;i_massCut++){
        sum_sig=0.0;
        sum_bkg=0.0;
   
           TTree *Outtree1 = (TTree*)input1->Get(Sig_name);
           int nentries1 = (int)Outtree1->GetEntries();
           Outtree1->SetBranchAddress("SF",&SF);
           Outtree1->SetBranchAddress("mll",&mll);
           for (int in=0;in<nentries1;in++){
           Outtree1->GetEntry(in);
              if(fabs(mll-zpmass)<cut[i_massCut]){
                sum_sig+=SF*lumi;
               }
           }

        for(int i=0;i<4;i++){
           TTree *Outtree2 = (TTree*)input1->Get(Bkg_name[i]);
           int nentries2 = (int)Outtree2->GetEntries();
           Outtree2->SetBranchAddress("SF",&SF);
           Outtree2->SetBranchAddress("mll",&mll);
           for (int in=0;in<nentries2;in++){
           Outtree2->GetEntry(in);
              if(fabs(mll-zpmass)<cut[i_massCut]){
                sum_bkg+=SF*lumi;
               }
            }
          }

 
       std::cout<<massCut[i_massCut]<<std::endl;
       std::cout<<sum_sig<<" "<<sum_bkg<<std::endl;
       std::cout<<TMath::Sqrt(2*((sum_sig+sum_bkg)*TMath::Log(1+(sum_sig/sum_bkg))-sum_sig))<<" "<<std::endl;

     }


   }


