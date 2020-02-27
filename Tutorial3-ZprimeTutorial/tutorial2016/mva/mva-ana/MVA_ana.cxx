/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: ClassApplication                                                   *
 *                                                                                *
 * Test suit for comparison of Reader and standalone class outputs                *
 **********************************************************************************/

#include <vector>
#include "MVA_ana.h"


MVA_ana::MVA_ana(){
        for(Int_t i=0; i<3; ++i) {
                ntrzp[i]=-1;
                ntrdy[i]=-1;
                ntrtt[i]=-1;
                ntrww[i]=-1;
                ntrst[i]=-1;
                trcut[i]=-1;
                trsig[i]=-1;
        }
}

TTree* MVA_ana::makeSelectedTree(TTree* origtree,TCut cuts){

cerr<<"The number of events in the original tree "<<origtree->GetName()<<"="<<origtree->GetEntries()<<"\n";
cerr<<"applying selection: "<<cuts.Print()<<"\n";
cerr<<"to tree: "<<origtree->GetName()<<"\n";
TEventList* evlist=new TEventList("cursel");
origtree->Draw(">>cursel",cuts);
//TTree::Draw can be used to fill a TEventList object (list of entry numbers)
//instead of histogramming one variable.
cerr<<"number of entries in the event list: "<<evlist->GetN()<<"\n";
//TFile *tmpoutput = TFile::Open("tmp.root","RECREATE");
TTree* t_clone=origtree->CloneTree(0);
//t_clone->SetDirectory(origtree->GetDirectory());
t_clone->SetDirectory(0);

for(int iev=0;iev<evlist->GetN();iev++){
//cerr<<"iev: "<<iev<<"\t-> ";
origtree->GetEntry(evlist->GetEntry(iev));
//cerr<<" got entry; filling\t";
t_clone->Fill();
//cerr<<" filled\n";
}
delete evlist;

//cerr<<"number of entries in the return tree: "<<t_clone->GetEntries()<<"\n";

return t_clone;
}



void MVA_ana::turnontmva(TTree *tree){
tree->Branch("BDT",&BDT,"BDT/F");
tree->Branch("BDTG",&BDTG,"BDTG/F");
tree->Branch("Fisher",&Fisher,"Fisher/F");
tree->Branch("LikelihoodD",&LikelihoodD,"LikelihoodD/F");
tree->Branch("MLP",&MLP,"MLP/F");
}


void MVA_ana::filltmva(){
BDT=valtmva[0];
BDTG=valtmva[1];
Fisher=valtmva[2];
LikelihoodD=valtmva[3];
MLP=valtmva[4];
}

void MVA_ana::findtraincut(TH1 *ha,TH1 *hb,TH1 *hc,TH1 *hd,TH1 *he,Double_t init_v,Double_t final_v,Int_t imva){
        Double_t cutbb,habb,hbbb,hcbb,hdbb,hebb;
        Double_t sig,maxsig;
        Int_t i;
        Double_t step =  final_v - init_v;
        step = step/float(100);
        maxsig=0.0;

        for(i=0; i<100; ++i ){
        cutbb=init_v + step*float(i);
        habb=ha->Integral(i+1,100,"");
        hbbb=hb->Integral(i+1,100,"");
        hcbb=hc->Integral(i+1,100,"");
        hdbb=hd->Integral(i+1,100,"");
        hebb=he->Integral(i+1,100,"");

        if (i==0){
                        std::cout << "before cuts "<< std::endl;
                        std::cout <<    "ntrzp"<<","<<"ntrdy"<<","<<"ntrtt"<<","<<"ntrww"<<","<<"ntrst"<<endl;
                        std::cout << habb << "    "<<hbbb << "    "<<hcbb << "    "<<hdbb << "    "<<hebb <<std::endl;
                        }
               //sig=(habb)/TMath::Sqrt(hbbb+hcbb+hdbb+hebb+0.0001);
                sig=TMath::Sqrt(2*((habb+hbbb+hcbb+hdbb+hebb)*TMath::Log(1+(habb/(hbbb+hcbb+hdbb+hebb)))-habb)); 

                if(sig > maxsig) {
                maxsig=sig;
                ntrzp[imva]=habb;
                ntrdy[imva]=hbbb;
                ntrtt[imva]=hcbb;
                ntrww[imva]=hdbb;
                ntrst[imva]=hebb;
                trcut[imva]=cutbb;
                trsig[imva]=maxsig;
                }}

}
 

void MVA_ana::Classification(TCut cuts, TString myMethodList, TString datafilename){
// this loads the library
   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;
//
//   Use["Cuts"]            = 1;
   Use["BDT"]             = 1;
   Use["BDTG"]            = 1;
   Use["Fisher"]          = 1;
   Use["LikelihoodD"]     = 1; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["MLP"]             = 1; // this is the recommended ANN

   // ---------------------------------------------------------------
   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

//   TString myMethodList=argv[1];

   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }
// Create a new root output file.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar" );
   factory->AddVariable("mll", "mll","GeV", 'D');
   factory->AddVariable("ptll", "ptll","GeV", 'D');
   factory->AddVariable("ptMu1", "ptMu1","GeV", 'D');
   factory->AddVariable("ptMu2", "ptMu2","GeV", 'D');
   factory->AddVariable("detall", "detall","", 'D');
   factory->AddVariable("dphill", "dphill","", 'D');

   // read training and test data
   // load the signal and background event samples from ROOT trees
   TFile *input(0);

   TString IFile =  "../../macro/" + datafilename;
   input = TFile::Open( IFile ); // check if file in local directory exists

   if (!input) {
         std::cout << "ERROR: could not open data file" << std::endl;
         exit(1);
      }
   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

        TTree *treea  = (TTree*)input->Get("ZP2000");
        TTree *treeb  = (TTree*)input->Get("DY");
        TTree *treec  = (TTree*)input->Get("ST");
        TTree *treed  = (TTree*)input->Get("TT");
        TTree *treee  = (TTree*)input->Get("WW");
   // global event weights per tree (see below for setting event-wise weights)

      factory->AddSignalTree    ( treea,    15.58/10000.);
      factory->AddBackgroundTree( treeb,     7.48/80000.);
      factory->AddBackgroundTree( treec,     0.089/10000.);
      factory->AddBackgroundTree( treed,     0.76/40000.);
      factory->AddBackgroundTree( treee,     0.338/40000.);
   /*   factory->AddSignalTree    ( treea,     weight);
      factory->AddBackgroundTree( treeb,     weight);
      factory->AddBackgroundTree( treec,     weight);
      factory->AddBackgroundTree( treed,     weight);
      factory->AddBackgroundTree( treee,     weight);
      factory->AddBackgroundTree( treef,     weight);
  
*/
    factory->PrepareTrainingAndTestTree( cuts, cuts,
                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=EqualNumEvents:!V" );

  // BDT
   /*if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=400:nEventsMin=1000:MaxDepth=10:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
*/
   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=100:nEventsMin=3000:MaxDepth=8:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20" );  
// Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=100:BoostType=Grad:Shrinkage=0.3:UseBaggedGrad:GradBaggingFraction=0.6:SeparationType=GiniIndex:nCuts=20:NNodesMax=6" );
  // Fisher discriminant   
   if (Use["Fisher"])
      factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=60:NsmoothMVAPdf=10" );
   if (Use["LikelihoodD"])
      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", 
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" ); 
   if (Use["MLP"])
      factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=500:HiddenLayers=N+5:TestRate=10:EpochMonitoring" );

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethodsForClassification();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   input->Close();

}

void MVA_ana::turnon(TTree *tree){
tree->SetBranchAddress("mll", &mll);
tree->SetBranchAddress("ptll",&ptll);
tree->SetBranchAddress("ptMu1",&ptMu1);
tree->SetBranchAddress("ptMu2",&ptMu2);
tree->SetBranchAddress("detall",&detall);
tree->SetBranchAddress("dphill",&dphill);
}

void MVA_ana::fill(std::vector<double>* input){
//(*input)[0]=weight;
(*input)[0]=mll;
(*input)[1]=ptll;
(*input)[2]=ptMu1;
(*input)[3]=ptMu2;
(*input)[4]=detall;
(*input)[5]=dphill;
}


void MVA_ana::ClassApplication(TCut cuts, TString myMethodList, TString datafilename, TString rawfile, TString tmvafile){

  TFile *seltarget  = new TFile( rawfile,"RECREATE" );
  TFile *tmvtarget  = new TFile( tmvafile,"RECREATE" );

   
   cout << "==> start ClassApplication" << endl;

   const int Nmvas = 5;

   const char* bulkname[Nmvas] = {"BDT","BDTG","Fisher","LikelihoodD","MLP"}; 
 //const char* bulkname[Nmvas] = {"BDT"};

   bool iuse[Nmvas] = { Nmvas*kFALSE };
   //cout << bulkname[1] <<endl;
   // interpret input list
   if (myMethodList != "") {
      TList* mlist = TMVA::gTools().ParseFormatLine( myMethodList, " :," );
      for (int imva=0; imva<Nmvas; imva++) if (mlist->FindObject( bulkname[imva] )) iuse[imva] = kTRUE;
      delete mlist;
   }
  //cout << iuse[1] << endl;
   // create a set of variables and declare them to the reader
   // - the variable names must corresponds in name and type to 
   // those given in the weight file(s) that you use
   std::vector<std::string> inputVars;
   std::vector<double>* inputVec = new std::vector<double>(15);
//   inputVars.push_back("weight");
   inputVars.push_back("mll");
   inputVars.push_back("ptll");
   inputVars.push_back("ptMu1");
   inputVars.push_back("ptMu2");
   inputVars.push_back("detall");
   inputVars.push_back("dphill");

   // preload standalonetcostp class(es)
   string dir    = "weights/";
   string prefix = "TMVAClassification";

   for (int imva=0; imva<Nmvas; imva++) {
      if (iuse[imva]) {
         TString cfile = dir + prefix + "_" + bulkname[imva] + ".class.C++";

         cout << "=== Macro        : Loading class  file: " << cfile << endl;

         // load the classifier's standalone class
         gROOT->LoadMacro( cfile );
      }
   }
   cout << "=== Macro        : Classifier class loading successfully terminated" << endl;

   // define classes
   IClassifierReader* classReader[Nmvas] = { Nmvas*0 };

   // ... and create them (and some histograms for the output)
   int nbin = 100;
   TH1 *hista[Nmvas],*histb[Nmvas],*histc[Nmvas],*histd[Nmvas],*histe[Nmvas];

   for (int imva=0; imva<Nmvas; imva++) {
      if (iuse[imva]) {
         cout << "=== Macro        : Testing " << bulkname[imva] << endl;
         if (bulkname[imva] == "BDT"          ) {
            classReader[imva] = new ReadBDT          ( inputVars );
            hista[imva] = new TH1F( "MVA_BDT_zp",           "MVA_BDT",           nbin, -1, 1 );
            histb[imva] = new TH1F( "MVA_BDT_dy",           "MVA_BDT",           nbin, -1, 1 );
            histc[imva] = new TH1F( "MVA_BDT_tt",           "MVA_BDT",           nbin, -1, 1 );
            histd[imva] = new TH1F( "MVA_BDT_ww",           "MVA_BDT",           nbin, -1, 1 );
            histe[imva] = new TH1F( "MVA_BDT_st",           "MVA_BDT",           nbin, -1, 1 );
         }
         if (bulkname[imva] == "BDTG"          ) {
            classReader[imva] = new ReadBDTG          ( inputVars );
            hista[imva] = new TH1F( "MVA_BDTG_zp",           "MVA_BDTG",         nbin, -1, 1 );
            histb[imva] = new TH1F( "MVA_BDTG_dy",           "MVA_BDTG",         nbin, -1, 1 );
            histc[imva] = new TH1F( "MVA_BDTG_tt",           "MVA_BDTG",         nbin, -1, 1 );
            histd[imva] = new TH1F( "MVA_BDTG_ww",           "MVA_BDTG",         nbin, -1, 1 );
            histe[imva] = new TH1F( "MVA_BDTG_st",           "MVA_BDTG",         nbin, -1, 1 );
         }
         if (bulkname[imva] == "Fisher"       ) {
            classReader[imva] = new ReadFisher       ( inputVars );
            hista[imva] = new TH1F( "MVA_Fisher_zp",        "MVA_Fisher",        nbin, -4, 4 );
            histb[imva] = new TH1F( "MVA_Fisher_dy",        "MVA_Fisher",        nbin, -4, 4 );
            histc[imva] = new TH1F( "MVA_Fisher_tt",        "MVA_Fisher",        nbin, -4, 4 );
            histd[imva] = new TH1F( "MVA_Fisher_ww",        "MVA_Fisher",        nbin, -4, 4 );
            histe[imva] = new TH1F( "MVA_Fisher_st",        "MVA_Fisher",        nbin, -4, 4 );
         }
        if (bulkname[imva] == "LikelihoodD") {
            classReader[imva] = new ReadLikelihoodD( inputVars );
            hista[imva] = new TH1F( "MVA_LikelihoodD_zp", "MVA_LikelihoodMIX", nbin,  0, 1 );
            histb[imva] = new TH1F( "MVA_LikelihoodD_dy", "MVA_LikelihoodMIX", nbin,  0, 1 );
            histc[imva] = new TH1F( "MVA_LikelihoodD_tt", "MVA_LikelihoodMIX", nbin,  0, 1 );
            histd[imva] = new TH1F( "MVA_LikelihoodD_ww", "MVA_LikelihoodMIX", nbin,  0, 1 );
            histe[imva] = new TH1F( "MVA_LikelihoodD_st", "MVA_LikelihoodMIX", nbin,  0, 1 );
         }
        if (bulkname[imva] == "MLP"          ) {
            classReader[imva] = new ReadMLP          ( inputVars );
            hista[imva] = new TH1F( "MVA_MLP_zp",           "MVA_MLP",           nbin, -1.2, 1.2 );
            histb[imva] = new TH1F( "MVA_MLP_dy",           "MVA_MLP",           nbin, -1.2, 1.2 );
            histc[imva] = new TH1F( "MVA_MLP_tt",           "MVA_MLP",           nbin, -1.2, 1.2 );
            histd[imva] = new TH1F( "MVA_MLP_ww",           "MVA_MLP",           nbin, -1.2, 1.2 );
            histe[imva] = new TH1F( "MVA_MLP_st",           "MVA_MLP",           nbin, -1.2, 1.2 );
         }
      }
   }
   cout << "=== Macro        : Class creation was successful" << endl;
   

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
   TString IFile =  "../../macro/" + datafilename;
   input = TFile::Open( IFile ); 

   if (!input) {
      cout << "ERROR: could not open data file" << endl;
      exit(1);
   }
   
   //
   // prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //

        TTree *treea  = (TTree*)input->Get("ZP2000");
        TTree *treeb  = (TTree*)input->Get("DY");
        TTree *treec  = (TTree*)input->Get("TT");
        TTree *treed  = (TTree*)input->Get("WW");
        TTree *treee  = (TTree*)input->Get("ST");

        TTree *clone;

        TTree *seltreea = new TTree("selzptr", "Selected events from zp");
        TTree *seltreeb = new TTree("seldytr", "Selected events from dy");
        TTree *seltreec = new TTree("seltttr", "Selected events from tt");
        TTree *seltreed = new TTree("selwwtr", "Selected events from ww");
        TTree *seltreee = new TTree("selsttr", "Selected events from st");

        TTree *tmvtreea= new TTree("tmva","the results of tmva");
        TTree *tmvtreeb= new TTree("tmvb","the results of tmva");
        TTree *tmvtreec= new TTree("tmvc","the results of tmva");
        TTree *tmvtreed= new TTree("tmvd","the results of tmva");
        TTree *tmvtreee= new TTree("tmve","the results of tmva");

        MVA_ana::turnontmva(tmvtreea);
        MVA_ana::turnontmva(tmvtreeb);
        MVA_ana::turnontmva(tmvtreec);
        MVA_ana::turnontmva(tmvtreed);
        MVA_ana::turnontmva(tmvtreee);
 

        cout << "=== Macro        : Loop over signal sample" << endl;

        // the references to the variables
        clone=makeSelectedTree(treea,cuts);
        seltreea=clone->CloneTree(0);
        seltreea->SetDirectory(0);
        tmvtreea->SetDirectory(0);
        MVA_ana::turnon(clone);
        for (Long64_t ievt=0; ievt<clone->GetEntries();ievt++) {
                if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;
                clone->GetEntry(ievt);
                seltreea->Fill();
                MVA_ana::fill(inputVec);
                //cout<< "loop over all booked classifiers" <<endl;
                for (int imva=0; imva<Nmvas; imva++) {
                         valtmva[imva]=0.0;
                         if (iuse[imva]) {
                                //cout << weight << endl;
                               // cout<<" retrive the classifier responses" <<endl;           
                                double retval = classReader[imva]->GetMvaValue( *inputVec );
                               //cout<< retval <<endl;
                                hista[imva]->Fill( retval, 15.58/10000.);
                                valtmva[imva]=retval;
                         }
                }
        MVA_ana::filltmva();
        tmvtreea->Fill();
        }

        seltarget->cd();
        seltreea->Write();
        tmvtarget->cd();
        tmvtreea->Write();

   clone=makeSelectedTree(treeb,cuts);
   seltreeb=clone->CloneTree(0);
   seltreeb->SetDirectory(0);
   tmvtreeb->SetDirectory(0);
   MVA_ana::turnon(clone);
   for (Long64_t ievt=0; ievt<clone->GetEntries();ievt++) {
      if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;
      clone->GetEntry(ievt);
      seltreeb->Fill();
          MVA_ana::fill(inputVec);
      // loop over all booked classifiers
      for (int imva=0; imva<Nmvas; imva++) {
         valtmva[imva]=0.0;
         if (iuse[imva]) {
            // retrive the classifier responses            
            double retval = classReader[imva]->GetMvaValue( *inputVec );
            histb[imva]->Fill( retval, 7.48/80000. );
            valtmva[imva]=retval;
         }
      }
      MVA_ana::filltmva();
      tmvtreeb->Fill();
   }

   seltarget->cd();
   seltreeb->Write();
   tmvtarget->cd();
   tmvtreeb->Write();

   clone=makeSelectedTree(treec,cuts);
   seltreec=clone->CloneTree(0);
   seltreec->SetDirectory(0);
   tmvtreec->SetDirectory(0);
   MVA_ana::turnon(clone);
   for (Long64_t ievt=0; ievt<clone->GetEntries();ievt++) {
      if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;
      clone->GetEntry(ievt);
      seltreec->Fill();
          MVA_ana::fill(inputVec);
      // loop over all booked classifiers
      for (int imva=0; imva<Nmvas; imva++) {
         valtmva[imva]=0.0;
         if (iuse[imva]) {
            // retrive the classifier responses            
            double retval = classReader[imva]->GetMvaValue( *inputVec );
            histc[imva]->Fill( retval, 0.089/10000. );
            valtmva[imva]=retval;
         }
      }
      MVA_ana::filltmva();
      tmvtreec->Fill();
   }
//   seltreea->Print();

   seltarget->cd();
   seltreec->Write();
   tmvtarget->cd();
   tmvtreec->Write();


    clone=makeSelectedTree(treed,cuts);
   seltreed=clone->CloneTree(0);
   seltreed->SetDirectory(0);
   tmvtreed->SetDirectory(0);
   MVA_ana::turnon(clone);
   for (Long64_t ievt=0; ievt<clone->GetEntries();ievt++) {
      if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;
      clone->GetEntry(ievt);
      seltreed->Fill();
          MVA_ana::fill(inputVec);
      // loop over all booked classifiers
      for (int imva=0; imva<Nmvas; imva++) {
         valtmva[imva]=0.0;
         if (iuse[imva]) {
            // retrive the classifier responses            
            double retval = classReader[imva]->GetMvaValue( *inputVec );
            histd[imva]->Fill( retval, 0.76/40000. );
            valtmva[imva]=retval;
         }
      }
      MVA_ana::filltmva();
      tmvtreed->Fill();
   }
//   seltreea->Print();

   seltarget->cd();
   seltreed->Write();
   tmvtarget->cd();
   tmvtreed->Write();

    clone=makeSelectedTree(treee,cuts);
   seltreee=clone->CloneTree(0);
   seltreee->SetDirectory(0);
   tmvtreee->SetDirectory(0);
   MVA_ana::turnon(clone);
   for (Long64_t ievt=0; ievt<clone->GetEntries();ievt++) {
      if (ievt%1000 == 0) cout << "=== Macro        : ... processing event: " << ievt << endl;
      clone->GetEntry(ievt);
      seltreee->Fill();
          MVA_ana::fill(inputVec);
      // loop over all booked classifiers
      for (int imva=0; imva<Nmvas; imva++) {
         valtmva[imva]=0.0;
         if (iuse[imva]) {
            // retrive the classifier responses            
            double retval = classReader[imva]->GetMvaValue( *inputVec );
            histe[imva]->Fill( retval,  0.338/40000. );
            valtmva[imva]=retval;
         }
      }
      MVA_ana::filltmva();
      tmvtreee->Fill();
   }
//   seltreea->Print();

   seltarget->cd();
   seltreee->Write();
   tmvtarget->cd();
   tmvtreee->Write();
 
 
   for (int imva=0; imva<Nmvas; imva++) {
      if (iuse[imva]) {
         cout << "=== Macro        : Testing " << bulkname[imva] << endl;
         if (bulkname[imva] == "BDT"          ) {
                        MVA_ana::findtraincut(hista[imva],histb[imva],histc[imva],histd[imva],histe[imva],-1.0,1.0,imva);
         }
         if (bulkname[imva] == "BDTG"          ) {
                        MVA_ana::findtraincut(hista[imva],histb[imva],histc[imva],histd[imva],histe[imva],-1.0,1.0,imva);
         }
         if (bulkname[imva] == "Fisher"       ) {
                        MVA_ana::findtraincut(hista[imva],histb[imva],histc[imva],histd[imva],histe[imva],-4.0,4.0,imva);
         }
         if (bulkname[imva] == "LikelihoodD") {
                        MVA_ana::findtraincut(hista[imva],histb[imva],histc[imva],histd[imva],histe[imva],0.0,1.0,imva);
         }
         if (bulkname[imva] == "MLP"          ) {
                        MVA_ana::findtraincut(hista[imva],histb[imva],histc[imva],histd[imva],histe[imva],-1.2,1.2,imva);
         }
      }
   }


            cout << "=== Macro        : Event loop done! " << endl;

   tmvtarget->cd();

   for (int imva=0; imva<Nmvas; imva++) {
      if (iuse[imva]) {
         hista[imva]->Write();
         histb[imva]->Write();
         histc[imva]->Write();
         histd[imva]->Write();
         histe[imva]->Write();
      }
   }
   cout << "=== Macro        : Created target file: " << seltarget->GetName() << endl;

   seltarget->Close();
   tmvtarget->Close();

   for (int imva=0; imva<Nmvas; imva++) {
      if (iuse[imva]) {
      std::cout <<bulkname[imva]<<std::endl;
          std::cout <<  "ntrzp"<<","<<"ntrdy"<<","<<"ntrtt"<<","<<"ntrww"<<","<<"ntrst"<<","<<"trcut"<<","<<"trsig"<<endl;
          std::cout <<  ntrzp[imva]<<","<<ntrdy[imva]<<","<<ntrtt[imva]<<","<<ntrww[imva]<<","<<ntrst[imva]<<","<<trcut[imva]<<","<<trsig[imva]<<endl;
      }
   }

   delete seltarget;
   delete tmvtarget;
   delete inputVec;
   input->Close();
   cout << "==> ClassApplication is done!" << endl << endl;
}













