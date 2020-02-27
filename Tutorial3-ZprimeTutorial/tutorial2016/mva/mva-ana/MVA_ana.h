#include "Rtypes.h"
#include <iostream>
#include "TMath.h"

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#endif

Bool_t ReadDataFromAsciiIFormat = kFALSE;

class MVA_ana {

public:
        MVA_ana();
        ~MVA_ana(){};

        void findtraincut(TH1 *ha,TH1 *hb,TH1 *hc,TH1 *hd,TH1 *he,Double_t init_v,Double_t final_v,Int_t imva);
        void findtestperformance(TH1 *ha,TH1 *hb,TH1 *hc,TH1 *hd,TH1 *he,Double_t init_v,Double_t final_v,Int_t imva);
        void Classification(TCut cuts, TString myMethodList="", TString datafilename);
        void ClassApplication(TCut cuts, TString myMethodList, TString datafilename, TString selfl, TString tmvfl);
        void Generalization(TCut cuts, TString myMethodList, TString datafilename, TString selfl, TString tmvfl);
        TTree* makeSelectedTree(TTree* origtree,TCut cuts);
        void turnon(TTree *tree);
        void turnontmva(TTree *tree);
        void filltmva();
        void fill(std::vector<double>* input);
        inline void readincuts(Double_t *hcuts){for(Int_t i=0; i<16; ++i){trcut[i]=hcuts[i];}};
        inline void outputcuts(Double_t *hcuts){for(Int_t i=0; i<16; ++i){hcuts[i]=trcut[i];}};

private:
        Double_t ntrzp[16],ntrdy[16],ntrtt[16],ntrww[16],ntrst[16],trcut[16],trsig[16];

        Double_t sf,mll,ptll,ptMu1,ptMu2,detall,dphill;
        
        Float_t Cuts,MLP,MLPBFGS,SVM,CFMlpANN,Fisher,Likelihood,LikelihoodD,LikelihoodPCA,LD,HMatrix,FDA_MT,FDA_MC,FDA_GA,BDT,BDTD,BDTG,valtmva[16];

};

