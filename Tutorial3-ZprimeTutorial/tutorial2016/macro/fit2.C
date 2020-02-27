//root -b -q fit2.C\(\"ZP2000\"\) 

//RooVoigtian is an efficient implementation of the convolution of a Breit-Wigner with a Gaussian, 
//https://root.cern.ch/root/html/RooVoigtian.html
//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPreExerciseFourthSet#Fitting_a_Convolution_of_Gaussia

void fit2(char* dataFileName)
{
	gSystem->Load("libRooFit");
	using namespace RooFit;
    

        TString buffer = dataFileName;
        cout<<buffer<<endl;
        TString file =  "../ntuple/ntuple_"+buffer + ".root";
        TString rlt_file = buffer+".pdf";

	TFile * fp = new TFile(file);
	TTree * tree = (TTree *)fp->Get("tree");


        double inibin=2200.;
        double finbin=3500.;

 

 
//-------Method 1-------------------------
//	RooRealVar x("mll","mll", inibin,finbin);
//	RooDataSet data("data","Dataset With X",tree, x);
//-------Method 2-------------------------
	double invmll;
	tree->SetBranchAddress("mll",&invmll);
	RooRealVar x("x","x", inibin,finbin);
	RooDataSet data("data","Dataset With X",x);
	for(Int_t j=0; j<tree->GetEntries(); j++){
		tree->GetEntry(j); 
		if(j%1000==0) {cout<<"Processing the events: "<<j<<endl;}
		x=invmll;
                if(invmll >inibin && invmll < finbin) { data.add(x);}
	}
//---------------------------------------------

    RooRealVar mean("mean1","mean of gaussian 1",3000,2400,3500) ;
    RooRealVar width("width","width",480.0, 320.0, 600.0);
    RooRealVar sigma("sigma","sigma",200.0, 100.0, 300.0);
    RooVoigtian sig("sig","sig",x,mean,width,sigma);
 


    RooRealVar nsig("nsig","nsig",2500,0,10000);
    RooExtendPdf model("esig","esig",sig,nsig);
    model.fitTo(data,Extended(kTRUE));



 
 
        TCanvas *cfit=new TCanvas("cfit","cfit");
    //    cfit->SetLogy();
	RooPlot* frame = x.frame();
	data.plotOn(frame,Binning(40));
    	model.plotOn(frame);
    //    data.statOn(frame);
 
        frame->Draw();
        cfit->Print(rlt_file);

        cout<<"mean="<<mean<<endl;
        cout<<"width="<<width<<endl;
        cout<<"sigma="<<sigma<<endl;
        cout<<"nsig="<<nsig<<endl;
        cout<<" Fit chi square/dof = "<<frame->chiSquare(3)<<endl; 
}



 
