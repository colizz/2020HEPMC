#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"

using namespace Pythia8;

int main() {

	// Create file
	TFile* outFile = new TFile("mg267p8mlm_hepmc.root", "RECREATE");

	// Book Tree.
	double ptz, ez;
	TTree *h = new TTree("h","h");
 	h->Branch("ptz",&ptz);
 	h->Branch("ez",&ez);

  	// specify an input file
	HepMC::IO_GenEvent ascii_in("mg267p8z012mlm.hepmc",std::ios::in);

  	// get the first event
  	HepMC::GenEvent* evt = ascii_in.read_next_event();
	int EveN=0, vbarcode=-1;
        HepMC::GenEvent::vertex_iterator tempt_v;
        HepMC::GenVertex::particles_out_const_iterator tempt_p;
        while ( evt ) {
  		EveN++;
                int isw=0;
		if(EveN%1==0) cout<<"event_number= "<<EveN<<endl;
                cout<<"------------"<<endl;
		for ( HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();v != evt->vertices_end(); ++v ){
                          if(isw!=0) cout<<"--"<<endl;
		          for ( HepMC::GenVertex::particles_out_const_iterator p = (*v)->particles_out_const_begin(); p!=(*v)->particles_out_const_end(); ++p ){
                                isw=0; 
				if((*p)->pdg_id()==23){
					isw=1; 
                                        vbarcode	=	(*p)->end_vertex()->barcode();                             
                                        cout<<" "<<vbarcode<<" "<<(*p)->status()<<endl;
					if((*p)->status()==62)	{
						ptz      =       (*p)->momentum().perp();
						ez	=	(*p)->momentum().e();
    						h->Fill();               
                		
					 }	
		                 }
                           }  
                }
  		

  	delete evt;
  	ascii_in >> evt;
	}

  	h->Write();
  	//delete outFile;
  	return 0;
}

