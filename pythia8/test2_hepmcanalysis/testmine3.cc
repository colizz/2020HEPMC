#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
//http://twiki.ific.uv.es/twiki/pub/Atlas/RelRecipeForAraTopQuarkAnalysis/JetPartonSelectorTool.cxx
using namespace Pythia8;


std::vector<const HepMC::GenParticle*> getMothers(const HepMC::GenParticle* child)
{
    std::vector<const HepMC::GenParticle*> particleList;

    if ( !(child->production_vertex()) && (child->status() != 1) )   {//initial state
	return particleList;
    }
    else {                
	if ( !(child->production_vertex()) )   { //something is wrong in HepMC
	    return particleList;
	}
    }
    HepMC::GenVertex::particle_iterator firstMother, thisMother, lastMother;
    const HepMC::GenParticle* mother;   
    firstMother = child->production_vertex()->particles_begin(HepMC::parents);
    lastMother = child->production_vertex()->particles_end(HepMC::parents);
    //Main loop:  iterate over all mothers, call method recursively.
    for( thisMother = firstMother; thisMother != lastMother++; ++thisMother)   {
	mother = (*thisMother);
	particleList.push_back(mother);
    }
    return particleList;
}

const HepMC::GenParticle* getfirstcopy(const HepMC::GenParticle* child)
{
    if ( !(child->production_vertex()) && (child->status() != 1) )   {//initial state
	return child;
    }
    else {                
	if ( !(child->production_vertex()) )   { //something is wrong in HepMC
	    return child;
	}
    }
    HepMC::GenVertex::particle_iterator firstMother, thisMother, lastMother;
    firstMother = child->production_vertex()->particles_begin(HepMC::parents);
    lastMother = child->production_vertex()->particles_end(HepMC::parents);
    if((*firstMother)->pdg_id()== (child)->pdg_id()) {
       return getfirstcopy(*firstMother);
    }
    else {
       return child;
    }
}

std::vector<const HepMC::GenParticle*> getDaughters(const HepMC::GenParticle* mother)
{
    std::vector<const HepMC::GenParticle*> particleList;
//    std::vector<int> particleList;
    
    if ( !(mother->end_vertex()) && (mother->status() != 2) )  {//final state
	return particleList;
    }
    else {                
	if ( !(mother->end_vertex()) )   { //something is wrong in HepMC
	    return particleList;
	};
    };
    HepMC::GenVertex::particle_iterator firstChild, thisChild, lastChild;
    const HepMC::GenParticle* child;        

    firstChild = mother->end_vertex()->particles_begin(HepMC::children);
    lastChild = mother->end_vertex()->particles_end(HepMC::children);
    //Main loop: iterate over all children, call method recursively
    for( thisChild = firstChild; thisChild != lastChild++; ++thisChild)   {
	child = (*thisChild);
	particleList.push_back(child);
    }
    return particleList;
}

const HepMC::GenParticle* getlastcopy(const HepMC::GenParticle* mother)
{
    if ( !(mother->end_vertex()) && (mother->status() != 2) )  {//final state
	return mother;
    }
    else {                
	if ( !(mother->end_vertex()) )   { //something is wrong in HepMC
	    return mother;
	};
    };
    HepMC::GenVertex::particle_iterator firstChild, thisChild, lastChild;
    firstChild = mother->end_vertex()->particles_begin(HepMC::children);
    lastChild = mother->end_vertex()->particles_end(HepMC::children);
    if((*firstChild)->pdg_id()== (mother)->pdg_id()) {
       return getlastcopy(*firstChild);
    }
    else {
       return mother;
    }
}

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
        while ( evt ) {
  		EveN++;
                cout<<"------------"<<endl;
		          for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p!=evt->particles_end(); ++p ){
				if((*p)->pdg_id()==23){

                                        const HepMC::GenParticle* mother1 = getfirstcopy(*p);
                                        cout<<"mother1 "<<mother1->pdg_id()<<" "<<mother1->status()<<endl;
                                        std::vector<const HepMC::GenParticle*> mothers = getMothers(mother1);
                                        cout<<"mothers "<<mothers.size()<<" "<<mothers[0]->pdg_id()<<" "<<mothers[1]->pdg_id()<<endl;


                                        const HepMC::GenParticle* daughter1 = getlastcopy(*p);
                                        cout<<"daughter1 "<<daughter1->pdg_id()<<" "<<daughter1->status()<<endl;
                                        std::vector<const HepMC::GenParticle*> daughters = getDaughters(daughter1);
                                        cout<<"daughters "<<daughters.size()<<" "<<daughters[0]->pdg_id()<<" "<<daughters[1]->pdg_id()<<endl;    
                                         

						ptz      =     daughter1->momentum().perp();   
						ez	=      daughter1->momentum().e();
    						h->Fill();               
                                                continue;
		                 }
                      }
  		

  	delete evt;
  	ascii_in >> evt;
	}

  	h->Write();
  	//delete outFile;
  	return 0;
}

