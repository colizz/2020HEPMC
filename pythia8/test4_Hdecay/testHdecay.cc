#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"
//http://twiki.ific.uv.es/twiki/pub/Atlas/RelRecipeForAraTopQuarkAnalysis/JetPartonSelectorTool.cxx
//http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
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
	TFile* outFile = new TFile("h.root", "RECREATE");

	// Book Tree.
        double ntot=0.;
	double hbb=0., hww=0., hzz=0., hgg=0., htautau=0., haa=0., hcc=0., hmm=0., hza=0.;
	TTree *h = new TTree("h","h");
 	h->Branch("hbb",&hbb);
 	h->Branch("hww",&hww);
 	h->Branch("hzz",&hzz);
 	h->Branch("hgg",&hgg);
 	h->Branch("htautau",&htautau);
 	h->Branch("haa",&haa);
 	h->Branch("hcc",&hcc);
 	h->Branch("hmm",&hmm);
 	h->Branch("hza",&hza);
 	h->Branch("ntot",&ntot);
  	// specify an input file
	HepMC::IO_GenEvent ascii_in("hdecay.hepmc",std::ios::in);

  	// get the first event
  	HepMC::GenEvent* evt = ascii_in.read_next_event();
	int EveN=0, vbarcode=-1;
        while ( evt ) {
  		EveN++;

		          for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin(); p!=evt->particles_end(); ++p ){
				if((*p)->pdg_id()==25){
                                        const HepMC::GenParticle* daughter1 = getlastcopy(*p);
                                        std::vector<const HepMC::GenParticle*> daughters = getDaughters(daughter1);
                                        //cout<<"daughters "<<daughters.size()<<" "<<daughters[0]->pdg_id()<<" "<<daughters[1]->pdg_id()<<endl;   
                                        if(abs(daughters[0]->pdg_id())==5 & abs(daughters[1]->pdg_id())==5) {
                                         hbb++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==24 & abs(daughters[1]->pdg_id())==24) {
                                         hww++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==23 & abs(daughters[1]->pdg_id())==23) {
                                         hzz++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==21 & abs(daughters[1]->pdg_id())==21) {
                                         hgg++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==15 & abs(daughters[1]->pdg_id())==15) {
                                         htautau++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==22 & abs(daughters[1]->pdg_id())==22) {
                                         haa++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==4 & abs(daughters[1]->pdg_id())==4) {
                                         hcc++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==13 & abs(daughters[1]->pdg_id())==13) {
                                         hmm++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==23 & abs(daughters[1]->pdg_id())==22) {
                                         hza++;
                                        }
                                        else if(abs(daughters[0]->pdg_id())==22 & abs(daughters[1]->pdg_id())==23) {
                                         hza++;
                                        }
                                        else{
                                        cout<<"daughters "<<daughters.size()<<" "<<daughters[0]->pdg_id()<<" "<<daughters[1]->pdg_id()<<endl;   
                                        }
                                        break;
		                    }
                           }
  	delete evt;
  	ascii_in >> evt;
	}
        ntot=double(EveN);  
        cout<<"hbb="<<hbb<<" "<<" hww="<<hww<<" hzz="<<hzz<<" hgg="<<hgg<<" htautau="<<htautau<<" haa="<<haa<<" hcc="<<hcc<<" hmm="<<hmm<<" hza="<<hza<<endl;
        cout<<"ntot="<<ntot<<" "<<hbb+hww+hzz+hgg+htautau+haa+hcc+hmm+hza<<endl;		
        cout<<"hbb="<<hbb/ntot<<" "<<" hww="<<hww/ntot<<" hzz="<<hzz/ntot<<" hgg="<<hgg/ntot<<" htautau="<<htautau/ntot<<" haa="<<haa/ntot<<" hcc="<<hcc/ntot<<" hcc="<<hcc/ntot<<" hza="<<hza/ntot<<endl;
        h->Fill();   
  	h->Write();
  	//delete outFile;
  	return 0;
}

