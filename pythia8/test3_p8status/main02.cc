// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
#include "TTree.h"
#include "TFile.h"
using namespace Pythia8;

//main73.cc
int getfirstcopy( Event& event, int i) {
      if (event[event[i].mother1()].id() == event[i].id()) {
          return getfirstcopy(event, event[i].mother1());
      }
      else {
          return i;//event[i].mother1();
      }
}

int getlastcopy( Event& event, int i) {
      if (event[event[i].daughter1()].id() == event[i].id()) {
          return getlastcopy(event, event[i].daughter1());
      }
      else {
          return i;//event[i].daughter1();
      }
}

int main() {

  // Create file
  TFile* outFile = new TFile("p8.root", "RECREATE");
  double ptz, ez;
  int status; 
  TTree *t2 = new TTree("t2","t2");
  t2->Branch("status",&status);

  TTree *t = new TTree("t","t");
  t->Branch("ptz",&ptz);
  t->Branch("ez",&ez);
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;

  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = dylo.lhe");

  pythia.readString("Init:showChangedSettings=on");
  pythia.readString("Next:numberShowInfo=0");
  pythia.readString("Next:numberShowProcess=0");
  pythia.readString("Next:numberShowEvent=0");
  pythia.readString("Main:timesAllowErrors = 100");
  pythia.readString("Tune:preferLHAPDF = 2");
  pythia.readString("Beams:setProductionScalesFromLHEF = off");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("ParticleDecays:allowPhotonRadiation = on");

  pythia.readString("Tune:pp 5");
  pythia.readString("Tune:ee 3");
//http://home.thep.lu.se/~torbjorn/pythia81html/Welcome.html Tunes
  pythia.readString("MultipartonInteractions:pT0Ref=2.1006");
  pythia.readString("MultipartonInteractions:ecmPow=0.21057");
  pythia.readString("MultipartonInteractions:expPow=1.6089");
  pythia.readString("MultipartonInteractions:a1=0.00");
  pythia.readString("ColourReconnection:range=3.31257");
//https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_X/Configuration/Generator/python/Pythia8CommonSettings_cfi.py

  pythia.readString("PartonLevel:MPI = on");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");

/*  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("PhaseSpace:mHatMin = 80.");
  pythia.readString("PhaseSpace:mHatMax = 120.");
  pythia.readString("TimeShower:pTmaxMatch=2");
  pythia.readString("SpaceShower:pTmaxMatch = 2");
*/ 

  pythia.init();

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;
    // Loop over particles in event. Find last Z0 copy. Fill its pT.
    int iZ = 0;
    cout<<"------------"<<endl;
    for (int i = 0; i < pythia.event.size(); ++i)
      if (pythia.event[i].id() == 23) { 
          iZ = i; 

          cout<<"m1 "<<pythia.event[iZ].mother1()<<" "<<pythia.event[pythia.event[iZ].mother1()].id()<<" "<<pythia.event[iZ].mother2()<<" "<<pythia.event[pythia.event[iZ].mother2()].id()<<endl;  
          cout<<"d1 "<<pythia.event[iZ].daughter1()<<" "<<pythia.event[pythia.event[iZ].daughter1()].id()<<" "<<pythia.event[iZ].daughter2()<<" "<<pythia.event[pythia.event[iZ].daughter2()].id()<<endl;  
          status=pythia.event[i].status();
          t2->Fill();
      }

    cout<<"my "<<getfirstcopy(pythia.event, iZ)<<" "<<getlastcopy(pythia.event, iZ)<<endl;

    ptz=pythia.event[iZ].pT();
    ez=pythia.event[iZ].e();
  // End of event loop. Statistics. Histogram. Done.
    t->Fill();
  }
  pythia.stat();
  t->Write();
  t2->Write();
  return 0;
}
