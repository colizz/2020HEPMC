#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"

using namespace Pythia8;

int main() {

   HepMC::Pythia8ToHepMC ToHepMC;
   HepMC::IO_GenEvent ascii_io("1.hep", std::ios::out);

  Pythia pythia;
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = ttbar.lhe");

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


  pythia.init();

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  for (int iEvent = 0; ; ++iEvent) {
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) break;
      if (++iAbort < nAbort) continue;
      break;
    }

    // Construct new empty HepMC event and fill it.
    // Units will be as chosen for HepMC build, but can be changed
    // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    ToHepMC.fill_next_event( pythia, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

  // End of event loop. Statistics.
  }
  pythia.stat();

  // Done.
  return 0;
}
