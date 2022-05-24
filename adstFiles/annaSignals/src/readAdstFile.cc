#include "readAdstFile.h"


readAdstFile::readAdstFile() {
  theRecEvent = new RecEvent();

  sglPmt.resize(3);
  muonTauPmt.resize(3);
  for(int i=0; i<3; i++) {
    sglPmt[i].clear();
    muonTauPmt[i].clear();
  }

  totSgl.clear();
  totSglErr.clear();
  sglGpsSec.clear();
}

void readAdstFile::setAdstFileName(const char* file) {
  fileName = file;
}

void readAdstFile::checkNumberEvts() { 
  RecEventFile adstFile(fileName); 
  cout << adstFile.GetNEvents() << " events in " 
    << fileName << endl;
}

void readAdstFile::readSignals(int st) {
  // Reading events and extracting signals in VEM 
  // for st
  readEvents(st);
}

void readAdstFile::readEvents(int st) {
  RecEventFile adstFile(fileName);
  adstFile.SetBuffers(&theRecEvent);
  int nCand = 0;
  bool readSt = false;
  vector < double > tmpSgl(3);
  vector < double > tmpTau(3);

  while ( adstFile.ReadNextEvent() == RecEventFile::eSuccess ) {
    const auto &sdEvent = theRecEvent->GetSDEvent();
    nCand = sdEvent.GetNumberOfCandidates();
    // Reading each station involved in this event
    for ( int st_i=0; st_i<nCand; st_i++ ) {
      // Looking if st_i is a desired station
      readSt = false;
      if ( sdEvent.GetStation(st_i)->GetId() == st )
        readSt = true;
      if ( !readSt )
        continue;
      // Reading signals in VEM per PMT
      for(int pmt=1; pmt<4; pmt++) { 
        const Traces& tmpTr = sdEvent.GetStation(st_i)->GetPMTTraces(eTotalTrace,pmt);
        tmpSgl[pmt-1] = tmpTr.GetVEMSignal();
        tmpTau[pmt-1] = tmpTr.GetMuonPulseDecayTime();
      }
      for(int pmt=0; pmt<3; pmt++) {
        sglPmt[pmt].push_back( tmpSgl[pmt] );        
        muonTauPmt[pmt].push_back( tmpTau[pmt] );
      }
      totSgl.push_back( sdEvent.GetStation(st_i)->GetTotalSignal() );
      totSglErr.push_back ( sdEvent.GetStation(st_i)->GetTotalSignalError() );
      sglGpsSec.push_back( sdEvent.GetGPSSecond() );
    }
  }
}

const vector < vector < double > > readAdstFile::getSglPmt() {
  return sglPmt;
}

const vector < vector < double > > readAdstFile::getMuonTauPmt() {
  return muonTauPmt;
}

const vector < double > readAdstFile::getTotSgl() {
  return totSgl;
}

const vector < double > readAdstFile::getTotSglErr() {
  return totSglErr;
}

const vector < double > readAdstFile::getSglGpsSec() {
  return sglGpsSec;
}
