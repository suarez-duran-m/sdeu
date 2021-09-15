#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

using namespace std;


int main (int argc, char *argv[]) {
  if ( argc < 4 ) {
    cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <output>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <output>: output file (IoAuger/IoSd depending on input)" << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl;
    exit(0);
  }
  const char* stationsFileName = argv[1];
  const char* outputName = argv[2];
  AugerIoSd input(argc-3, argv+3);
  const unsigned int totalNrEvents = input.NumberOfEvents();
  
  ifstream stationsFile(stationsFileName, ios::in);
  if (!stationsFile.is_open()){
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  vector<unsigned int> stationsIds;
  while (stationsFile.good()) {
    unsigned int st = 0;
    stationsFile >> st;
    if (st)
      stationsIds.push_back(st);
  }
  
  if (stationsIds.empty()){
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
  
  cout << "======> I am  selecting events containing stations: "; 
  for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
         iter!= stationsIds.end(); ++iter)
    cout << *iter << " ";
  cout << " <====== " << endl;
  
  IoSd *sdFile = NULL;
  AugerFile *adFile = NULL;
  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
  EventPos pos;
  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {

    TEcEvent calibEv(pos);

    ++nrEventsRead;
    if (nrEventsRead%1000 == 0){
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
    for (unsigned int i = 0; i < calibEv.fCalibStations.size(); i++){
      const TCalibStation& calibStation = calibEv.fCalibStations[i];

      // take just the UUB stations
      if (!calibStation.IsUUB)
        continue;
              
      if (!calibStation.histo())
        continue;
      
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
             iter!= stationsIds.end(); ++iter){
        if (calibStation.Id == *iter){
          found = true;
          break;
        }
      }
    }
    if (found) {
      ++nrEvents;
      IoSdEvent event(pos);
      if (event.RootClassName == "AugerEvent") {
        if (!adFile)
          adFile = new AugerFile(outputName, AugerFile::eWrite);
        adFile->Write(*(event.RawAugerEvent()), false);
      }     
      else { 
        if (!sdFile)
          sdFile = new IoSd(outputName, "w");
        sdFile->Write(event, "");
      }
    }  
  }
  
  if (sdFile)
    sdFile->Close();
  if (adFile) {
    adFile->Close();
  }
  return 0;
}
