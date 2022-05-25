#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TGraphErrors.h>

using namespace std;


// ========================== 
// ******** The MAIN ********
// ==========================
int main (int argc, char *argv[]) {
  if ( argc < 5 ) {
    cout << endl
      << "Usage: " << argv[0] << " <list_stations>  <PMT>  <Month> <files>" << endl
      << "  <list_stations>: file with a list of stations" << endl
      << "  <PMT>: ID for the PMT you want to analyse" << endl 
      << "  <files>: IoSd or IoAuger files to be read" << endl
			<< " " << endl
			<< "In case you want the distribution of all events for a specific Station, " << endl
			<< "just make sure the stationsFile conteins a single station." << endl
			<< endl;
    exit(0);
  }	

  const char* stationsFileName = argv[1];
  const char* whichpmt = argv[2];
  AugerIoSd input(argc-3, argv+3);
  const unsigned int totalNrEvents = input.NumberOfEvents();
  ifstream stationsFile(stationsFileName, ios::in);

  if (!stationsFile.is_open()) {
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }

  vector<unsigned int> stationsIds;
  vector< unsigned int > utc2read;
  while (stationsFile.good()) {
    unsigned int st = 0;
    double gps = 0.;
    double utc = 0.;
    stationsFile >> st >> gps >> utc;
    if (st) {
      stationsIds.push_back(st);
      utc2read.push_back(utc);
    }
  }
  
  if (stationsIds.empty()) {
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
  
  TString nameStati = to_string( stationsIds[0] );
  TString pmtname = whichpmt;
  int pmtId= atoi( pmtname );
  if ( pmtId > 0 && pmtId < 4 )
     pmtname = "PMT"+to_string( pmtId );
  else {
    cout << "==================================================" << endl;
    cout << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
    cout << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " << endl;
    cout << "==================================================" << endl;
    exit(0);
	}
  
  cerr << "You have selected PMT: " << pmtname << endl;

  unsigned int previusEvent = 0; // Avoiding read the same event
  unsigned int nrEventsRead = 0;
  unsigned int nrEvents = 0;
  bool found = false;
  //int st2PosVec = 0;
  fstream outTraceLg; 
  fstream outTraceHg;

  pmtname = "results/traceLgSt"+to_string( stationsIds[0] )
    +"Pmt"+to_string( pmtId )+".dat";
  outTraceLg.open(pmtname, ios_base::out);
  pmtname = "results/traceHgSt"+to_string( stationsIds[0] )
    +"Pmt"+to_string( pmtId )+".dat";
  outTraceHg.open(pmtname, ios_base::out);

  EventPos pos;
  // Moving through events
  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    nrEventsRead++;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }
    IoSdEvent event(pos);
    if ( event.Id == previusEvent )
      continue;
    previusEvent = event.Id;
    // Moving through stations inside events
    for (unsigned int i = 0 ; i < event.Stations.size(); ++i) {
      found = false;
      // Searching for desire station        
      for ( unsigned int st_i=0; st_i<stationsIds.size(); st_i++ ) {        
        // Two seconds forward respect GPS time, who know's why
        if ( event.utctime() == utc2read[st_i]+2 )
          if (event.Stations[i].Id == stationsIds[st_i] )
            found = true;
      }
      if ( !found )
        continue;
      
      if ( event.Stations[i].IsUUB ) {
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1 << endl;
        if (event.Stations[i].Error==256) {
          outTraceLg << event.Id << " " << event.utctime();
          outTraceHg << event.Id << " " << event.utctime();
          for ( unsigned int k=0; k<event.Stations[i].UFadc->NSample; k++ ) {
            outTraceLg << " ";
            outTraceHg << " ";
            outTraceLg << event.Stations[i].UFadc->GetValue(pmtId-1,1,k);
            outTraceHg << event.Stations[i].UFadc->GetValue(pmtId-1,0,k);
          }
          outTraceLg << endl;
          outTraceHg << endl;
        }
      }
    }
  }
  outTraceLg.close();
  outTraceHg.close();  
 
	return 0;
}
