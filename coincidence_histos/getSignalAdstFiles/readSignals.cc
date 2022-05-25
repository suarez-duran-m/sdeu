#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <RecEvent.h>
#include <RecEventFile.h>
#include <DetectorGeometry.h>
#include <Traces.h>
#include <TraceType.h>

#include <TFile.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

using namespace std;

/********************/
/* Global variables */
/********************/

// Loaded SPMT starting
const int loadStart = 1324857618;

int main ( int argc, char** argv) {
  if ( argc < 2 ) { cout << endl
    << "Usage: " << argv[0] << " <list_stations> <PMT> <adstfiles>" << endl
      << "  <list_stations>: file with a list of stations" << endl
      << "  <PMT>: ID for the PMT you want to analyse" << endl
      << "  <adstfiles>: File with list of adst files to read" << endl;
    exit(0);
  }
  cout << endl << endl;
  
  const char* stationsFileName = argv[1];
  const char* whichpmt = argv[2];
  ifstream stationsFile(stationsFileName, ios::in);

  vector<unsigned int> stationsIds;
  vector< unsigned int > gps2read;
  while (stationsFile.good()) {
    unsigned int st = 0;
    double gps = 0.;
    double utc = 0.;
    stationsFile >> st >> gps >> utc;
    if (st) {
      stationsIds.push_back(st);
      gps2read.push_back(gps);
    }
  }

  if (stationsIds.empty()) {
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }

  int st_posVec = 0;
  bool readSt = false;

  // Reading ADST files
  for ( int adst_i = 3; adst_i<argc; adst_i++ ) {
    RecEventFile inputFile(argv[adst_i]);
    RecEvent *theRecEvent = new RecEvent();
    inputFile.SetBuffers(&theRecEvent);
    cout << "Opened " << argv[adst_i]
      << " with " << inputFile.GetNEvents()
      << " events " << endl;

    // Reading events 
    while ( inputFile.ReadNextEvent() == RecEventFile::eSuccess ) {
      const auto &sdEvent = theRecEvent->GetSDEvent();
      int nCand = sdEvent.GetNumberOfCandidates();
      // Reading for event related stations
      for ( int st_i=0; st_i<nCand; st_i++ ) {
        readSt = false;
        // Looking if there is a desired station
        for ( int st2r_i=0; st2r_i<stationsIds.size(); st2r_i++ ) {
          //if ( sdEvent.GetGPSSecond() == gps2read[st_i] )            
            if ( sdEvent.GetStation(st_i)->GetId() == stationsIds[st2r_i] ) {
              readSt = true;
              st_posVec = st2r_i;
              cout << "MSD " << sdEvent.GetGPSSecond() << " " << gps2read[st_i] << endl;
            }          
        }
        if ( !readSt )
          continue;
        
        // Printing signals
        /*
        cout << sdEvent.GetGPSSecond() << " "
          << sdEvent.GetStation(st_i)->GetTotalSignal() << " " 
          << sdEvent.GetStation(st_i)->GetTotalSignalError() << " " 
          << endl;
          */
      }
    }
  }
  return 0;
}

