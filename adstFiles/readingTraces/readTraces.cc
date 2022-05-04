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

int main ( int argc, char** argv) {
  if ( argc < 2 ) { 
    cout << endl 
      << "Usage: " << argv[0] << " <file_adstfiles> <output>" << endl
      << "  <list_stations>: file with a list of stations" << endl
      << "  <PMT>: ID for the PMT you want to analyse" << endl
      << "  <file_adstfiles>: File with list of adst files to read" << endl
      << endl;
    exit(0);
  }
  
  const char* stationsFileName = argv[1];
  const char* whichpmt = argv[2];
  ifstream stationsFile(stationsFileName, ios::in);
  int pmtId= atoi( whichpmt );

  vector<unsigned int> st2read;
  vector< unsigned int > gps2read;
  unsigned int st = 0;
  double gps = 0.;
  double utc = 0.;
  while (stationsFile.good()) {
    stationsFile >> st >> gps >> utc;
    if (st) {
      st2read.push_back(st);
      gps2read.push_back(gps);
    }
  }

  // Creating output files
  TString pmtname = whichpmt;
  fstream outTrace;
  pmtname = "PMT"+to_string( pmtId );
  pmtname = "results/trHgSt"+to_string( st2read[0] )
    +"Pmt"+to_string( pmtId )+".dat";
  outTrace.open(pmtname, ios_base::out);

  if (st2read.empty()) {
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
 
  bool readSt = false;

  // Reading ADST files
  for ( int adst_i = 3; adst_i<argc; adst_i++ ) {
    RecEventFile inputFile( argv[adst_i] );
    RecEvent *theRecEvent = new RecEvent();
    inputFile.SetBuffers( &theRecEvent );
    cout << "Opened " << argv[adst_i]
      << " with " << inputFile.GetNEvents()
      << " events " << endl;
    // Reading events    
    while ( inputFile.ReadNextEvent() == RecEventFile::eSuccess ) {
      const auto &sdEvent = theRecEvent->GetSDEvent();
      int nCand = sdEvent.GetNumberOfCandidates();
      // Reading stations in the event
      for ( int st_i=0; st_i<nCand; st_i++ ) {
        readSt = false;
        // Looking for desired station
        for ( int st2r_i=0; st2r_i<st2read.size(); st2r_i++ ) {
          if ( sdEvent.GetGPSSecond() == gps2read[st2r_i] )
            if ( sdEvent.GetStation(st_i)->GetId() == st2read[st2r_i] )
              readSt = true;
        }
        if ( !readSt )
          continue;
        if ( sdEvent.HasCalibrationHistograms() ) {
          // Getting the vector with all traces
          const vector<Traces>& traces = sdEvent.GetStation(st_i)->GetPMTTraces();
          //const vector < float >& tr_i = sdEvent.GetStation(st_i)->GetVEMTrace( pmtId );
          //const vector <short unsigned int >& tr_i = sdEvent.GetStation(st_i)->GetLowGainTrace( pmtId );
          const vector <short unsigned int >& tr_i = sdEvent.GetStation(st_i)->GetHighGainTrace( pmtId );
          outTrace << sdEvent.GetEventId() << " " << sdEvent.GetGPSSecond();
          for ( vector<short unsigned int>::const_iterator i=tr_i.begin();
              i != tr_i.end(); ++i)
            outTrace << " " << *i;
          outTrace << endl;
        }
      }
    }
  } 

  return 0;
}
