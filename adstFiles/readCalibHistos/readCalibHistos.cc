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

#include <iomanip> 
#include <fstream>

using namespace std;


int main ( int argc, char** argv) {
  if ( argc < 2 ) { cout << endl
      << "Usage: " << argv[0] << " <st_time_list> <adst_file>" << endl
      << "<st_time_list>: ASCII with St. IDs and GPS time to read" << endl
      << "<adst_file>: ADST File to read"
      << endl;
    exit(0);
  }
  cout << endl << endl;
  // Reading the ASCII file
  const char* stationsFileName = argv[1];
  ifstream stationsFile(stationsFileName, ios::in);
  if (!stationsFile.is_open()){
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  vector< unsigned int > st2read;
  vector< int > gps2read;
  unsigned int st = 0;
  double gps = 0.;
  double utc = 0.;
  while (stationsFile.good()) {
     st = 0;
     stationsFile >> st >> gps >> utc;
     if (st) {
       st2read.push_back(st);
       gps2read.push_back(gps);
     }
  }
  // Creating file to store charge histograms
  fstream file;
  string outFileHistName = "nan";
  bool readSt = false;
  int nEvtFound = 0;
  file.open("events.dat", ios_base::out);
  bool stFound = false;

  // Reading ADST files
  for ( int adst_i = 2; adst_i<argc; adst_i++ ) {
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
      stFound = false;
      // Reading station in the event
      for ( int st_i=0; st_i<nCand; st_i++ ) {
        readSt = false;
        //cout << sdEvent.GetStation(st_i)->GetId() << endl;
        if ( sdEvent.GetStation(st_i)->GetId() == 1205 )
          stFound = true;
        // Looking if there is a desired station
        for ( int st2r_i=0; st2r_i<st2read.size(); st2r_i++ ) {
          //cout << sdEvent.GetGPSSecond() << " " << (int)gps2read[st2r_i] << endl;
          if ( sdEvent.GetGPSSecond() == gps2read[st2r_i] ) {
            if ( sdEvent.GetStation(st_i)->GetId() == st2read[st2r_i] ) {
              readSt = true;
              file << st2read[st2r_i] << " "
                << sdEvent.GetEventId() << " "
                << sdEvent.GetGPSSecond() << endl;
            }
          }
        }
        if ( !readSt )
          continue;
        nEvtFound++;
        //cout << sdEvent.GetEventId() << endl;
        /*
        // Reading and storing Charge histogram
        if ( sdEvent.HasCalibrationHistograms() ) {
          // Fetching events traces
          const vector<Traces>& traces = 
            sdEvent.GetStation(st_i)->GetPMTTraces();
          for (vector<Traces>::const_iterator trIter = traces.begin();
              trIter != traces.end(); ++trIter) {
            if ( trIter->GetType() != eTotalTrace )
              continue; 
            // Excluding SSD-PMT and SPMT
            if ( trIter->GetPMTId() == 4 || trIter->GetPMTId() == 5 )
              continue;
            const CalibHistogram& calibHisto = trIter->GetChargeHistogram();
            vector <double> bins = calibHisto.GetBinning();
            vector <int> cnts = calibHisto.GetValues();
            outFileHistName = Form("hst_St%d_PMT%d_gps%d.dat",
                sdEvent.GetStation(st_i)->GetId(), trIter->GetPMTId(), 
                sdEvent.GetGPSSecond());
            file.open(outFileHistName, ios_base::out);
            for ( int bin_i=0; bin_i<bins.size(); bin_i++ )
              file << bins[bin_i] << " " << cnts[bin_i] << endl;
            file.close();
          }
        }
        */
      }
      cout << "MSD " << sdEvent.GetEventId() << " " << stFound << endl;
    }
  }
  cout << nEvtFound << endl;
  file.close();
  return 0;
}
