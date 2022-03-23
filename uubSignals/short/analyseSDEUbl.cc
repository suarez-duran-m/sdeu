#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TText.h>
#include <TLine.h>

#include "fitcharge.h"

using namespace std;

double getmean( vector<double> arr ){
  double mean = 0.;
  for (auto & i : arr)
    mean += i;
  return mean/arr.size();
}

float getrms( int arr[], float meanarr ){
  float rms = 0.;
  for (int i=0; i<100; i++)
    rms += (arr[i] - meanarr)*(arr[i] - meanarr);
  return rms/100.;
}

int main (int argc, char *argv[]) {
   if ( argc < 4 ) {
     cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <output>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <output>: output file with whatever you want inside" << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl;
    exit(0);
  }
  // Open files with raw data
  const char* stationsFileName = argv[1];
  AugerIoSd input(argc-3, argv+3);
  const unsigned int totalNrEvents = input.NumberOfEvents();  
  ifstream stationsFile(stationsFileName, ios::in);
  if (!stationsFile.is_open()){
    cout << "Could not open file: " << stationsFileName << endl;
    exit(0);
  }
  // Reading stations to study from file
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
 
  // For Qpk calculation
  fitcharge fittingQpk;
  TH1F *receivedChHisto = new TH1F();
  TH1F *receivedChCrr = new TH1F(); 
  vector < vector < vector < double > > > qpk(stationsIds.size());
  vector < vector < vector < double > > > qpkErr(stationsIds.size());
  vector < vector < vector < double > > > qpkTime(stationsIds.size());
  // Re-size for PMTs
  for ( unsigned int st_i=0; st_i<stationsIds.size(); st_i++ ) {
    qpk[st_i].resize(3);
    qpkErr[st_i].resize(3);
    qpkTime[st_i].resize(3);
  }

  unsigned int previusEvent = 0;
  int sameUtc = 0;
  int stPosVect = 0;

  unsigned int nrEventsRead = 0;
  EventPos pos;

  // Reading Events
  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0)
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;    
    bool found = false;
    IoSdEvent event(pos);
    // Searching is the station triggered in this event
    for (unsigned int st_i = 0 ; st_i < event.Stations.size(); ++st_i) {
      found = false;
      for ( unsigned int stList_i=0; stList_i<stationsIds.size(); stList_i++ )
        if (event.Stations[st_i].Id == stationsIds[stList_i]) {
          found = true;
          stPosVect = stList_i;
        } 
      if (!found)
        continue;
      // Asking if it is UUB
      if ( event.Stations[st_i].IsUUB && event.Id != previusEvent 
          && sameUtc != event.utctime() ) {
        previusEvent = event.Id;
        sameUtc = event.utctime();
        cout << "# Event " << event.Id << " Station " << event.Stations[st_i].Id 
          << endl;
        IoSdEvent event(pos);
        cout << "MSD " << event.Stations[st_i].Error << endl;
        // Filter for error 
        if ( !(event.Stations[st_i].Error==256) )
          continue;
        // Reading by PMT
        for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
          receivedChHisto = event.Stations[st_i].HCharge(pmt_i);
          fittingQpk.setChCrr(*receivedChHisto, 
              event.Stations[st_i].Histo->Offset[pmt_i+6], 
              Form("%d-%d-%ld", stPosVect, pmt_i, event.utctime()));
          receivedChCrr = fittingQpk.getChCrr();
          fittingQpk.getFitCh(*receivedChCrr, 30 );
          if ( fittingQpk.vemPosCh > 0. ) {
            qpk[stPosVect][pmt_i].push_back( fittingQpk.vemPosCh );
            qpkErr[stPosVect][pmt_i].push_back( fittingQpk.vemPosChErr );
            qpkTime[stPosVect][pmt_i].push_back( event.utctime() );
          }
        }
      }
    }
  }

  // Printing Qpk values  
  fstream output;
  output.open("qpks.dat", ios_base::out);

  for ( unsigned int st_i=0; st_i<qpk.size(); st_i++ )
    for ( unsigned int pmt_i=0; pmt_i<3; pmt_i++ )
      for (auto & qpk_i : qpk[st_i][pmt_i] ) 
        output << stationsIds[st_i] << " " << pmt_i+1 << " " << qpk_i << endl;
  output.close();
  
  return 0;
}
