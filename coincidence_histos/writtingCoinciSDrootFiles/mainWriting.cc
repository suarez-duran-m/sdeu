#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "include/rawCoinciHistoData.h"
#include <IoAuger.h>

using namespace std;

/********************/
/* Global variables */
/********************/


int main ( int argc, char *argv[]) {
  int minArg = 4;
  if ( argc < minArg ) {
    cout << endl << "=========================" << endl << endl
      << "Usage: " << argv[0] << " sd*packs_file ad_output_name ad_files" << endl;
    cout << endl 
      << "sd*packs_file: file with coincidence histos" << endl
      << "sd_output_name: filename for new ad*.root with coincidence histos" << endl
      << "sd_files: ad*.root files where coincidence histos will be added." << endl 
      << endl;
    exit(0);
  } 

  const char *fileWithCoinci = argv[1];
  const char *outPutSDfileName = argv[2];
  rawCoinciHistoData rawCoinciHistoData;
  rawCoinciHistoData.readData( fileWithCoinci );
  // 
  // Charging coincidence histograms
  vector < int > utcChisto = rawCoinciHistoData.getUtcEvtWidthChisto();
  vector < int > stChisto = rawCoinciHistoData.getStWidthChisto();
  vector < vector < vector < int > > > cQhisto = rawCoinciHistoData.getCQhisto();
  vector < vector < vector < int > > > cHeigth = rawCoinciHistoData.getCheight();
  // Reading ad*.root files
  //
  IoSd sdFileInput( argc-(minArg-1), argv+(minArg-1) );
  const unsigned int totalNrEvents = sdFileInput.NumberOfEvents();  
  IoSd sdOutPutFile(outPutSDfileName, "w");
  IoSdHistoCoinci *newCoinciHistos;
  unsigned int nrEventsRead = 0;
  EventPos pos;
  // Moving through events
  //
  for (pos = sdFileInput.FirstEvent(); pos < sdFileInput.LastEvent(); 
      pos = sdFileInput.NextEvent()) {
    IoSdEvent event(pos);    
    ++nrEventsRead;
    
    if (nrEventsRead%1000 == 0)
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
    //
    // Checking if current evt-UTC match with some CHisto-UTC
    vector < int > indexForUtcMatch;
    bool utcMatched = false;
    for ( unsigned int utc_pos=0; utc_pos<utcChisto.size(); utc_pos++ )
      if ( fabs(utcChisto[utc_pos] - event.utctime()) < 3 ) {
        indexForUtcMatch.push_back( utc_pos );
        utcMatched = true;
      }
    //
    // Ignoring unmatched events
    if ( !utcMatched ) {
      sdOutPutFile.Write(event);
      continue;
    }
    //
    // Moving through stations inside the event        
    for (unsigned int evtSt_i = 0; evtSt_i < event.Stations.size(); ++evtSt_i) {
      if ( !event.Stations[evtSt_i].IsUUB )
        continue;
      //
      // Looking if current st matches with CHisto      
      for ( int i=0; i<indexForUtcMatch.size(); i++ ) {
        int match_i = indexForUtcMatch[i];
        if ( !(event.Stations[evtSt_i].Id == stChisto[match_i]) )
          continue;
        //
        // Storing coincidence histos
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id << endl;
        newCoinciHistos = event.Stations[evtSt_i].SetHistoCoinci(
            cHeigth[match_i], cQhisto[match_i] );
        for ( int j=0; j<10; j++ ) // from NB_HISTO_CALIB
          newCoinciHistos->Offset[j] = event.Stations[evtSt_i].Histo->Offset[j];
        event.Stations[evtSt_i].HasCoinciHisto = 1;
        event.Stations[evtSt_i].HistoCoinci = newCoinciHistos;
        // Histo->Entries it is not the entries into the
        // respective histo
        event.Stations[evtSt_i].HistoCoinci->Entries =
          event.Stations[evtSt_i].Histo->Entries;
        int cnt = 0;
        if ( event.Stations[evtSt_i].Id == 871 && (event.utctime() - 315964782) == 1343001246 ) {
          for( int i=0; i<100; i++ ) //cQhisto[match_i][0][:100] )
            cnt += cQhisto[match_i][0][i]; //i;
          cout << "MSD match " << event.Id << " " 
            << event.utctime() - utcChisto[match_i] << " "
            << event.utctime() - 315964782 << " " << cnt << endl;
        }
      }
    }
    cout << "MSD " << event.Id << endl;
    sdOutPutFile.Write(event);
  }
  
  sdOutPutFile.Close();
  
  return 0;
}
