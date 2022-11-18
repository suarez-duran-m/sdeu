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
  int minArg = 6;
  if ( argc < minArg ) {
    cout << endl << "=========================" << endl << endl
      << "Usage: " << argv[0] << " sd*packs_file ad_output_name ad_files" << endl;
    cout << endl 
      << "sd*packs_file: file with coincidence histos" << endl
      << "ad_output_name: filename for new ad*.root with coincidence histos" << endl
      << "ad_files: ad*.root files where coincidence histos will be added." << endl 
      << endl;
    exit(0);
  } 

  const char *outPutADfileName = argv[4];
  // 
  // Charging coincidence histograms
  vector < int > utcChisto;
  vector < int > stChisto;
  vector < vector < vector < int > > > cQhisto;
  vector < vector < vector < int > > > cHeigth;
  for(int i=1; i<4; i++) {
    char *fileWithCoinci = argv[i];
    rawCoinciHistoData rawCoinciHistoData;
    rawCoinciHistoData.readData( fileWithCoinci );
    for(int i=0; i<rawCoinciHistoData.getUtcEvtWidthChisto().size(); i++) {
      utcChisto.push_back( rawCoinciHistoData.getUtcEvtWidthChisto()[i] );
      stChisto.push_back( rawCoinciHistoData.getStWidthChisto()[i] );
      cQhisto.push_back( rawCoinciHistoData.getCQhisto()[i] );
      cHeigth.push_back( rawCoinciHistoData.getCheight()[i] );
    }
    cout << "MSD0 rawUTCs: " << rawCoinciHistoData.getUtcEvtWidthChisto().size() << endl;
    cout << "MSD0 utcChis: " << utcChisto.size() << endl;
    rawCoinciHistoData.SetClear();
    cout << "MSD1 rawUTCs: " << rawCoinciHistoData.getUtcEvtWidthChisto().size() << endl;
    cout << "MSD1 utcChis: " << utcChisto.size() << endl;
  }    
  //
  // Reading ad*.root files  
  AugerFile adFileInput( argc-(minArg-1), argv+(minArg-1) );
  const unsigned int totalNrEvents = adFileInput.NumberOfEvents();
  AugerFile adOutPutFile;
  adOutPutFile.Open(outPutADfileName, AugerFile::eWrite);
  IoSdHistoCoinci *newCoinciHistos;
  unsigned int nrEventsRead = 0;
  AugerEvent theAugerEvent;
  // Moving through events
  //
  while ( adFileInput.ReadNext( theAugerEvent ) == AugerFile::eSuccess ) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0)
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
    // Ignoring FD events
    // 
    if ( !theAugerEvent.HasSd() ) {
      adOutPutFile.Write(theAugerEvent, false);
      continue;
    }
    //
    // Checking if current evt-UTC match with some CHisto-UTC
    IoSdEvent &event = theAugerEvent.Sd();
    vector < int > indexForUtcMatch;
    bool utcMatched = false;
    for ( unsigned int utc_pos=0; utc_pos<utcChisto.size(); utc_pos++ ) {
      if ( fabs(utcChisto[utc_pos] - event.utctime()) < 3 ) {
        indexForUtcMatch.push_back( utc_pos );
        utcMatched = true;
      }
    }
    //
    // Ignoring unmatched events
    if ( !utcMatched ) {
      adOutPutFile.Write(theAugerEvent, false);
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
        //cout << "# Event " << event.Id << " Station " << event.Stations[i].Id << endl;
        cout << "MSDtime " << event.utctime() << " " << event.Stations[i].Id << endl;
        // 
        // Writing histo into ascii file for cross-check
        /*
        if ( event.Stations[evtSt_i].Id==1191 ) {
          for ( int pmt_i=0; pmt_i<3; pmt_i++ ) {
            fstream outPutHistoCrossCheck;
            TString outputName = Form("%ld_%d_%d", event.utctime()-315964782, event.Stations[evtSt_i].Id, pmt_i);
            outPutHistoCrossCheck.open(outputName+".dat", ios_base::out);
            //outPutHistoCrossCheck.open("histos_crossCheck/"+outputName+".dat", ios_base::out);
            for ( auto &i : cQhisto[match_i][pmt_i] )
              outPutHistoCrossCheck << i << endl;
            outPutHistoCrossCheck.close();
          }
        }
        */
        //
        // Storing coincidence histos into current event        
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
      }
    }
    //cout << "MSD " << event.Id << endl;
    adOutPutFile.Write(theAugerEvent, false);
  }
  adOutPutFile.Close();
  return 0;
}
