/* 
 * ============================================
 * Code to understand the Offset value in terms
 * of: Baseline-HBase and Baseline-Calib and 
 * temperature.
 * ============================================
 *
 */

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

#include "readHistos.h"

using namespace std;

// ========================== 
// ******** The MAIN ********
// ==========================
int main (int argc, char *argv[]) {
   if ( argc < 4 ) {
		 cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <PMT>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
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
  
  TString nameStati = to_string( stationsIds[0] );
  TString pmtname = whichpmt;
  int pmtId= atoi( pmtname );
  if ( pmtId > 0 && pmtId < 4 ){
     pmtname = "PMT"+to_string( pmtId );
  }
  else if ( pmtId == 4 )
    pmtname = "SPMT";
  else if ( pmtId == 5 )
    pmtname = "PMTSSD";
  else{
    cout << "==================================================" << endl;
    cout << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
    cout << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " 
      << "4 For SPMT; " << "5 For PMTSSD" << endl;
    cout << "==================================================" << endl;
    exit(0);
	}
  
  cerr << "You have selected " << pmtname << endl;

	unsigned int totSt = stationsIds.size();

	if ( totSt > 1 ) {
		cerr << "This codes if for single station." << endl;
		exit(0);
	}

	pmtname += "St"+to_string( stationsIds[0] );

  TFile hfile("ubOffset"+pmtname+".root","RECREATE","");

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
	bool readSglEvt = false;

	readHistos fitHist;
	TString tmpName;	

	double baselineHbase = 0.; // Save Baseline-Hbase
	double baselineCalib = 0.; // Save Baseline-Calib
	int offSetCh = 0; // Save offset from IoSdHisto::Histo for Ch
	int offSetPk = 0; // Save offset from IoSdHisto::Histo for Pk 
	unsigned int entryEvt = 0; // Save EventId as entry
	
	TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("baselineHbase", &baselineHbase, "baselineHbase/D");
	treeHist->Branch("baselineCalib", &baselineCalib, "baselineCalib/D");
	treeHist->Branch("offSetCh",&offSetCh,"offSetCh/I");
	treeHist->Branch("offSetPk",&offSetPk,"offSetPk/I");
	treeHist->Branch("entryEvt",&entryEvt,"entryEvt/I");

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
    IoSdEvent event(pos);

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i){
      found = false;
      for ( vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter )
        if ( event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;

			if (event.Stations[i].Error)
				continue;

			cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
				<< " " << nrEventsRead-1 << endl;

			//readSglEvt = true;
			offSetCh = event.Stations[i].Histo->Offset[pmtId-1+6];
			offSetPk = event.Stations[i].Histo->Offset[pmtId-1+3];

			baselineHbase = event.Stations[i].HBase(pmtId-1)->GetMean();
			baselineCalib = event.Stations[i].Calib->Base[pmtId-1];

			entryEvt = event.UTCTime;
			treeHist->Fill();

			if (readSglEvt)
				break; // Read only one event
		}
	}

  hfile.Write();
  hfile.Close();
	return 0;
}