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

#include "stats.h"
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

  TFile hfile("uubCalibHist"+pmtname+".root","RECREATE","");

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
	bool readSglEvt = false;

	readHistos fitHist;
	TString tmpName;

  TH1F *rawPk = new TH1F ("rawPk","",150, 0, 150); // Receive Pk Raw from IoSdHisto::Peak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge
	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *corrCh = new TH1F(); // Receive Ch corrected for offset
	TH1F *corrPk = new TH1F(); // Receive Pk corrected for offset
	TGraphErrors *fllHCh = new TGraphErrors();
	TGraphErrors *fllHPk = new TGraphErrors();

	double baselineHbase = 0.;
	double baselineCalib = 0.;
	int offSetCh = 0; // Receive offset from IoSdHisto::Histo
	int offSetPk = 0; // Receive offset from IoSdHisto::Histo
	double ap = 0; // Receive area/peak
	unsigned int entryEvt = 0; // Get the Event Id for the respective entry
	
	TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("rawPk", "TH1F", &rawPk);
	treeHist->Branch("receCh", "TH1F", &receCh);
	treeHist->Branch("recePk", "TH1F", &recePk);
	treeHist->Branch("corrPk", "TH1F", &corrPk);
	treeHist->Branch("corrCh", "TH1F", &corrCh);
	treeHist->Branch("fllHCh", "TGraphErrors", &fllHCh);
	treeHist->Branch("fllHPk", "TGraphErrors", &fllHPk);

	treeHist->Branch("baselineHbase", &baselineHbase, "baselineHbase/D");
	treeHist->Branch("baselineCalib", &baselineCalib, "baselineCalib/D");
	treeHist->Branch("offSetCh",&offSetCh,"offSetCh/I");
	treeHist->Branch("offSetPk",&offSetPk,"offSetPk/I");
	treeHist->Branch("ap",&ap,"ap/D");
	treeHist->Branch("entryEvt",&entryEvt,"entryEvt/I");

  EventPos pos;

	double baseCorrApp = 0.;

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
      
      if ( event.Stations[i].IsUUB ) {
				cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
					<< " " << nrEventsRead-1 << endl;

				if (event.Stations[i].Error==256) 
				{
					//readSglEvt = true;         
          /*
          cout << "MSD: " 
            << event.Id << " " 
            << event.Stations[i].Histo->Peak[pmtId-1][0] << endl;
            */
          for (int bn = 0; bn<150; bn++)
            rawPk->SetBinContent( bn, event.Stations[i].Histo->Peak[pmtId-1][bn] );

					tmpName.Form("%d", nrEventsRead-1);
					receCh = event.Stations[i].HCharge(pmtId-1);
					recePk = event.Stations[i].HPeak(pmtId-1);
					offSetCh = event.Stations[i].Histo->Offset[pmtId-1+6];
					offSetPk = event.Stations[i].Histo->Offset[pmtId-1+3];
					fitHist.getGraph = true;
					baseCorrApp = event.Stations[i].HBase(pmtId-1)->GetMean();
					//baseCorrApp = event.Stations[i].Calib->Base[pmtId-1];

					// ==================
					// *** Fitting Pk ***
					fitHist.getPkCrrOst(*recePk, baseCorrApp, tmpName+"pk");
					corrPk = fitHist.getPkCorr();
					fitHist.getFitPk(*corrPk, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
					if ( fitHist.fitPkOk )
						fllHPk = fitHist.getFitGraphPk();

					// ==================
					// *** Fitting Ch ***
					fitHist.getChCrrOst(*receCh, baseCorrApp*0., tmpName+"ch");
					corrCh = fitHist.getChCorr();
					fitHist.getFitCh(*corrCh, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);
					if ( fitHist.fitChOk )
						fllHCh = fitHist.getFitGraphCh();

					// ==================
					// *** Saving A/P ***
					if ( fitHist.fitChOk && fitHist.fitPkOk )
          {
						ap = fitHist.vemPosCh/fitHist.vemPosPk;
            if ( ap > 8 )
              cout << "MSD " << event.Id << " " << ap << endl;
          }
					else
						ap = 0.;

					baselineHbase = event.Stations[i].HBase(pmtId-1)->GetMean();
					baselineCalib = event.Stations[i].Calib->Base[pmtId-1];

					entryEvt = event.UTCTime;
					treeHist->Fill();
				}
			}
		}
		if (readSglEvt)
			break; // Read only one event
	}

  hfile.Write();
  hfile.Close();
	return 0;
}
