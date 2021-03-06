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
   if ( argc < 5 ) {
		 cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <PMT>  <month> <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <PMT>: ID for the PMT you want to analyse" << endl
         << "  <month>: Month of the date taken " << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl
				 << " " << endl
				 << "In case you want the distribution of all events for a specific Station, " << endl
				 << "just make sure the stationsFile conteins a single station." << endl
				 << endl;
 		 exit(0);
	 }
	
  const char* stationsFileName = argv[1];
  const char* whichpmt = argv[2];
  const char* whichmonth = argv[3];
  AugerIoSd input(argc-4, argv+4);
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
  TString monthdata = whichmonth;
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
  
	cerr << endl;
	cerr << "========================" << endl;
  cerr << " You have selected " << pmtname << endl;
	cerr << endl;

  TFile hfile("ubFitaop"+pmtname+monthdata+".root","RECREATE","");

	unsigned int totSt = stationsIds.size();
  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
	bool readSglEvt = false;

	sort(stationsIds.begin(), stationsIds.end());
	double stationsBins[totSt];
	for ( unsigned int i=0; i<totSt; i++)
	{
		 stationsBins[i] = stationsIds[i];
		 cout << i << " -> " << stationsIds[i] << endl;
	}

	readHistos fitHist;
	TString tmpName;

	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge
	TH1F *apHbase = new TH1F("apHbase","",totSt, 0, totSt);//It stores A/P from HBase for St.
	TH1F *apCalib = new TH1F("apCalib","",totSt, 0, totSt);//It stores A/P from Calib for St.
  
	TH1F *chisHbasePk = new TH1F("chisHbasePk","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisCalibPk = new TH1F("chisCalibPk","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisHbaseCh = new TH1F("chisHbaseCh","",totSt, 0, totSt);//It stores chis from fitted charge for St.
	TH1F *chisCalibCh = new TH1F("chisCalibCh","",totSt, 0, totSt);//It stores chis from fitted charge for St.
	TH1F *chisHbasePkbad = new TH1F("chisHbasePkbad","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisCalibPkbad = new TH1F("chisCalibPkbad","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisHbaseChbad = new TH1F("chisHbaseChbad","",totSt, 0, totSt);//It stores chis from fitted charge for St.
	TH1F *chisCalibChbad = new TH1F("chisCalibChbad","",totSt, 0, totSt);//It stores chis from fitted charge for St.

	TH1F *eventStat = new TH1F("eventStat","",totSt, 0, totSt);//It stores the number of events per St.
	unsigned int evtTime = 0; //[totSt]; // Get the Event Id for the respective entry and St.

  TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("receCh", "TH1F", &receCh);
	treeHist->Branch("recePk", "TH1F", &recePk);
	treeHist->Branch("apHbase", "TH1F", &apHbase);
	treeHist->Branch("apCalib", "TH1F", &apCalib);

	treeHist->Branch("chisHbasePk", "TH1F", &chisHbasePk);
	treeHist->Branch("chisCalibPk", "TH1F", &chisCalibPk);
	treeHist->Branch("chisHbaseCh", "TH1F", &chisHbaseCh);
	treeHist->Branch("chisCalibCh", "TH1F", &chisCalibCh);
	treeHist->Branch("chisHbasePkbad", "TH1F", &chisHbasePkbad);
	treeHist->Branch("chisCalibPkbad", "TH1F", &chisCalibPkbad);
	treeHist->Branch("chisHbaseChbad", "TH1F", &chisHbaseChbad);
	treeHist->Branch("chisCalibChbad", "TH1F", &chisCalibChbad);
  
	treeHist->Branch("eventStat", "TH1F", &eventStat);
	treeHist->Branch("evtTime",&evtTime,"evtTime/I");

  EventPos pos;
	double blCorrHbase = 0.;
	double blCorrCalib = 0.;
  TH1F *tmp = new TH1F();
	int currentDaySt[totSt];
  fitHist.getGraph = false;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
    IoSdEvent event(pos);

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i)
		{
      found = false;
      for ( vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter )
        if ( event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue; 
      
      if ( event.Stations[i].Error ) 
        continue;
      cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
        << " " << nrEventsRead-1 << endl;
 
      for ( unsigned int stid = 0; stid<totSt; stid++ )
        if ( stationsBins[stid] == event.Stations[i].Id && 
            currentDaySt[stid] != event.UTCTime )
        {
          // ================
          // *** Baseline ***
          blCorrHbase = event.Stations[i].HBase(pmtId-1)->GetMean();
          blCorrCalib = event.Stations[i].Calib->Base[pmtId-1];
          
          // =================
          // *** For HBase ***
          recePk = event.Stations[i].HPeak(pmtId-1);
          receCh = event.Stations[i].HCharge(pmtId-1);
          tmpName.Form("%d%d%d", event.UTCTime,stid, nrEventsRead-1);
              
          fitHist.getPkCrrOst(*recePk, blCorrHbase, tmpName+"Hbpk");
          tmp = fitHist.getPkCorr();
          fitHist.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
          
          fitHist.getChCrrOst(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
          tmp = fitHist.getChCorr();
          fitHist.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);

    			// ============================
		  		// *** Saving A/P for HBase ***
		    	if ( fitHist.fitChOk && fitHist.fitPkOk )
          {
            apHbase->SetBinContent( stid, fitHist.vemPosCh/fitHist.vemPosPk );
            chisHbasePk->SetBinContent( stid, fitHist.chisPeak );
            chisHbaseCh->SetBinContent( stid, fitHist.chisCharge );
          }
			    else
          {
            apHbase->SetBinContent( stid, 0. );
            chisHbasePkbad->SetBinContent( stid, fitHist.chisPeak );
            chisHbaseChbad->SetBinContent( stid, fitHist.chisCharge );
          }

          // =================
          // *** For Calib ***
          fitHist.getPkCrrOst(*recePk, blCorrCalib, tmpName+"Clpk");
          tmp = fitHist.getPkCorr();
          fitHist.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
          
          fitHist.getChCrrOst(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Clch");
          tmp = fitHist.getChCorr();
          fitHist.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);

          // ============================
		    	// *** Saving A/P for Calib ***
			  	if ( fitHist.fitChOk && fitHist.fitPkOk )
          {
            apCalib->SetBinContent( stid, fitHist.vemPosCh/fitHist.vemPosPk );
            chisCalibPk->SetBinContent( stid, fitHist.chisPeak );
            chisCalibCh->SetBinContent( stid, fitHist.chisCharge );
          }
          else
          {
            apCalib->SetBinContent( stid, 0. );
            chisCalibPkbad->SetBinContent( stid, fitHist.chisPeak );
            chisCalibChbad->SetBinContent( stid, fitHist.chisCharge );
          } 

          eventStat->Fill( stid, 1 );
          evtTime = event.UTCTime;
          treeHist->Fill();
          currentDaySt[stid] = event.UTCTime;
          break;
        }
    }
		if (readSglEvt)
			break; // Read only one event
	}

  hfile.Write();
  hfile.Close();
	return 0;
}
