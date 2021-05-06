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

#include "fitpeak.h"
#include "fitcharge.h"

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
  
	cerr << endl;
	cerr << "========================" << endl;
  cerr << " You have selected " << pmtname << endl;
	cerr << endl;

  TFile hfile("kkuubAoPsgl"+pmtname+"St1223.root","RECREATE","");

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

	fitpeak fitPk;
	fitcharge fitCh;
	TString tmpName;

	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
  //TH1F *smoothPk = new TH1F (); // Just to check how the smooth likes
  //TH1F *smoothCh = new TH1F (); // Just to check how the smooth likes
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge
	TH1F *apHbase = new TH1F("apHbase","",totSt, 0, totSt);//It stores A/P from HBase for St.
	TH1F *apCalib = new TH1F("apCalib","",totSt, 0, totSt);//It stores A/P from Calib for St.
  
	TH1F *chisHbasePk = new TH1F("chisHbasePk","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisCalibPk = new TH1F("chisCalibPk","",totSt, 0, totSt);//It stores chis from fitted peak for St.
	TH1F *chisCh = new TH1F("chisCh","",totSt, 0, totSt);//It stores chis from fitted charge for St.

  TH1F *okFitPkHb = new TH1F("okFitPkHb","",totSt, 0, totSt);//It stores the number of Fit Ok for charge per St. 
  TH1F *okFitPkCa = new TH1F("okFitPkCa","",totSt, 0, totSt);//It stores the number of Fit Ok for charge per St.
  TH1F *okFitCh = new TH1F("okFitCh","",totSt, 0, totSt);//It stores the number of Fit Ok for charge per St.
  TH1F *nslopsHbPk = new TH1F("nslopsHbPk","",totSt, 0, totSt);//It stores number slops in Pk Histo. per St, for Hbase.
  TH1F *nslopsCaPk = new TH1F("nslopsCaPk","",totSt, 0, totSt);//It stores number slops Pk Histo. per St, for Calib.
  TH1F *nslopsCh = new TH1F("nslopsCh", "",totSt, 0, totSt);

	TH1F *eventStat = new TH1F("eventStat","",totSt, 0, totSt);//It stores the number of events per St.
	unsigned int evtTime = 0; //[totSt]; // Get the Event Id for the respective entry and St.

  TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("receCh", "TH1F", &receCh);
	//treeHist->Branch("smoothCh", "TH1F", &smoothCh);
	treeHist->Branch("recePk", "TH1F", &recePk);
	//treeHist->Branch("smoothPk", "TH1F", &smoothPk);
	treeHist->Branch("apHbase", "TH1F", &apHbase);
	treeHist->Branch("apCalib", "TH1F", &apCalib);
  
	treeHist->Branch("chisHbasePk", "TH1F", &chisHbasePk);
	treeHist->Branch("chisCalibPk", "TH1F", &chisCalibPk);
	treeHist->Branch("chisCh", "TH1F", &chisCh);
	
  treeHist->Branch("okFitPkHb", "TH1F", &okFitPkHb);
  treeHist->Branch("okFitPkCa", "TH1F", &okFitPkCa);
	treeHist->Branch("okFitCh", "TH1F", &okFitCh);

  treeHist->Branch("nslopsHbPk","TH1F",&nslopsHbPk);
  treeHist->Branch("nslopsCaPk","TH1F",&nslopsCaPk);
  treeHist->Branch("nslopsCh","TH1F",&nslopsCh);
  
	treeHist->Branch("eventStat", "TH1F", &eventStat);
	treeHist->Branch("evtTime",&evtTime,"evtTime/I");

  EventPos pos;
	double blCorrHbase = 0.;
	double blCorrCalib = 0.;
  TH1F *tmp = new TH1F();
	int currentDaySt[totSt];
  fitPk.getGraph = false;

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
      
      if ( event.Stations[i].IsUUB ) 
			{
				cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
					<< " " << nrEventsRead-1 << endl;

				if (event.Stations[i].Error==256) 
				{
          //readSglEvt = false;
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
              tmpName.Form("%d%d%d", event.UTCTime,stid, nrEventsRead-1);
              
              fitPk.getCrrSmooth(*recePk, blCorrHbase, tmpName+"Hbpk");
              tmp = fitPk.getPkCorrSmooth();
              //smoothPk = fitPk.getPkCorr();
              //recePk = fitPk.getPkCorr2();

              if ( fitPk.getValidHisto(*tmp) >= 3 )
              {
                nslopsHbPk->SetBinContent( stid, fitPk.getValidHisto(*tmp) );
                fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
              }
              if ( fitPk.fitPkOk )
                okFitPkHb->Fill( stid, 1 );

              receCh = event.Stations[i].HCharge(pmtId-1);
              fitCh.getChCrrSmooth(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
              tmp = fitCh.getChCorrSmooth();
              //smoothCh = fitCh.getChCorrSmooth();
              //receCh = fitCh.getChCorr2();
              if ( fitCh.getValidHisto(*tmp) >= 2 )
              {
                nslopsCh->SetBinContent( stid, fitCh.getValidHisto(*tmp) );
                fitCh.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);
              }
              if ( fitCh.fitChOk )
                okFitCh->Fill( stid, 1 );

    					// ============================
		    			// *** Saving A/P for HBase ***
				    	if ( fitCh.fitChOk && fitPk.fitPkOk )
              {
						    apHbase->SetBinContent( stid, fitCh.vemPosCh/fitPk.vemPosPk );
                //chisHbasePk->SetBinContent( stid, fitPk.chisPeak );
              }
					    else
						    apHbase->SetBinContent( stid, 0. );
             
              // =================
              // *** For Calib ***
              fitPk.getCrrSmooth(*recePk, blCorrCalib, tmpName+"Capk");
              tmp = fitPk.getPkCorrSmooth();
              if ( fitPk.getValidHisto(*tmp) >= 3 )
              {
                nslopsCaPk->SetBinContent( stid, fitPk.getValidHisto(*tmp) );
                fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
              }
              if ( fitPk.fitPkOk )
                okFitPkCa->Fill( stid, 1 );  

    					// ============================
		    			// *** Saving A/P for Calib ***
				    	if ( fitCh.fitChOk && fitPk.fitPkOk )
              {
						    apCalib->SetBinContent( stid, fitCh.vemPosCh/fitPk.vemPosPk );
                //chisCalibPk->SetBinContent( stid, fitPk.chisPeak );
              }
					    else
						    apCalib->SetBinContent( stid, 0. );

              eventStat->Fill( stid, 1 );
							evtTime = event.UTCTime;
							treeHist->Fill();
							currentDaySt[stid] = event.UTCTime;
							break;
						}
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
