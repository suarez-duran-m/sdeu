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

#include "../stats.h"
#include "../readHistos.h"

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
	double stationsBins[totSt];
	unsigned int checkDiffOffsetC[totSt]; // Check if diff offset and firstbin is 0 Charge
	unsigned int checkDiffOffsetP[totSt]; // Check if diff offset and firstbin is 0 Peak 
		
  sort(stationsIds.begin(), stationsIds.end());
  for ( unsigned int i=0; i<totSt; i++ ) {
    stationsBins[i] = stationsIds[i];
		checkDiffOffsetC[i] = 0;
		checkDiffOffsetP[i] = 0;
    cout << stationsBins[i] << " -> " << stationsIds[i] << endl;
  }

	//unsigned int totDays = 121; // From 1st December, 2020 to 31st March, 2021 

  TFile hfile("uubCalibHistAll"+pmtname+".root","RECREATE","");

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
	bool readSglEvt = false;

	TH1F *setCh = new TH1F (); // Receive Charge from IoSdStation::HCharge
	TH1F *setPk = new TH1F (); // Receive Peak from IoSdStation::HPeak
	TH1F *setChHisto = new TH1F ("setChHisto", "", 600, 0, 600);//Receive Charge from IoSdHisto::Charge
	TH1F *setPkHisto = new TH1F ("setPkHisto", "", 150, 0, 150); // Receive Peak from IoSdHisto::Peak
	TH1F *base = new TH1F (); // Receive Baselinefrom IoSdStation::Histo

	TH1F *diffOffsetStationP = new TH1F("diffOffsetStationP", 
			"Difference between Offset and first bin UUB Peak Histogram "+pmtname,
			totSt, 0, totSt);
	TH1F *diffOffsetStationC = new TH1F("diffOffsetStationC", 
			"Difference between Offset and first bin UUB Charge Histogram "+pmtname,
			totSt, 0, totSt);


	int offSetCh = 0; // Receive offset from IoSdHisto::Histo
	int offSetPk = 0; // Receive offset from IoSdHisto::Histo
	unsigned int entryEvt = 0; // Get the Event Id for the respective entry
	
	TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("setCh", "TH1F", &setCh);
	treeHist->Branch("setPk", "TH1F", &setPk);
	treeHist->Branch("setChHisto", "TH1F", &setChHisto);
	treeHist->Branch("setPkHisto", "TH1F", &setPkHisto);
	treeHist->Branch("base", "TH1F", &base);
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
      
      if ( event.Stations[i].IsUUB ) {
				cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
					<< " " << nrEventsRead-1 << endl;

				if (event.Stations[i].Error==256) { //0+256
					//readSglEvt = true;
					setCh = event.Stations[i].HCharge(pmtId-1);
					setPk = event.Stations[i].HPeak(pmtId-1);
					offSetCh = event.Stations[i].Histo->Offset[pmtId-1+6];
					offSetPk = event.Stations[i].Histo->Offset[pmtId-1+3];
					base = event.Stations[i].HBase(pmtId-1);

					for ( unsigned int id=0; id<totSt; id++ )
						if ( stationsBins[id] == event.Stations[i].Id ) {
							if ( setPk->GetBinLowEdge(1) - offSetPk != 0 )
								checkDiffOffsetP[id] += fabs(setPk->GetBinLowEdge(1) - offSetPk);
							if ( setCh->GetBinLowEdge(1) - offSetCh != 0 )
								checkDiffOffsetC[id] += fabs(setCh->GetBinLowEdge(1) - offSetCh);
						}
					for ( unsigned b=1; b<601; b++ ) {
						setChHisto->Fill(b, event.Stations[i].Histo->Charge[pmtId-1][b]);
						if (b < 150)
							setPkHisto->Fill(b, event.Stations[i].Histo->Peak[pmtId-1][b]);
					}
					entryEvt = event.Id;		
					treeHist->Fill();
				}
			}
		}
		if (readSglEvt)
			break; // One event
	}
	for ( unsigned int idi=0; idi<totSt; idi++ ) {
		diffOffsetStationP->Fill( idi, checkDiffOffsetP[idi] );
		diffOffsetStationC->Fill( idi, checkDiffOffsetC[idi] );
	}

  hfile.Write();
  hfile.Close();
	return 0;
}
