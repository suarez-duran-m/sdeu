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

	TH1F *setCh = new TH1F (); // Receive Charge from IoSdStation::HCharge
	TH1F *setPk = new TH1F (); // Receive Peak from IoSdStation::HPeak
	TH1F *setChHisto = new TH1F ("setChHisto", "", 600, 0, 600); // Receive Charge from IoSdHisto::Charge
	TH1F *setPkHisto = new TH1F ("setPkHisto", "", 150, 0, 150); // Receive Peak from IoSdHisto::Peak
	TH1F *base = new TH1F (); //"base", "", 1000, 0, 1000); // Receive Baselinefrom IoSdHisto::Histo

	TH1F *pkCorrBl = new TH1F(); // Correcion for Baseline
	TH1F *pkCorrOff = new TH1F(); // Correction for Offset 

	int offSetCh = 0; // Receive offset from IoSdHisto::Histo
	int offSetPk = 0; // Receive offset from IoSdHisto::Histo
	unsigned int entryEvt = 0; // Get the Event Id for the respective entry
	
	TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("setCh", "TH1F", &setCh);
	treeHist->Branch("setPk", "TH1F", &setPk);
	treeHist->Branch("setChHisto", "TH1F", &setChHisto);
	treeHist->Branch("setPkHisto", "TH1F", &setPkHisto);
	treeHist->Branch("base", "TH1F", &base);
	treeHist->Branch("pkCorrBl", "TH1F", &pkCorrBl);
	treeHist->Branch("pkCorrOff", "TH1F", &pkCorrOff);
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
					for ( unsigned b=1; b<601; b++ ) {
						setChHisto->Fill(b, event.Stations[i].Histo->Charge[pmtId-1][b]);
						if (b < 150)
							setPkHisto->Fill(b, event.Stations[i].Histo->Peak[pmtId-1][b]);
					}
					if ( nrEventsRead-1==1681 ) {
						double x[151];
						TH1F *tmp;
						for ( unsigned b=0; b<150; b++ )
							x[b] = setPk->GetBinLowEdge(b+1)-284; // Correction Bl

						x[150] = x[149] + setPk->GetBinWidth(1);
						tmp = new TH1F("tmp", "ok", 150, x);
						for ( unsigned b=0; b<150; b++ )
							tmp->SetBinContent(b+1, setPk->GetBinContent(b+1));
						pkCorrBl = tmp;

						TH1F *tmp1;
						for ( unsigned b=0; b<150; b++ )
							x[b] = setPk->GetBinLowEdge(b+1)-273; // Correction Offset
						x[150] = x[149] + setPk->GetBinWidth(1);
						tmp1 = new TH1F("tmp1", "ok", 150, x);
						for ( unsigned b=0; b<150; b++ )
							tmp1->SetBinContent(b+1, setPk->GetBinContent(b+1));
						pkCorrOff = tmp1;
					}

					entryEvt = event.Id;				
					treeHist->Fill();
				}
			}
		}
		if (readSglEvt)
			break; // One event
	}

  hfile.Write();
  hfile.Close();
	return 0;
}
