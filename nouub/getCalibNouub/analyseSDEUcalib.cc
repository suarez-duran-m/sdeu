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
	
	stats getmrms;
	readHistos rdHist;

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

  const unsigned int nblbins = 100;
  double meanf = 0.;
  double meanl = 0.;
  double rmsf = 0.;
  //double rmsl = 0.;

  TFile hfile("calibHist"+pmtname+".root","RECREATE","");

  unsigned int totSt = stationsIds.size();
	double chok = 0.;
	double pkok = 0.;

  vector < int > blpmth;
	vector < vector < unsigned int > > stckEvt;
  stckEvt.resize(totSt);
  vector <  double > stckVch;
  stckVch.resize(totSt);
  vector <  double > stckVpk;
  stckVpk.resize(totSt);
	vector < double > stckAP;
	stckAP.resize(totSt);
	
  sort(stationsIds.begin(), stationsIds.end());
  Double_t stationsBins[totSt];
  for ( unsigned int i=0; i<totSt; i++){
    stationsBins[i] = stationsIds[i];
		stckVch[i] = 0.;
		stckVpk[i] = 0.;
		stckAP[i] = 0.;
		for ( int j=0; j<3; j++ )
			stckEvt[i].push_back(0);
    cout << i << " -> " << stationsIds[i] << endl;
  }

  unsigned int previusEvent = 0;
  unsigned int sameUtc = 0;
	unsigned int nday = 0;
	unsigned int cday = 1599048000;
	unsigned int dday = 86400;
	unsigned int totDays = 90; // From 1st September, 2020 to 30th November, 2020

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;

	TH2F hCh ("hCh", "Average VEM-Charge per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt); 
	TH2F hPk ("hPk", "Average VEM-Peak per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt); 
	TH2F hap ("hap", "Average VEM-Charge/VEM-Peak per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt);

	unsigned int histBins = 12000;
	TH1F *setCh = new TH1F ("setCh", "", histBins, 0, histBins );
	TH1F *setPk = new TH1F ("setPk", "", histBins, 0, histBins );
	
	TH1F sglSt ("sglSt", "", totalNrEvents, 0, totalNrEvents);

	//unsigned int nentry = 0;
	
	TTree *treeCh = new TTree("Charge","");
	TGraphErrors *fllHCh = new TGraphErrors ();
	treeCh->Branch("fllHCh", "TGraphErrors", &fllHCh);
	TTree *treePk = new TTree("Peak","");
	TGraphErrors *fllHPk = new TGraphErrors ();
	treePk->Branch("fllHPk", "TGraphErrors", &fllHPk);

	if ( totSt==1 )
		rdHist.getGraph = true;

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
		TEcEvent event(pos);
    //IoSdEvent event(pos);

    if ( event.Id == previusEvent ) //&& sameUtc == event.utctime() )
      continue;

    previusEvent = event.Id;
    sameUtc = event.utctime();

    if ( sameUtc > cday ) {
      for ( unsigned int id=0; id<totSt; id++) {
        if ( stckEvt[id][0] > 0 )
					hCh.Fill( nday, id, stckVch[id] / stckEvt[id][0] );
				if ( stckEvt[id][1] > 0 )
					hPk.Fill( nday, id, stckVpk[id] / stckEvt[id][1] );
				if ( stckEvt[id][2] > 0 )
					hap.Fill( nday, id, stckAP[id] / stckEvt[id][2] );
				stckVch[id] = 0.;
				stckVpk[id] = 0.;
				stckAP[id] = 0.;
				for ( int j=0; j<3; j++ )
					stckEvt[id][j] = 0;
      }
      cday += dday;
      nday++;
    }

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i) {
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter )
				if (event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;

			cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
				<< " " << nrEventsRead-1 << " " << sameUtc << endl;
	
			if(event.fCalibStations[i].Error)
				continue;
			if(event.fCalibStations[i].fSigInVEM < 0)
				continue;

			blpmth.clear();

			for (unsigned int ii=0; ii<event.Stations[i].Fadc->NSample; ii++)
				blpmth.push_back( 0 );

			for (unsigned int k=0;k<event.Stations[i].Fadc->NSample; k++)
				blpmth[k] = event.Stations[i].Fadc->Trace[pmtId-1][0][k];

			for ( unsigned int id=0; id<totSt; id++ )
				if ( stationsBins[id] == event.Stations[i].Id ) {
					chok = 0.;
					pkok = 0.;
          meanf = getmrms.getMean(blpmth, nblbins, true);
          meanl = getmrms.getMean(blpmth, nblbins, false);
					rmsf = getmrms.getRms(blpmth, meanf, nblbins, true);
          if ( fabs(meanl-meanf) < 2*rmsf ) {
						setCh = event.Stations[i].HCharge(pmtId-1);
						rdHist.getFullFit( *setCh, true, 0.1, 10 ); // 20% EMpeak and 10 bins forward from this last.
						if ( rdHist.fitChOk ) {
							stckVch[id] += rdHist.vemPos;
							stckEvt[id][0]++;
							chok = rdHist.vemPos;
							if ( totSt==1 ) {
								fllHCh = rdHist.getFitGraph();
								treeCh->Fill();
							}
						}
						setCh->Reset();
						setPk = event.Stations[i].HPeak(pmtId-1);
						rdHist.getFullFit( *setPk, false, 0.1, 5 ); // 20% EMpeak and 5 bins forward from this last.
						//nentry++;
						if ( rdHist.fitPkOk ) {
							stckVpk[id] += rdHist.vemPos;
							stckEvt[id][1]++;
							pkok = rdHist.vemPos;
							if ( totSt==1 ) {
								fllHPk = rdHist.getFitGraph();
								treePk->Fill();
							}
						}
						if ( rdHist.fitChOk && rdHist.fitPkOk ) {
							stckAP[id] += chok/pkok;
							stckEvt[id][2]++;
							sglSt.Fill( nrEventsRead, chok/pkok );
						}
						setPk->Reset();
						}
						break;
					}
		}
	}

  hfile.Write();
  hfile.Close();
	return 0;
}
