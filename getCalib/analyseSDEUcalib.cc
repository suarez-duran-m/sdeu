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
	unsigned int totSt = stationsIds.size();
	double chok = 0.;
	double pkok = 0.;
	double base = 0.;

	if ( totSt==1 )
		pmtname += "St"+to_string( stationsIds[0] );

  TFile hfile("uubCalibHist"+pmtname+".root","RECREATE","");

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
	unsigned int cday = 1606867200;
	unsigned int dday = 86400;
	unsigned int totDays = 90; // From 1st December, 2020 to 28th February, 2021

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
	double ubUub = .38; //8.33/25.0; // 8.33/25.0 From UUB to UB

	TH1F *charge = new TH1F ();// ("charge", "Charge for single station", totalNrEvents, 0, totalNrEvents);

	TH2F hCh ("hCh", "Average VEM-Charge per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt); // Average Charge per day - per station
	TH2F hPk ("hPk", "Average VEM-Peak per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt); // Average Peak per day - per station
	TH2F hap ("hap", "Average VEM-Charge/VEM-Peak per day for "+pmtname, totDays, 0, totDays, totSt, 0, totSt); // Average A/P per day - per station
	TH2F apDist ("apDist", "Average VEM-Charge/VEM-Peak per day for "+pmtname, totSt, 0, totSt, totalNrEvents, 0, totalNrEvents); // A/P per station for all events.

	TH1F sglSt ("sglSt", "", totalNrEvents, 0, totalNrEvents); // A/P for all events at specific station
	TH1F stckHc ("stckHc", "Average Charge Histogram", 5000, 0, 5000); // Charge for all events at specific station
	TH1F stckHp ("stckHp", "Average Peak Histogram", 5000, 0, 5000); // Peak for all events at specific station
	double tmpBin = 0.;

	TH1F *setCh; // Receive Charge from raw data
	TH1F *setPk; // Receive Peak from raw data
	TH1F *offSetCh = new TH1F (); // For Charge offset
	TH1F *offSetPk = new TH1F (); // For Peak offset

	unsigned int nentry = 0;
	
	TTree *treeCh = new TTree("Charge","");
	TGraphErrors *fllHCh = new TGraphErrors ();
	treeCh->Branch("fllHCh", "TGraphErrors", &fllHCh);
	treeCh->Branch("charge", "TH1F", &charge);
	treeCh->Branch("offSetCh", "TH1F", &offSetCh);
	treeCh->Branch("offSetPk", "TH1F", &offSetPk);

	TTree *treePk = new TTree("Peak","");
	TGraphErrors *fllHPk = new TGraphErrors ();
	treePk->Branch("fllHPk", "TGraphErrors", &fllHPk);
	TTree *forEntries = new TTree("forEntries","");
	forEntries->Branch("nentry",&nentry,"nentry/I");

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
    IoSdEvent event(pos);

    if ( event.Id == previusEvent )
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

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i){
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter)
        if (event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;     
      
      if ( event.Stations[i].IsUUB ) {
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1 << " " << sameUtc
          << endl;
 
        IoSdEvent event(pos);
				
				blpmth.clear();

        for (unsigned int ii=0; ii<event.Stations[i].UFadc->NSample; ii++) 
          blpmth.push_back( 0 );
        
        if (event.Stations[i].Error==256) { //0+256
          for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++)
						blpmth[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
          for ( unsigned int id=0; id<totSt; id++ )
            if ( stationsBins[id] == event.Stations[i].Id ) {
							chok = 0.;
							pkok = 0.;
              meanf = getmrms.getMean(blpmth, nblbins, true);
              meanl = getmrms.getMean(blpmth, nblbins, false);
							rmsf = getmrms.getRms(blpmth, meanf, nblbins, true);
							rdHist.fitPkOk = false;
							rdHist.fitChOk = false;
							if ( fabs(meanl-meanf) < 2*rmsf ) {
								setCh = event.Stations[i].HCharge(pmtId-1);
								cerr << event.Stations[i].Histo->Offset[pmtId-1+6] << endl;
								if ( totSt==1 ) // To check the fit hist by hist
									for( int b=0; b<setCh->GetXaxis()->GetNbins(); b++ ) {
										tmpBin = setCh->GetBinCenter(b);
										stckHc.Fill( tmpBin, setCh->GetBinContent(b) );
									}
								nentry++;
								rdHist.getFullFit( *setCh, true, 0.1, 30, base ); // 10% of EMpeak, 30 bins forward from this last.
								if ( rdHist.fitChOk ) {
									stckVch[id] += rdHist.vemPos;
									stckEvt[id][0]++;
									chok = rdHist.vemPos;
									if ( totSt==1 ) { // To check the fit hist by hist	
 										fllHCh = rdHist.getFitGraph();
										charge = (TH1F*)setCh->Clone();
										treeCh->Fill();
									}
								}
								setCh->Reset();
								setPk = event.Stations[i].HPeak(pmtId-1);
								rdHist.getFullFit( *setPk, false, 0.1, 5, base ); // 10% of EMpeak, 5 bins forward from this last.
								if ( totSt==1 ) // To check the fit hist by hist
									for( int b=0; b<setPk->GetXaxis()->GetNbins(); b++ )
										stckHp.Fill( setPk->GetBinCenter(b)-setPk->GetBinCenter(0), setPk->GetBinContent(b) );
								if ( rdHist.fitPkOk ) { 
									stckVpk[id] += rdHist.vemPos;
									stckEvt[id][1]++;
									pkok = rdHist.vemPos;
									if ( totSt==1 ) { // To check the fit hist by hist
										fllHPk = rdHist.getFitGraph();
										treePk->Fill();
									}
								}
								if ( rdHist.fitChOk ) {
									stckAP[id] += chok/pkok;
									stckEvt[id][2]++;
									apDist.Fill( id, nrEventsRead, ubUub*chok/pkok ); 
									if ( totSt==1 )
										sglSt.Fill( nrEventsRead, ubUub*chok/pkok );
								}
								setPk->Reset();
							}
              break;
            }
        }
      }
    }    
  }

	forEntries->Fill();
  hfile.Write();
  hfile.Close();
	return 0;
}
