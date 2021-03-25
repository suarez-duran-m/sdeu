#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TCanvas.h>

#include "stats.h"

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

  const unsigned int nblbins = 100;
	stats getmrms;
  double meanf = 0.;
  double meanl = 0.;
  double rmsf = 0.;
  double rmsl = 0.;

  TFile hfile("zooTraces100bins"+pmtname+".root","RECREATE","");

	
	vector < int > lightning; 
	lightning.push_back( 1608768000 ); // Dec. 24th, 2020
	lightning.push_back( 1609113600 ); // Dec. 28th, 2020
	lightning.push_back( 1609459200 ); // Jan. 1st, 2021
	lightning.push_back( 1609632000 ); // Jan. 3th, 2021
	lightning.push_back( 1609718400 ); // Jan. 4th, 2021
	lightning.push_back( 1609891200 ); // Jan. 6th, 2021
	lightning.push_back( 1613001600 ); // Feb. 11th, 2021
	lightning.push_back( 1613433600 ); // Feb. 16th, 2021

  unsigned int totSt = stationsIds.size();

	vector < vector < int > > stckEvt;
  stckEvt.resize(totSt);
  vector < vector < double > > stckdRMS;
  stckdRMS.resize(totSt);
	
  sort(stationsIds.begin(), stationsIds.end());
  Double_t stationsBins[totSt];
  for ( unsigned int i=0; i<totSt; i++){
    stationsBins[i] = stationsIds[i];
		for ( unsigned int j=0; j<8; j++ ) {
			if ( j<4 )
				stckEvt[i].push_back(0.);
			stckdRMS[i].push_back(0.); // 0 for h and okf; 1 for h and okl.; ... 4 for l and okf.. 
		}
    cout << i << " -> " << stationsIds[i] << endl;
  }

  vector < int > *blpmth = new vector < int >;
  vector < int > *blpmtl = new vector < int >;
	vector < int > *blStpmt2 = new vector < int >;
	vector < int > *blStpmt3 = new vector < int >;

  unsigned int previusEvent = 0;
  unsigned int sameUtc = 0;
	unsigned int nday = 0;
	unsigned int cday = 1606867200;
	unsigned int dday = 86400;
	unsigned int totDays = 90; // From 1st December, 2020 to 28th February, 2021

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;

  TH2D okDiffRmsh ("okDiffRmsh","Difference of the RMS for the first and last 100 bins, selected traces, for "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2D discDiffRmsh ("discDiffRmsh","Difference of the RMS for the first and last 100 bins, discarted traces, for "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);

  TH2D okDiffRmsl ("okDiffRmsl","Difference of the RMS for the first and last 100 bins, selected traces, for "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2D discDiffRmsl ("discDiffRmsl","Difference of the RMS for the first and last 100 bins, discarted traces, for "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

	TH2I stTraces ("stTraces","Traces, for "+pmtname+" HG", totalNrEvents, 0, totalNrEvents, 2048, 0, 2048);
	TH2I stTracesNo ("stTracesNo","Traces, for "+pmtname+" HG", totalNrEvents, 0, totalNrEvents, 2048, 0, 2048);

	TH1D stRms1 ("stRms1", "PMT1", 210, -10, 10); 
	TH1D stRms2 ("stRms2", "PMT2", 210, -10, 10);
	TH1D stRms3 ("stRms3", "PMT3", 210, -10, 10);

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0){
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
    IoSdEvent event(pos);

    if ( event.Id == previusEvent ) //&& sameUtc == event.utctime() )
      continue;

    previusEvent = event.Id;
    sameUtc = event.utctime();

		if ( event.utctime() >= lightning[0] && event.utctime() <= lightning[1] )
			continue;
		else if ( event.utctime() >= lightning[2] && event.utctime() <= lightning[3] )
			continue;
		else if ( event.utctime() >= lightning[4] && event.utctime() <= lightning[5] )
			continue;
		else if( event.utctime() >= lightning[6] && event.utctime() <= lightning[7] )
			continue;

			//cerr << nday << " " << event.utcdate() << " " << nday*dday + 1606867200 << " " << cday << " " << (event.utctime() - 1606867200)/dday +1 << endl;

    if ( sameUtc > cday ){
      for ( unsigned int id=0; id<totSt; id++) {
        if ( stckEvt[id][0] > 0 )
					okDiffRmsh.Fill( nday , id, (stckdRMS[id][0] - stckdRMS[id][1])/stckEvt[id][0] );
				if ( stckEvt[id][1] > 0 )
					discDiffRmsh.Fill( nday , id, (stckdRMS[id][2] - stckdRMS[id][3])/stckEvt[id][1] );
				if ( stckEvt[id][2] > 0 )
					okDiffRmsl.Fill( nday , id, (stckdRMS[id][4] - stckdRMS[id][5])/stckEvt[id][2] );
				if ( stckEvt[id][3] > 0 )
					discDiffRmsl.Fill( nday , id, (stckdRMS[id][6] - stckdRMS[id][7])/stckEvt[id][3] );
				for ( int k=0; k<8; k++ ){
					if ( k<4 )
						stckEvt[id][k] = 0;
					stckdRMS[id][k] = 0;
				}
			}
			cday += dday;
			nday++;
		}
		if ((event.utctime() - 1606867200)/dday +1 > nday )
			continue;



    for (unsigned int i = 0 ; i < event.Stations.size(); ++i){
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter)
				if (event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;
      if ( event.Stations[i].IsUUB ){
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1 << " " << sameUtc
          << endl;
 
        IoSdEvent event(pos);
				
        blpmth->resize(event.Stations[i].UFadc->NSample);
        blpmtl->resize(event.Stations[i].UFadc->NSample);
				blStpmt2->resize(event.Stations[i].UFadc->NSample);
				blStpmt3->resize(event.Stations[i].UFadc->NSample);
        for (unsigned int i=0; i<blpmth->size(); i++) {
          (*blpmth)[i] = 0;
          (*blpmtl)[i] = 0;
          (*blStpmt2)[i] = 0;
          (*blStpmt3)[i] = 0;
        }

        if (event.Stations[i].Error==256) { //0+256
          for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++){
            (*blpmth)[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
            (*blpmtl)[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,1,k) );
						if ( event.Stations[i].Id==1851 ) {
							(*blStpmt2)[k] = ( event.Stations[i].UFadc->GetValue(1,0,k) );
							(*blStpmt3)[k] = ( event.Stations[i].UFadc->GetValue(2,0,k) );
						}
						//if ( event.Id==61252751 && event.Stations[i].Id==1818 )
							//sglEvt.Fill( k, event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
          }
          for ( unsigned int id=0; id<totSt; id++ )
            if ( stationsBins[id] == event.Stations[i].Id ) {
              meanf = getmrms.getMean( blpmth, nblbins, true );
              meanl = getmrms.getMean( blpmth, nblbins, false );
							rmsf = getmrms.getRms( blpmth, meanf, nblbins, true );
							rmsl = getmrms.getRms( blpmth, meanl, nblbins, false );
							if ( fabs(meanl-meanf) < 2*rmsf ){
								stckEvt[id][0]++;
								stckdRMS[id][0] += rmsf;
								stckdRMS[id][1] += rmsl;
								if ( event.Stations[i].Id==1851 )
									for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++)
										stTraces.Fill( nrEventsRead-1, k, event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
							}
							else{
								stckEvt[id][1]++;
								stckdRMS[id][2] += rmsf;
								stckdRMS[id][3] += rmsl;
								if ( event.Stations[i].Id==1851 )
									for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++)
										stTracesNo.Fill( nrEventsRead-1, k, event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
							}
							if ( event.Stations[i].Id==1851 ){
								stRms1.Fill( rmsl );
								//stRms1.Fill( rmsf - rmsl );
								meanf = getmrms.getMean( blStpmt2, nblbins, true );
								meanl = getmrms.getMean( blStpmt2, nblbins, false );
								rmsf = getmrms.getRms( blStpmt2, meanf, nblbins, true );
								rmsl = getmrms.getRms( blStpmt2, meanl, nblbins, false );
								stRms2.Fill( rmsl );
								//stRms2.Fill( rmsf - rmsl );
								meanf = getmrms.getMean( blStpmt3, nblbins, true );
								meanl = getmrms.getMean( blStpmt3, nblbins, false );
								rmsf = getmrms.getRms( blStpmt3, meanf, nblbins, true );
								rmsl = getmrms.getRms( blStpmt3, meanl, nblbins, false );
								stRms3.Fill( rmsl );
								//stRms3.Fill( rmsf - rmsl );
							}	
							meanf = getmrms.getMean( blpmtl, nblbins, true );
							meanl = getmrms.getMean( blpmtl, nblbins, false );
							rmsf = getmrms.getRms( blpmtl, meanf, nblbins, true );
							rmsl = getmrms.getRms( blpmtl, meanl, nblbins, false );
							if (fabs(meanl-meanf) < 2*rmsf ){
								stckEvt[id][2]++;
								stckdRMS[id][4] += rmsf;
								stckdRMS[id][5] += rmsl;
							}
							else{
								stckEvt[id][3]++;								
								stckdRMS[id][6] += rmsf;
								stckdRMS[id][7] += rmsf;
							}
              break;
            }
        }
      }
    }    
  }

  hfile.Write();
  hfile.Close();

   return 0;
}
