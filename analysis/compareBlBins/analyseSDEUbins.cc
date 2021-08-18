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


using namespace std;

double getmean( vector<int> *arr, unsigned nb, bool fok ){
  double mean = 0.;
	int lb = arr->size() - 1;
    for (unsigned int i=0; i<nb; i++){
      if ( fok )
        mean += (*arr)[i];
      else
        mean += (*arr)[lb-i];
    }
  return mean/nb;
}

double getrms( vector<int> *arr, double meanarr, unsigned int nb, bool fok ){
  double rms = 0.;
	int lb = arr->size() - 1;
  for (unsigned i=0; i<nb; i++){
    if ( fok )
      rms += ( (*arr)[i] - meanarr )*( (*arr)[i] - meanarr );
    else
      rms += ( (*arr)[lb-i] - meanarr )*( (*arr)[lb-i] - meanarr );
  }
  return sqrt(rms/nb);
}

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
	bool especSt = false;
	TString especStname;
	if ( stationsIds.size()==1 ){
		cerr << "Analysis on " << stationsIds[0] << " station." << endl;
		especSt = true;
		especStname = to_string( stationsIds[0] );
	}
  
  TString pmtname = whichpmt;
  unsigned int pmtId = atoi( pmtname );
  if ( pmtId > 0 && pmtId < 4 ){
     pmtname = "PMT"+to_string( pmtId );
  }
  else if ( pmtId == 4 )
    pmtname = "SPMT";
  else if ( pmtId == 5 )
    pmtname = "PMTSSD";
  else{
    cerr << "==================================================" << endl;
    cerr << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
    cerr << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " 
      << "4 For SPMT; " << "5 For PMTSSD" << endl;
    cerr << "==================================================" << endl;
    exit(0);
  }
  cerr << "You have selected " << pmtname << endl;

	TString outputName;
	if ( !especSt )
		outputName = "bl100bins"+pmtname+".root";
	else
		outputName = "bl100bins"+especStname+".root";
		
	TFile hfile (outputName,"RECREATE","");

  int totSt = stationsIds.size();

  vector < int > stckEvt;
  stckEvt.resize(totSt);
  vector < vector < float > > stckMean;
  stckMean.resize(totSt);
  vector < vector < float > > stckRMS;
  stckRMS.resize(totSt);

  sort(stationsIds.begin(), stationsIds.end());
  Double_t stationsBins[totSt];
  for ( int i=0; i<totSt; i++){
    stationsBins[i] = stationsIds[i];
    stckEvt[i] = 0;
    for ( int j=0; j<20; j++ ){
      stckMean[i].push_back(0);
      stckRMS[i].push_back(0.);
    }
    cout << i+1 << " -> " << stationsIds[i] << endl;
  }
  
  unsigned int nblbins = 100;
  unsigned int nday = 0;
  unsigned int totDays = 90; // From 1st December, 2020 to 28th February, 2021
  unsigned int cday = 1606867200; //2nd December, 2020;
  unsigned int dday = 86400; 
  double tmpMean = 0;

	double tmpMeanStf = 0;
	double tmpMeanStl = 0;
	double tmpMeanStf2 = 0;
	double tmpMeanStl2 = 0;
	double tmpMeanStf3 = 0;
	double tmpMeanStl3 = 0;
	TH1F stdiffDist("stdiffDist", "Distribution of The Difference of The Mean for the Station "+especStname+"'s PMT1 HG",31, -15, 15);
	TH1F stdiffDist2("stdiffDist2", "Distribution of The Difference of The Mean for the Station "+especStname+"'s PMT2 HG",31, -15, 15);
	TH1F stdiffDist3("stdiffDist3", "Distribution of The Difference of The Mean for the Station "+especStname+"'s PMT3 HG",31, -15, 15);
  
  // For Low Gain
  TH2F pmtlmeanf("pmtlmeanf", "Mean for First 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmtlmeanl("pmtlmeanl", "Mean for Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmtlrmsf("pmtlrmsf", "RMS for First 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmtlrmsl("pmtlrmsl", "RMS for Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmtldiffmean("pmtldiffmean", "Difference of the mean of First and Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmtldiffrms("pmtldiffrms", "Difference of the RMS of First and Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

  // For High Gain
  TH2F pmthmeanf("pmthmeanf", "Mean for First 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmthmeanl("pmthmeanl", "Mean for Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmthrmsf("pmthrmsf", "RMS for First 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmthrmsl("pmthrmsl", "RMS for Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmthdiffmean("pmthdiffmean", "Difference of the Mean for First and Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmthdiffrms("pmthdiffrms", "Difference of the RMS for First and Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);



  vector < int > *blpmth = new vector < int >;
  vector < int > *blpmth2 = new vector < int >;
  vector < int > *blpmth3 = new vector < int >;
  vector < int > *blpmtl = new vector < int >;

  unsigned int previusEvent = 0;
  unsigned int sameUtc = 0;

  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;
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
    
    if ( !especSt && sameUtc > cday ){
      for ( int id=0; id<totSt; id++){
        if ( stckEvt[id] > 0 ){
          pmthmeanf.Fill( nday , id, stckMean[id][0]/stckEvt[id] );
          pmthmeanl.Fill( nday , id, stckMean[id][1]/stckEvt[id] );
          pmthrmsf.Fill( nday , id, stckRMS[id][0]/stckEvt[id] );
          pmthrmsl.Fill( nday , id, stckRMS[id][1]/stckEvt[id] );
          pmthdiffmean.Fill( nday, id, stckMean[id][0]/stckEvt[id] - stckMean[id][1]/stckEvt[id] );
          pmthdiffrms.Fill( nday, id, stckRMS[id][0]/stckEvt[id] - stckRMS[id][1]/stckEvt[id] );  

          pmtlmeanf.Fill( nday , id, stckMean[id][2]/stckEvt[id] );
          pmtlmeanl.Fill( nday , id, stckMean[id][3]/stckEvt[id] );
          pmtlrmsf.Fill( nday , id, stckRMS[id][2]/stckEvt[id] );
          pmtlrmsl.Fill( nday , id, stckRMS[id][3]/stckEvt[id] );
          pmtldiffmean.Fill( nday, id, stckMean[id][2]/stckEvt[id] - stckMean[id][3]/stckEvt[id] );
          pmtldiffrms.Fill( nday, id, stckRMS[id][2]/stckEvt[id] - stckRMS[id][3]/stckEvt[id] );

          for ( int k=0; k<4; k++ ){
            stckMean[id][k] = 0;
            stckRMS[id][k] = 0;
            stckEvt[id] = 0;
          }
        }
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
      
      if ( event.Stations[i].IsUUB ){
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1 << " " << sameUtc
          << endl;
 
        IoSdEvent event(pos);
				blpmth->resize(event.Stations[i].UFadc->NSample);
				blpmth2->resize(event.Stations[i].UFadc->NSample);
				blpmth3->resize(event.Stations[i].UFadc->NSample);
				blpmtl->resize(event.Stations[i].UFadc->NSample);
			
				for (unsigned int i=0; i<blpmth->size(); i++){
					(*blpmth)[i] = 0;
					(*blpmth2)[i] = 0;
					(*blpmth3)[i] = 0;
					(*blpmtl)[i] = 0;
				}

        
        if ( event.Stations[i].Error==256 ) { //0+256
          for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++){
            (*blpmth)[k] = event.Stations[i].UFadc->GetValue(pmtId-1,0,k);
            (*blpmtl)[k] = event.Stations[i].UFadc->GetValue(pmtId-1,1,k);

						if ( especSt ){
							(*blpmth2)[k] = event.Stations[i].UFadc->GetValue(1,0,k);
							(*blpmth3)[k] = event.Stations[i].UFadc->GetValue(2,0,k);
						}
          }

					if ( !especSt ){
						for ( int id=0; id<totSt; id++ )
							if ( stationsBins[id] == event.Stations[i].Id ){
								stckEvt[id]++;
  	            tmpMean = getmean(blpmth, nblbins, true);
    	          stckMean[id][0] += tmpMean;
      	        stckRMS[id][0] += getrms(blpmth, tmpMean, nblbins, true);
        	      tmpMean = getmean(blpmth, nblbins, false);
          	    stckMean[id][1] += tmpMean;
            	  stckRMS[id][1] += getrms(blpmth, tmpMean, nblbins, false);
								
								tmpMean = getmean(blpmtl, nblbins, true);
    	          stckMean[id][2] += tmpMean;
      	        stckRMS[id][2] += getrms(blpmtl, tmpMean, nblbins, true);
        	      tmpMean = getmean(blpmtl, nblbins, false);
          	    stckMean[id][3] += tmpMean;
            	  stckRMS[id][3] += getrms(blpmtl, tmpMean, nblbins, false);
								break;
							}
					}
					else{
						tmpMeanStf = getmean(blpmth, nblbins, true);
            tmpMeanStl = getmean(blpmth, nblbins, false);
            stdiffDist.Fill( tmpMeanStf-tmpMeanStl );
            tmpMeanStf2 = getmean(blpmth2, nblbins, true);
            tmpMeanStl2 = getmean(blpmth2, nblbins, false);
            stdiffDist2.Fill( tmpMeanStf2-tmpMeanStl2 );
            tmpMeanStf3 = getmean(blpmth3, nblbins, true);
            tmpMeanStl3 = getmean(blpmth3, nblbins, false);
            stdiffDist3.Fill( tmpMeanStf3-tmpMeanStl3 );
					}
				}
      }
    }    
  }

	cerr << nday << endl;

  hfile.Write();
  hfile.Close();

  return 0;
}
