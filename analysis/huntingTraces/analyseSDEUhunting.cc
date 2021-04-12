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

double getmean( vector<int> *arr, unsigned int nb, bool lok ){
  double mean = 0.;
  int lb = arr->size() - 1;
    for  (int i=0; i<nb; i++){
      if ( !lok )
        mean += (*arr)[i];
      else
        mean += (*arr)[lb-i];
    }
  return mean/nb;
}

double getrms( vector<int> *arr, double meanarr, unsigned int nb, bool lok ){
  double rms = 0.;
  int lb = arr->size() - 1;
  for (int i=0; i<nb; i++){
    if ( lok == 0 )
      rms += ((*arr)[i] - meanarr)*((*arr)[i] - meanarr);
    else
      rms += ((*arr)[lb-i] - meanarr)*((*arr)[lb-i] - meanarr);
  }
  return sqrt(rms/nb);
}

double getlowfr( int arr[], double mean, double rms, unsigned int nb, bool lok ){
  //int lowfr = 0;
  int sum = 0;
  double cond = 100.*(mean-6.*rms);
  for ( int i=0; i<nb; i++ )
    if ( !lok )
      if ( (arr[i] - mean) < 0 )
        sum += arr[i];
    else
      if ( (arr[2047-i] - mean) < 0 )
        sum += arr[2047-i];
  if ( sum < cond )
    return 1;
  else
    return 0;
}
// ========================== 
// ******** The MAIN ********
// ==========================
int main (int argc, char *argv[]) {
   if ( argc < 4 ) {
     cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <output>  <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <output>: output file with whatever you want inside" << endl
         << "  <files>: IoSd or IoAuger files to be read" << endl;
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
  cout << "You have selected " << pmtname << endl;

  unsigned int nblbins = 100;
  double meanf = 0.;
  double meanl = 0.;
  double rmsf = 0.;
  double rmsl = 0.;
  double diffm = 0.;
  double diffr = 0.;

  TFile hfile("zooTraces100bins"+pmtname+".root","RECREATE","");

  TH2F badTrac ("badTrac","Traces not Ok "+pmtname+"HG", totalNrEvents, 0, totalNrEvents, 25*nblbins, 0, 25*nblbins);
  TH2F gooTrac ("gooTrac","Traces not Ok "+pmtname+"HG", totalNrEvents, 0, totalNrEvents, 25*nblbins, 0, 25*nblbins);
  TH1F badTrEv ("badTrEv","Traces not Ok "+pmtname+"HG", totalNrEvents, 0, totalNrEvents);

  int totSt = stationsIds.size();
  int lowfr = 0;

  vector < int > stckEvt;
  stckEvt.resize(totSt);
  vector < vector < int > > stckOk;
  stckOk.resize(totSt);
  vector < vector < double > > stckUn;
  stckUn.resize(totSt);
  vector < vector < int > > stckFr;
  stckFr.resize(totSt);
  vector < vector < int > > stckDrms;
  stckDrms.resize(totSt);

  sort(stationsIds.begin(), stationsIds.end());
  Double_t stationsBins[totSt];
  for ( int i=0; i<totSt; i++){
    stationsBins[i] = stationsIds[i];
    stckEvt[i] = 0;
    for ( int j=0; j<4; j++ ){
      stckOk[i].push_back(0);
      stckUn[i].push_back(0.);
      stckFr[i].push_back(0.);
      stckDrms[i].push_back(0.);
    }
    cout << i << " -> " << stationsIds[i] << endl;
  }
  
  // For Low Gain
  TH1F pmtlok("pmtlok", "Fraction of Events Ok for "+pmtname+" LG", totSt, 0, totSt);
  TH1F pmtlun("pmtlun", "Fraction of Events Undershoot for "+pmtname+" LG", totSt, 0, totSt);
  TH1F pmtlfr("pmtlfr", "Fraction of Events LowFre for "+pmtname+" LG", totSt, 0, totSt);
  /*
  TH2F pmtlmeanl("pmtlmeanl", "Mean for Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmtlrmsf("pmtlrmsf", "RMS for First 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmtlrmsl("pmtlrmsl", "RMS for Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmtldiffmean("pmtldiffmean", "Difference of the mean of First and Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmtldiffrms("pmtldiffrms", "Difference of the RMS of First and Last 100 bins of Baseline "+pmtname+" LG", totDays, 0, totDays, totSt, 0, totSt);
  */

  // For High Gain
  TH1F pmthok("pmthok", "Fraction of Events Ok for "+pmtname+" HG", totSt, 0, totSt);
  TH1F pmthun("pmthun", "Fraction of Events Undershoot for "+pmtname+" HG", totSt, 0, totSt);
  TH1F pmthfr("pmthfr", "Fraction of Events LowFre for "+pmtname+" HG", totSt, 0, totSt);

  /*
  TH2F pmthmeanl("pmthmeanl", "Mean for Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmthrmsf("pmthrmsf", "RMS for First 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmthrmsl("pmthrmsl", "RMS for Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);

  TH2F pmthdiffmean("pmthdiffmean", "Difference of the Mean for First and Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  TH2F pmthdiffrms("pmthdiffrms", "Difference of the RMS for First and Last 100 bins of Baseline "+pmtname+" HG", totDays, 0, totDays, totSt, 0, totSt);
  */

  vector < int > *blpmth = new vector < int >;
  vector < int > *blpmtl = new vector < int >;

  unsigned int previusEvent = 0;
  int sameUtc = 0;

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
        for (int i=0; i<blpmth->size(); i++){
          (*blpmth)[i] = 0;
          (*blpmtl)[i] = 0;
        }

        if (event.Stations[i].Error==256) { //0+256
          for (unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++){
            (*blpmth)[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
            (*blpmtl)[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );
          }
          for ( int id=0; id<totSt; id++)
            if ( stationsBins[id] == event.Stations[i].Id ){
              stckEvt[id]++;
              meanf = getmean(blpmth, nblbins, false);
              meanl = getmean(blpmth, nblbins, true);
              rmsf = getrms(blpmth, meanf, nblbins, false);
              rmsl = getrms(blpmth, meanl, nblbins, true);
              //lowfr = getlowfr(blpmth, meanl, rmsl, nblbins, true);
              diffm = meanf - meanl;
              diffr = rmsf - rmsl;
              if (fabs(meanl-meanf) < 2*rmsf){
                stckOk[id][0]++;
                if( id==8 )
                  for (int k=0; k<blpmth->size(); k++)
                    gooTrac.Fill( nrEventsRead, k, (*blpmth)[k] );
              }
              else if( id==9 ){
                //cerr << previusEvent << " " << nrEventsRead << endl;
                int tmp = blpmth->size();
                for (int k=0; k<tmp; k++){
                  badTrac.Fill( nrEventsRead, k, (*blpmth)[k] );
                  //badTrac.Fill( nrEventsRead, k+nblbins, (*blpmth)[tmp-k] );
                }
                badTrEv.Fill( nrEventsRead, event.Id );                
              }
              break;
            }

              /*
              meanf = getmean(blpmtl, nblbins, 0);
              meanl = getmean(blpmtl, nblbins, 1);
              rmsf = getrms(blpmtl, meanf, nblbins, 1);
              rmsl = getrms(blpmtl, meanl, nblbins, 1);
              lowfr = getlowfr(blpmtl, meanf, rmsf, nblbins, 1);
              diffm = meanf - meanl;
              diffr = rmsf - rmsl;
              if ( rmsl > 0.5*rmsf && rmsl < 1.5*rmsf ){
                if ( meanl > (meanf - 2*rmsf) && meanl < (meanf+2*rmsf) ){
                  stckOk[id][1]++;
                }
              }
              else if ( meanl < (meanf+2.5*rmsf) ){
                stckUn[id][1]++;
                cout << "low: " << id << " " << meanl << endl;
              }
              else if ( lowfr ){
                cout << "lowlowfr: " << id << " " << meanl << endl;
                stckFr[id][1]++;
              }
              */
        }
      }
    }    
  }

  for ( int id=0; id<totSt; id++ ){
    pmthok.Fill( id, (1.*stckOk[id][0])/(1.*stckEvt[id]) );
    pmtlok.Fill( id, (1.*stckOk[id][1])/(1.*stckEvt[id]) );
    pmthun.Fill( id, (1.*stckUn[id][0])/(1.*stckEvt[id]) );
    pmtlun.Fill( id, (1.*stckUn[id][1])/(1.*stckEvt[id]) );
    pmthfr.Fill( id, (1.*stckFr[id][0])/(1.*stckEvt[id]) );
    pmtlfr.Fill( id, (1.*stckFr[id][1])/(1.*stckEvt[id]) );
  }

  hfile.Write();
  hfile.Close();

   return 0;
}
