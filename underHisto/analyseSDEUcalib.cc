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


double getmean( vector<int> *arr, unsigned int nb, bool lok ){
  double mean = 0.;
  int lb = arr->size() - 1;
    for  (unsigned int i=0; i<nb; i++){
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
  for (unsigned int i=0; i<nb; i++){
    if ( lok == 0 )
      rms += ((*arr)[i] - meanarr)*((*arr)[i] - meanarr);
    else
      rms += ((*arr)[lb-i] - meanarr)*((*arr)[lb-i] - meanarr);
  }
  return sqrt(rms/nb);
}

// ========================== 
// ******** The MAIN ********
// ==========================
int main (int argc, char *argv[]) {
   if ( argc < 4 ) {
		 cout << endl
         << "Usage: " << argv[0] << " <stationsFile>  <PMT>  <Month> <files>" << endl
         << "  <stationsFile>: file with a list of stations" << endl
         << "  <PMT>: ID for the PMT you want to analyse" << endl
         << "  <Month>: Month in which you want to analyse" << endl
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

	if ( totSt==1 )
		pmtname += "St"+to_string( stationsIds[0] );
 
  string doMonth = string(whichmonth);
  pmtname +=  "Mth" + doMonth;
  TFile hfile("uubAoPtime"+pmtname+"chpk.root","RECREATE","");

	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

  TGraphErrors *pkHistFit = new TGraphErrors();
  double pkChi2 = 0.;
  double peak = 0.;

  TGraphErrors *chHistFit = new TGraphErrors();
  double chChi2 = 0.;
  double charge = 0.; 

  unsigned int nrEventsRead = 0;
  unsigned int nrEvents = 0;
  bool found = false;
  unsigned int nblbins = 100;
	double blCorrHbase = 0.;
  double meanf = 0.;
  double meanl = 0.;
  double rmsf = 0.;

  TH1F *tmp = new TH1F();
  TString tmpName;
  fitpeak fitPk;
  fitcharge fitCh;
  vector < int > *blpmth  = new vector < int >;

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) 
  {
    nrEventsRead++;
    if (nrEventsRead%1000 == 0) 
    {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    IoSdEvent event(pos);

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i)
    {
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter)
        if (event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;
      
      if ( event.Stations[i].IsUUB )
      {
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1
          << endl;
        if (event.Stations[i].Error==256)
        {
          blpmth->clear();
          blpmth->resize(event.Stations[i].UFadc->NSample);
          for ( unsigned int k=0;k<event.Stations[i].UFadc->NSample;k++ )
            (*blpmth)[k] = ( event.Stations[i].UFadc->GetValue(pmtId-1,0,k) );

          meanf = getmean(blpmth, nblbins, false);
          meanl = getmean(blpmth, nblbins, true);
          rmsf = getrms(blpmth, meanf, nblbins, false);
          if (fabs(meanl-meanf) < 2*rmsf ) // Reading event with stable baseline
          {
            // ================
            // *** Baseline ***
            blCorrHbase = event.Stations[i].HBase(pmtId-1)->GetMean(); // Extracting calib-baseline
            tmpName.Form("%d%d", event.UTCTime, nrEventsRead-1);

            recePk = event.Stations[i].HPeak(pmtId-1); // Receiving Peak histogram
            fitPk.getCrrSmooth(*recePk, blCorrHbase, tmpName+"Hbpk"); // Correcting for calib-baseline
            tmp = fitPk.getPkCorrSmooth();
            //for ( int k = 0; k<tmp->GetXaxis()->GetNbins(); k++ )
              //cout << k << " " << tmp->GetBinCenter(k) << " " << tmp->GetBinContent(k) << endl;
            fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]); // Fitting

            pkHistFit = fitPk.getFitGraphPk();
            pkChi2 = fitPk.chisPeak;
            peak = fitPk.vemPosPk;
            cerr << event.Id << endl;        
            
            receCh = event.Stations[i].HCharge(pmtId-1);
            fitCh.getChCrrSmooth(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
            tmp = fitCh.getChCorrSmooth();
            for ( int k = 0; k<tmp->GetXaxis()->GetNbins(); k++ )
              cout << k << " " << tmp->GetBinCenter(k) << " " << tmp->GetBinContent(k) << endl;
            fitCh.getFitCh(*tmp, 0.2, 30, event.Stations[i].Calib->VemCharge[pmtId-1]);             
          
            chHistFit = fitCh.getFitGraphCh();
            chChi2 = fitCh.chisCharge;
            charge =  fitCh.vemPosCh;
 
            exit(0);           

            
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
