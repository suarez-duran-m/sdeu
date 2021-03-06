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
  TFile hfile("uubAoP"+pmtname+".root","RECREATE","");

	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

  TGraphErrors *pkHistFit = new TGraphErrors();
  double pkChi2 = 0.;
  double pkNdf = 0.;
  double pkProb = 0.;
  double peak = 0.;

  TGraphErrors *chHistFit = new TGraphErrors();
  double chChi2 = 0.;
  double chNdf = 0.;
  double chProb = 0.;
  double charge = 0.; 

  unsigned int evtIdPk = 0; //Storing event Id
  unsigned int evtTimePk = 0; //Storing day-Unixtime
  unsigned int evtIdCh = 0; //Storing event Id
  unsigned int evtTimeCh = 0; //Storing day-Unixtime
  unsigned int previusEvent = 0; // Avoiding read the same event
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

  // ======================================
  // *** *** *** Creating Trees *** *** *** 
  TTree *treePeak = new TTree("PeakData","");
  TTree *treeCharge = new TTree("ChargeData","");

  //treePeak->Branch("","",&);
  treePeak->Branch("peakVal",&peak,"peak/D");
  treePeak->Branch("chi2",&pkChi2,"pkChi2/D");
  treePeak->Branch("ndf",&pkNdf,"pkNdf/D");
  treePeak->Branch("prob",&pkProb,"pkProb/D");
  treePeak->Branch("eventId",&evtIdPk,"evtIdPk/I");
  treePeak->Branch("timeEvnt",&evtTimePk,"evtTimePk/I");

  treeCharge->Branch("chargeVal",&charge,"charge/D");
  treeCharge->Branch("chi2",&chChi2,"chChi2/D");
  treeCharge->Branch("ndf",&chNdf,"chNdf/D");
  treeCharge->Branch("prob",&chProb,"chProb/D");
  treeCharge->Branch("eventId",&evtIdCh,"evtIdCh/I");
  treeCharge->Branch("timeEvnt",&evtTimeCh,"evtTimeCh/I");

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
    if ( event.Id == previusEvent )
      continue;

    previusEvent = event.Id;

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
            fitPk.getCrr(*recePk, blCorrHbase, tmpName+"Hbpk"); // Correcting for calib-baseline
            tmp = fitPk.getPkCorr(); // Receiving corrected histogram
            fitPk.getFitPk(*tmp); // Fitting

            pkHistFit = fitPk.getFitGraphPk();
            pkChi2 = fitPk.chisPeak;
            pkNdf = fitPk.ndfPeak;
            pkProb = fitPk.probPeak;
            peak = fitPk.vemPosPk;
            /*
            if ( peak==0 )
            {
              cout << "MSD " << " " << event.Id << " " << pkChi2 << " " << pkNdf << " " << pkChi2/pkNdf << endl;
              for ( int kk=0; kk<tmp->GetXaxis()->GetNbins(); kk++ )
                cout << kk << " " << tmp->GetBinCenter(kk) << " " << tmp->GetBinContent(kk) << endl;
              exit(0);
            }
            */
            
            if ( pkChi2/pkNdf > 4 ) //5.0e+08 ) //( pkChi2/pkNdf > 1.3 && pkChi2/pkNdf < 1.7 )
            {
              int tmpkk = 0;
              for ( int kk=0; kk<10; kk++ )
                tmpkk += tmp->GetBinContent(kk);
              if ( tmpkk > 0 )
              {
                cout << "MSD " << " " << event.Id << " " << pkChi2 << " " << pkNdf << " " << pkChi2/pkNdf << endl;
                for ( int kk=0; kk<tmp->GetXaxis()->GetNbins(); kk++ )
                  cout << kk << " " << tmp->GetBinCenter(kk) << " " << tmp->GetBinContent(kk) << endl;
                exit(0);
              }
            }
            
            receCh = event.Stations[i].HCharge(pmtId-1);
            fitCh.getChCrr(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
            tmp = fitCh.getChCrr();           
            fitCh.getFitCh(*tmp);
          
            chHistFit = fitCh.getFitGraphCh();
            chChi2 = fitCh.chisCharge;
            chNdf = fitCh.ndfCharge;
            chProb = fitCh.probCharge;
            charge =  fitCh.vemPosCh;
            /*
            if ( chChi2/chNdf > 4 ) //5.0e+08 ) //( pkChi2/pkNdf > 1.3 && pkChi2/pkNdf < 1.7 )
            {
              cout << "MSD " << " " << event.Id << " " << chChi2 << " " << chNdf << " " << chChi2/pkNdf << endl;
              for ( int kk=0; kk<tmp->GetXaxis()->GetNbins(); kk++ )
                cout << kk << " " << tmp->GetBinCenter(kk) << " " << tmp->GetBinContent(kk) << endl;
              exit(0);
            }
            */
            
            evtIdPk = event.Id;
            evtIdCh = event.Id;
            evtTimePk = event.utctime();
            evtTimeCh = event.utctime();

            treePeak->Fill();
            treeCharge->Fill();

            //exit(0);
            
		  			break; // Apply if the is running for a single station.
          }
        }
      }
    }
  }

  hfile.Write();
  hfile.Close();
  
	return 0;
}
