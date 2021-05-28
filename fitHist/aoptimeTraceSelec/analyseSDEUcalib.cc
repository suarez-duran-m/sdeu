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

  unsigned int totDays = 0;
	unsigned int cday = 0;

  if ( doMonth == "dec" )
  {
    totDays = 31;
    cday = 1606867200;
  }
  else if ( doMonth == "jan" )
  {
    totDays = 31;
    cday = 1609545600;
  }
  else if ( doMonth == "feb" )
  {
    totDays = 28;
    cday = 1612224000;
  }
  else if ( doMonth == "mar" )
  {
    totDays = 31;
    cday = 1614643200;
  }
  else if ( doMonth == "abr" )
  {
    totDays = 30;
    cday = 1617321600;
  }
  else
  {
    cerr << "Month no found" << endl;
    exit(0);
  }

  unsigned int previusEvent = 0;
  unsigned int sameUtc = 0;
	unsigned int nday = 0;
	unsigned int dday = 86400;
  unsigned int nrEvents = 0;
  unsigned int nrEventsRead = 0;

  vector < double > stckapHb; // It stores ap Hb per day
  vector < double > stckchHb; // It stores Ch Hb per day
  vector < double > stckpkHb; // It stores Pk Hb per day
	
	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

	TH1F *apHbase = new TH1F("apHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *pkHbase = new TH1F("pkHbase","",totDays, 0, totDays);//It stores VemP from HBase for St.
  TH1F *pkChis = new TH1F("pkChis", "", 500, 0, 50); // It stores Chis^2/NdF for Peak
	TH1F *chHbase = new TH1F("chHbase","",totDays, 0, totDays);//It stores VemQ from HBase for St.
  TH1F *chChis = new TH1F("chChis", "", 500, 0, 50); // It stores Chis^2/NdF for Charge
	TH1F *rmsHbase = new TH1F("rmsHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *rmspkHbase = new TH1F("rmspkHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *rmschHbase = new TH1F("rmschHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *timestmp = new TH1F("timestmp","",totDays, 0, totDays);//It stores A/P from Calib for St.

  TGraphErrors *pkHistFit = new TGraphErrors();
  double pkChi2 = 0.;
  TGraphErrors *chHistFit = new TGraphErrors();
  double chChi2 = 0.;

	unsigned int eventStat = 0; //It stores the number of events.
	unsigned int evtTime = 0; //It stores the day-Unixtime

  TTree *treeHist = new TTree("Histograms","");
	treeHist->Branch("eventStat", &eventStat, "eventStat/I");
	treeHist->Branch("evtTime",&evtTime,"evtTime/I");

  TTree *treeHistCh2 = new TTree("HistForChi2","");
  treeHistCh2->Branch("pkHistFit","TGraphErrors", &pkHistFit);
  treeHistCh2->Branch("pkChi2", &pkChi2, "pkChi2/D");
  treeHistCh2->Branch("chHistFit","TGraphErrors", &chHistFit);
  treeHistCh2->Branch("chChi2", &chChi2, "chChi2/D");


  unsigned int nblbins = 100;
	double blCorrHbase = 0.;
  double tmpaopHb = 0.;
  double tmpchHb = 0.;
  double tmppkHb = 0.;
  double tmprms = 0.;
  double tmprmspk = 0.;
  double tmprmsch = 0.;
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
    ++nrEventsRead;
    if (nrEventsRead%1000 == 0) 
    {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    bool found = false;
    IoSdEvent event(pos);

    if ( event.Id == previusEvent )
      continue;

    previusEvent = event.Id;
    sameUtc = event.utctime();

    if ( sameUtc > cday )
    {
      tmprms = 0;
      tmprmspk = 0.;
      tmprmsch = 0.;

      tmpaopHb = tmpaopHb/stckapHb.size();
      tmppkHb = tmppkHb/stckpkHb.size();
      tmpchHb = tmpchHb/stckchHb.size();

      for ( unsigned int kk=0; kk< stckapHb.size(); kk++ )
      {
        tmprms += (tmpaopHb - stckapHb[kk])*(tmpaopHb - stckapHb[kk]);
        tmprmspk += (tmppkHb - stckpkHb[kk])*(tmppkHb - stckpkHb[kk]);
        tmprmsch += (tmpchHb - stckchHb[kk])*(tmpchHb - stckchHb[kk]);
      }
      
      apHbase->SetBinContent( nday, tmpaopHb );
      pkHbase->SetBinContent( nday, tmppkHb );
      chHbase->SetBinContent( nday, tmpchHb );

      rmsHbase->SetBinContent( nday, sqrt( tmprms/stckapHb.size() ) );
      rmspkHbase->SetBinContent( nday, sqrt( tmprmspk/stckpkHb.size() ) );
      rmschHbase->SetBinContent( nday, sqrt( tmprmsch/stckchHb.size() ) );

      tmprms = 0.;
      timestmp->SetBinContent( nday, sameUtc );

      tmpaopHb = 0.;
      tmppkHb = 0.;
      tmpchHb = 0.;

      cday += dday;
      nday++;
      treeHist->Fill();
      stckapHb.clear();
      stckpkHb.clear();
      stckchHb.clear();
    }

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
          << " " << nrEventsRead-1 << " " << sameUtc
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
          if (fabs(meanl-meanf) < 2*rmsf)
          {

            // ================
            // *** Baseline ***
            blCorrHbase = event.Stations[i].HBase(pmtId-1)->GetMean();

            // =================
            // *** For HBase ***
            recePk = event.Stations[i].HPeak(pmtId-1);
            tmpName.Form("%d%d", event.UTCTime, nrEventsRead-1);
          
            fitPk.getCrrSmooth(*recePk, blCorrHbase, tmpName+"Hbpk");
            tmp = fitPk.getPkCorrSmooth();
            fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);

            if ( fitPk.chisPeak < 50. )
              pkChis->Fill( fitPk.chisPeak );
            else
              pkChis->Fill( 49. );
  
            //if ( fitPk.getValidHisto(*tmp) >= 3 )
              //fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
          
            receCh = event.Stations[i].HCharge(pmtId-1);
            fitCh.getChCrrSmooth(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
            tmp = fitCh.getChCorrSmooth();
            fitCh.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);
            if ( fitCh.chisCharge < 50. )
              chChis->Fill( fitCh.chisCharge );
            else
              chChis->Fill( 49. );

            pkHistFit = fitPk.getFitGraphPk();
            pkChi2 = fitPk.chisPeak;

            chHistFit = fitCh.getFitGraphCh();
            chChi2 = fitCh.chisCharge;

            treeHistCh2->Fill();
          
            //if ( fitCh.getValidHisto(*tmp) >= 2 )
              //fitCh.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);

      			// ============================
		      	// *** Saving A/P for HBase ***
				    if ( fitCh.fitChOk && fitPk.fitPkOk )
            {
  					  stckapHb.push_back( fitCh.vemPosCh/fitPk.vemPosPk );
              stckchHb.push_back( fitCh.vemPosCh );
              stckpkHb.push_back( fitPk.vemPosPk );
              tmpaopHb += fitCh.vemPosCh/fitPk.vemPosPk;
              tmpchHb += fitCh.vemPosCh;
              tmppkHb += fitPk.vemPosPk;
            }

            eventStat++;
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
