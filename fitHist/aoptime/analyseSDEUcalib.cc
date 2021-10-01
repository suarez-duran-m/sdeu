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
  vector < double > stckchHb; // It stores ap Hb per day
  vector < double > stckpkHb; // It stores ap Hb per day
  vector < double > stckapCa; // It stores ap Ca per day
	
	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

	TH1F *apHbase = new TH1F("apHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *pkHbase = new TH1F("pkHbase","",totDays, 0, totDays);//It stores VemP from HBase for St.
	TH1F *chHbase = new TH1F("chHbase","",totDays, 0, totDays);//It stores VemQ from HBase for St.
	TH1F *rmsHbase = new TH1F("rmsHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *rmspkHbase = new TH1F("rmspkHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *rmschHbase = new TH1F("rmschHbase","",totDays, 0, totDays);//It stores A/P from HBase for St.
	TH1F *apCalib = new TH1F("apCalib","",totDays, 0, totDays);//It stores A/P from Calib for St.
	TH1F *rmsCalib = new TH1F("rmsCalib","",totDays, 0, totDays);//It stores A/P from Calib for St.
	TH1F *timestmp = new TH1F("timestmp","",totDays, 0, totDays);//It stores A/P from Calib for St.
/*
	TH1F *chisHbasePk = new TH1F("chisHbasePk","",totDays, 0, totDays);//It stores chis from fitted peak for St.
	TH1F *chisCalibPk = new TH1F("chisCalibPk","",totDays, 0, totDays);//It stores chis from fitted peak for St.
	TH1F *chisCh = new TH1F("chisCh","",totDays, 0, totDays);//It stores chis from fitted charge for St.

  TH1F *okFitPkHb = new TH1F("okFitPkHb","",totDays, 0, totDays);//It stores number Fit Ok for charge per St. 
  TH1F *okFitPkCa = new TH1F("okFitPkCa","",totDays, 0, totDays);//It stores number Fit Ok for charge per St.
  TH1F *okFitCh = new TH1F("okFitCh","",totDays, 0, totDays);//It stores number Fit Ok for charge per St.

  TH1F *nslopsHbPk = new TH1F("nslopsHbPk","",totDays, 0, totDays);//It stores peaks in Pk Histo. per St, for Hbase.
  TH1F *nslopsCaPk = new TH1F("nslopsCaPk","",totDays, 0, totDays);//It stores pekas Pk Histo. per St, for Calib.
  TH1F *nslopsCh = new TH1F("nslopsCh", "",totDays, 0, totDays);
*/
	unsigned int eventStat = 0; //It stores the number of events.
	unsigned int evtTime = 0; //It stores the day-Unixtime

  TTree *treeHist = new TTree("Histograms","");
  /*
	treeHist->Branch("receCh", "TH1F", &receCh);
	treeHist->Branch("recePk", "TH1F", &recePk);
  */
	//treeHist->Branch("apHbase", "TH1F", &apHbase);
	//treeHist->Branch("apCalib", "TH1F", &apCalib);
  /*
	treeHist->Branch("chisHbasePk", "TH1F", &chisHbasePk);
	treeHist->Branch("chisCalibPk", "TH1F", &chisCalibPk);
	treeHist->Branch("chisCh", "TH1F", &chisCh);
	
  treeHist->Branch("okFitPkHb", "TH1F", &okFitPkHb);
  treeHist->Branch("okFitPkCa", "TH1F", &okFitPkCa);
	treeHist->Branch("okFitCh", "TH1F", &okFitCh);

  treeHist->Branch("nslopsHbPk","TH1F",&nslopsHbPk);
  treeHist->Branch("nslopsCaPk","TH1F",&nslopsCaPk);
  treeHist->Branch("nslopsCh","TH1F",&nslopsCh);
  */  
	treeHist->Branch("eventStat", &eventStat, "eventStat/I");
	treeHist->Branch("evtTime",&evtTime,"evtTime/I");

	double blCorrHbase = 0.;
	double blCorrCalib = 0.;
  double tmpaopHb = 0.;
  double tmpchHb = 0.;
  double tmppkHb = 0.;
  double tmpaopCa = 0.;
  double tmprms = 0.;
  double tmprmspk = 0.;
  double tmprmsch = 0.;
  TH1F *tmp = new TH1F();
  TString tmpName;
  fitpeak fitPk;
  fitcharge fitCh;

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
      tmpaopCa = tmpaopCa/stckapCa.size();
      tmppkHb = tmppkHb/stckpkHb.size();
      tmpchHb = tmpchHb/stckchHb.size();

      //for ( vector< double >::iterator evtap=stckapHb.begin(); evtap!=stckapHb.end(); evtap ++ )
      for ( int kk=0; kk< stckapHb.size(); kk++ )
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
      /*
      for ( vector< double >::iterator evtap=stckapCa.begin(); evtap!=stckapCa.end(); evtap ++ )
        tmprms += (tmpaopCa - *evtap)*(tmpaopCa - *evtap);
      
      apCalib->SetBinContent( nday, tmpaopCa );
      rmsCalib->SetBinContent( nday, sqrt( tmprms/stckapCa.size() ) );
      */

      timestmp->SetBinContent( nday, sameUtc );

      tmpaopHb = 0.;
      tmppkHb = 0.;
      tmpchHb = 0.;
      tmpaopCa = 0.;

      cday += dday;
      nday++;
      treeHist->Fill();
      stckapHb.clear();
      stckpkHb.clear();
      stckchHb.clear();
      stckapCa.clear();
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
          // ================
          // *** Baseline ***
          blCorrHbase = event.Stations[i].HBase(pmtId-1)->GetMean();
          blCorrCalib = event.Stations[i].Calib->Base[pmtId-1];

          // =================
          // *** For HBase ***
          recePk = event.Stations[i].HPeak(pmtId-1);
          tmpName.Form("%d%d", event.UTCTime, nrEventsRead-1);
          
          fitPk.getCrrSmooth(*recePk, blCorrHbase, tmpName+"Hbpk");
          tmp = fitPk.getPkCorrSmooth();

          if ( fitPk.getValidHisto(*tmp) >= 3 )
            fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
          
          receCh = event.Stations[i].HCharge(pmtId-1);
          fitCh.getChCrrSmooth(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hbch");
          tmp = fitCh.getChCorrSmooth();
          
          if ( fitCh.getValidHisto(*tmp) >= 2 )
            fitCh.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);

          //cerr << fitCh.vemPosCh << " " << fitPk.vemPosPk << " "
            //fitCh.vemPosCh/fitPk.vemPosPk << endl;

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

          /*
          
          // =================
          // *** For Calib ***
          fitPk.getCrrSmooth(*recePk, blCorrCalib, tmpName+"Capk");
          tmp = fitPk.getPkCorrSmooth();
          if ( fitPk.getValidHisto(*tmp) >= 3 )
          {
            //nslopsCaPk->SetBinContent( stid, fitPk.getValidHisto(*tmp) );
            fitPk.getFitPk(*tmp, 0.2, 10, event.Stations[i].Calib->VemPeak[pmtId-1]);
          }

    			// ============================
		    	// *** Saving A/P for Calib ***
				  if ( fitCh.fitChOk && fitPk.fitPkOk )
          {
					  stckapCa.push_back( fitCh.vemPosCh/fitPk.vemPosPk );
            tmpaopCa += fitCh.vemPosCh/fitPk.vemPosPk;
          }
          */

          eventStat ++;
					break;
        }
      }
    }
  }

  hfile.Write();
  hfile.Close();
  
	return 0;
}
