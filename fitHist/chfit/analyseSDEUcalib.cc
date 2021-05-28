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

  sort(stationsIds.begin(), stationsIds.end());
  //double stationsBins[totSt];
  for ( unsigned int i=0; i<totSt; i++)
  {
    //stationsBins[i] = stationsIds[i];
    cout << i << " -> " << stationsIds[i] << endl;
  }

	if ( totSt==1 )
		pmtname += "St"+to_string( stationsIds[0] );
 
  string doMonth = string(whichmonth);
  pmtname += doMonth;
  TFile hfile("uubCharge"+pmtname+".root","RECREATE","");

  unsigned int totDays = 0; // Total month's day
	unsigned int cday = 0; // Unixtime for 2nd month's day

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

  unsigned int previusEvent = 0; // It avoids read the same event twice
  unsigned int sameUtc = 0; // It allows to know when a full day has been read
	unsigned int nday = 0; // It counts the number of days
	unsigned int dday = 86400; // Day duration in sec

  TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HPeak
 
  // ====================
  // *** For Averages ***
  vector < vector < double > > stckvemchHb; // It stores Vem Ch Hb per day per St 
  vector < vector < double > > stckcntvemchHb; // It stores counts at vem Ch for Hb per St
  vector < vector < double > > stckChiHb; // It stores chi^2/NDF for fit Hb per St

  stckvemchHb.resize( totSt );
  stckcntvemchHb.resize( totSt );
  stckChiHb.resize( totSt );
  
  TH1F *vemchHb = new TH1F("vemchHb","",totSt, 0, totSt);//It stores average Vem Ch HB per St per day
  TH1F *cntvemchHb = new TH1F("cntvemchHb","",totSt, 0, totSt);//It stores average conts Vem Ch HB per St per day
  TH1F *chiHb = new TH1F("chiHb","",totSt, 0, totSt);//It stores average chi^2/NDF for fit Hb per St per day

  // ===============
  // *** For RMS ***
  
  TH1F *vemchHbrms = new TH1F("vemchHbrms","",totSt, 0, totSt);//It stores rms Vem Ch HB per St per day
  TH1F *cntvemchHbrms = new TH1F("cntvemchHbrms","",totSt, 0, totSt);//It stores rms conts Vem Ch HB per St per day
  TH1F *chiHbrms = new TH1F("chiHbrms","",totSt, 0, totSt);//It stores rms chi^2/NDF for fit Hb per St per day
 
  // ==================
  // *** Time label ***
	TH1F *timestmp = new TH1F("timestmp","", totSt, 0, totSt);//It stores Timestamp for Hb per St
	TH1F *eventStat = new TH1F("eventStat","", totSt, 0, totSt);//It stores the number of events per St per day

  // =====================
  // *** Creating Tree ***

  TTree *treeHist = new TTree("Histograms","");

	treeHist->Branch("vemchHb", "TH1F", &vemchHb);
  treeHist->Branch("cntvemchHb", "TH1F", &cntvemchHb);
  treeHist->Branch("chiHb", "TH1F", &chiHb);

  treeHist->Branch("vemchHbrms", "TH1F", &vemchHbrms);
  treeHist->Branch("cntvemchHbrms", "TH1F", &cntvemchHbrms);
  treeHist->Branch("chiHbrms", "TH1F", &chiHbrms);

  treeHist->Branch("timestmp", "TH1F", &timestmp);
  treeHist->Branch("eventStat", "TH1F", &eventStat);

  // =================
  // *** Temporals ***

  unsigned int nrEventsRead = 0;
  unsigned int nrEvents = 0;

  unsigned int nblbins = 100;
  double meanf = 0.;
  double meanl = 0.;
  double rmsf = 0.;

  int tmpcnt= 0; // To find the Is of station running on
  unsigned int currentSt = 0;// Id to identify the station running on
  TH1F *tmp = new TH1F();
  TString tmpName;
  fitcharge fitCh;
  vector < int > *blpmth  = new vector < int >;

  vector < double > tmpvemHb(totSt);
  vector < int > tmpcntHb(totSt);
  vector < double > tmpChiHb(totSt);

  vector < int > tmpeventStat(totSt);
  
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
      for ( unsigned int st = 0; st < stationsIds.size(); st++ )
      {  
        if ( stckvemchHb[st].size() > 0 )
        {
          tmpvemHb[st] = tmpvemHb[st] / stckvemchHb[st].size();
          tmpcntHb[st] = tmpcntHb[st] / stckvemchHb[st].size();
          tmpChiHb[st] = tmpChiHb[st] / stckvemchHb[st].size();

          vemchHb->SetBinContent( st, tmpvemHb[st] );
          cntvemchHb->SetBinContent( st, tmpcntHb[st] );
          chiHb->SetBinContent( st, tmpChiHb[st] );

          vemchHbrms->SetBinContent( st, fitCh.getRms( stckvemchHb, tmpvemHb, st) );
          cntvemchHbrms->SetBinContent(st,nday,fitCh.getRms(stckcntvemchHb,tmpcntHb,st));
          chiHbrms->SetBinContent( st, fitCh.getRms( stckChiHb, tmpChiHb, st ) );

          eventStat->SetBinContent( st, (1.*stckvemchHb[st].size())/(1.*tmpeventStat[st]) );
        }
        
        timestmp->SetBinContent( st, cday );

        tmpvemHb[st] = 0;
        tmpcntHb[st] = 0;
        tmpChiHb[st] = 0;

        tmpeventStat[st] = 0;

        stckvemchHb[st].clear();
        stckcntvemchHb[st].clear();
        stckChiHb[st].clear();
      }
      cday += dday;
      nday++;
      treeHist->Fill();
    }

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i)
    {
      found = false;
      tmpcnt = 0;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter )
      {
        if ( event.Stations[i].Id == *iter )
        {
          currentSt = tmpcnt;
          found = true;
        }
        tmpcnt++;
      }
      if ( !found )
        continue;

      cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
        << " " << nrEventsRead-1 << " " << sameUtc  << endl;

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
          // =================
          // *** For HBase ***
          receCh = event.Stations[i].HCharge(pmtId-1);
          tmpName.Form("%d%d%d", int(event.UTCTime), nrEventsRead-1, currentSt);
        
          fitCh.getCrrSmooth(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6]/20., tmpName+"Hb");
          tmp = fitCh.getChCorrSmooth();
          fitCh.getFitCh(*tmp, 0.3, 10, event.Stations[i].Calib->VemCharge[pmtId-1]);
  
      		// =======================
		      // *** Stacking Values ***
			  	if ( fitCh.fitChOk )
          {
            stckvemchHb[currentSt].push_back( fitCh.vemPosCh ); 
            stckcntvemchHb[currentSt].push_back( fitCh.cntvemCh );
            stckChiHb[currentSt].push_back( fitCh.chisCharge );
            tmpvemHb[currentSt] += fitCh.vemPosCh;
            tmpcntHb[currentSt] += fitCh.cntvemCh;
            tmpChiHb[currentSt] += fitCh.chisCharge;
          }
          tmpeventStat[currentSt]++;
		  		break;
        }
      }
    }
  }

  hfile.Write();
  hfile.Close();
  
	return 0;
}
