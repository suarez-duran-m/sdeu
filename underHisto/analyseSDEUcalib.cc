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
  if ( argc < 5 ) {
    cout << endl
      << "Usage: " << argv[0] << " <stationsFile>  <PMT>  <Month> <files>" << endl
      << "  <stationsFile>: file with a list of stations" << endl
      << "  <PMT>: ID for the PMT you want to analyse" << endl
      << "  <lrb>: N bins for fit (n towards left, n towards right)" << endl
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
  const char* nlrb = argv[3];
  const char* whichmonth = argv[4];
  AugerIoSd input(argc-5, argv+5);
  const unsigned int totalNrEvents = input.NumberOfEvents();
  ifstream stationsFile(stationsFileName, ios::in);

  if (!stationsFile.is_open()) {
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
  
  if (stationsIds.empty()) {
    cout << "Please specify the stations ids in the file " << endl;
    exit(0);
  }
  
  TString nameStati = to_string( stationsIds[0] );
  TString pmtname = whichpmt;
  TString strNblr = nlrb;
  int pmtId= atoi( pmtname );
  if ( pmtId > 0 && pmtId < 4 )
     pmtname = "PMT"+to_string( pmtId );
  else if ( pmtId == 4 )
    pmtname = "SPMT";
  else if ( pmtId == 5 )
    pmtname = "PMTSSD";
  else {
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
  pmtname +=  "lrb" + strNblr + doMonth + "2021";

  TFile hfile("uubChPk"+pmtname+".root","RECREATE","");
  //TFile hfile("kk.root", "RECREATE","");

	TH1F *recePk = new TH1F (); // Receive Pk from IoSdStation::HPeak
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

  TGraphErrors *pkHistFit = new TGraphErrors();
  TH1F *pkForFit = new TH1F();
  double pkChi2 = 0.;
  int pkNdf = 0.;
  double pkProb = 0.;
  double peak = 0.;
  double peakDeri = 0.;
  double pkLow = 0.;
  double pkHigh = 0.;
  double pkPar0 = 0.; 
  double pkPar1 = 0.; 
  double pkPar2 = 0.; 

  TGraphErrors *chHistFit = new TGraphErrors();
  TH1F *chForFit = new TH1F();
  double chChi2 = 0.;
  int chNdf = 0.;
  double chProb = 0.;
  double charge = 0.;
  double chargeDeri = 0.;
  double chLow = 0.;
  double chHigh = 0.;
  double chPar0 = 0.; 
  double chPar1 = 0.; 
  double chPar2 = 0.; 

  unsigned int evtIdPk = 0; //Storing event Id
  unsigned int evtTimePk = 0; //Storing day-Unixtime
  unsigned int evtIdCh = 0; //Storing event Id
  unsigned int evtTimeCh = 0; //Storing day-Unixtime
  unsigned int previusEvent = 0; // Avoiding read the same event
  unsigned int nrEventsRead = 0;
  unsigned int nrEvents = 0;
  bool found = false;
	double blCorrHbase = 0.;

  TH1F *tmp = new TH1F();
  TString tmpName;
  fitpeak fitPk;
  fitcharge fitCh;

  // ======================================
  // *** *** *** Creating Trees *** *** *** 
  TTree *treePeak = new TTree("PeakData","");
  TTree *treeCharge = new TTree("ChargeData","");

  treePeak->Branch("graph","TGraphErrors",&pkHistFit,32000,0);
  treePeak->Branch("peakForFit","TH1F",&pkForFit,32000,0);
  treePeak->Branch("peakVal",&peak,"peak/D");
  treePeak->Branch("peakValDer",&peakDeri,"peakDeri/D");
  treePeak->Branch("chi2",&pkChi2,"pkChi2/D");
  treePeak->Branch("ndf",&pkNdf,"pkNdf/I");
  treePeak->Branch("prob",&pkProb,"pkProb/D");
  treePeak->Branch("low",&pkLow,"pkLow/D");
  treePeak->Branch("high",&pkHigh,"pkHigh/D");
  treePeak->Branch("pkPar0",&pkPar0,"pkPar0/D");
  treePeak->Branch("pkPar1",&pkPar1,"pkPar1/D");
  treePeak->Branch("pkPar2",&pkPar2,"pkPar2/D");
  treePeak->Branch("eventId",&evtIdPk,"evtIdPk/I");
  treePeak->Branch("timeEvnt",&evtTimePk,"evtTimePk/I");

  treeCharge->Branch("graph","TGraphErrors",&chHistFit,32000,0);
  treeCharge->Branch("chargeForFit","TH1F",&chForFit,32000,0);
  treeCharge->Branch("chargeVal",&charge,"charge/D");
  treeCharge->Branch("chargeValDer",&chargeDeri,"chargeDeri/D");
  treeCharge->Branch("chi2",&chChi2,"chChi2/D");
  treeCharge->Branch("ndf",&chNdf,"chNdf/I");
  treeCharge->Branch("prob",&chProb,"chProb/D");
  treeCharge->Branch("low",&chLow,"chLow/D");
  treeCharge->Branch("high",&chHigh,"chHigh/D");
  treeCharge->Branch("chPar0",&chPar0,"chPar0/D");
  treeCharge->Branch("chPar1",&chPar1,"chPar1/D");
  treeCharge->Branch("chPar2",&chPar2,"chPar2/D");
  treeCharge->Branch("eventId",&evtIdCh,"evtIdCh/I");
  treeCharge->Branch("timeEvnt",&evtTimeCh,"evtTimeCh/I");

  EventPos pos;

  for (pos=input.FirstEvent(); pos<input.LastEvent(); pos=input.NextEvent()) {
    nrEventsRead++;
    if (nrEventsRead%1000 == 0) {
      cout << "====> Read " << nrEventsRead << " out of " << totalNrEvents << endl;
      cout << "      Wrote: " << nrEvents << " events" << endl;
    }

    IoSdEvent event(pos);
    if ( event.Id == previusEvent )
      continue;

    previusEvent = event.Id;

    for (unsigned int i = 0 ; i < event.Stations.size(); ++i) {
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
          iter!= stationsIds.end(); ++iter)
        if (event.Stations[i].Id == *iter )
          found = true;
      if ( !found )
        continue;
      
      if ( event.Stations[i].IsUUB ) {
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
          << " " << nrEventsRead-1
          << endl;

        //cout << "MSD: " << event.Stations[i].calib()->VemCharge[pmtId-1] << endl;

        if (event.Stations[i].Error==256) {
          tmpName.Form("%d%d", event.UTCTime, nrEventsRead-1);
          blCorrHbase = event.Stations[i].HBase(pmtId-1)->GetMean(); // Extracting calib-baseline
          
          recePk = event.Stations[i].HPeak(pmtId-1); // Receiving Peak histogram
          
          blCorrHbase = ( fabs(blCorrHbase - recePk->GetBinCenter(1) < 20 ) ? (recePk->GetBinCenter(1)):0 ); // From OffLine
          fitPk.getCrr(*recePk, blCorrHbase, tmpName+"Hbpk"); // Correcting for calib-baseline
          tmp = fitPk.getPkCorr(); // Receiving corrected histogram
          pkForFit = fitPk.getPkCorr();
          fitPk.getFitPk(*tmp); // Fitting

          pkHistFit = fitPk.getFitGraphPk();
          pkChi2 = fitPk.chisPeak;
          pkNdf = fitPk.ndfPeak;
          pkProb = fitPk.probPeak;
          peak = fitPk.vemPosPk;
          peakDeri = fitPk.vemFromDeri;
          pkLow = fitPk.rangXmin;
          pkHigh = fitPk.rangXmax;
          pkPar0 = fitPk.par0;
          pkPar1 = fitPk.par1;
          pkPar2 = fitPk.par2;
            
          receCh = event.Stations[i].HCharge(pmtId-1);
          fitCh.setChCrr(*receCh, event.Stations[i].Histo->Offset[pmtId-1+6], tmpName+"Hbch");
          tmp = fitCh.getChCrr();
          chForFit = fitCh.getChCrr();
          fitCh.getFitCh(*tmp, atoi(argv[3]) );
          
          chHistFit = fitCh.getFitGraphCh();
          chChi2 = fitCh.chisCharge;
          chNdf = fitCh.ndfCharge;
          chProb = fitCh.probCharge;
          charge =  fitCh.vemPosCh;
          chargeDeri =  fitCh.vemPosDeri;
          chLow = fitCh.rangXmin;
          chHigh = fitCh.rangXmax;
          chPar0 = fitCh.par0;
          chPar1 = fitCh.par1; 
          chPar2 = fitCh.par2;
          
          evtIdPk = event.Id;
          evtIdCh = event.Id;
          evtTimePk = event.UTCTime - 315964782;
          evtTimeCh = event.UTCTime - 315964782;

          treePeak->Fill();
          treeCharge->Fill();

          
          break; // Apply if the is running for a single station.
        }
      }
    }
  }

  hfile.Write();
  hfile.Close();
  
	return 0;
}
