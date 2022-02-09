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
      << "4 For PMTSSD; " << "5 For SPMT" << endl;
    cout << "==================================================" << endl;
    exit(0);
	}
  
  cerr << "You have selected " << pmtname << endl;
	unsigned int totSt = stationsIds.size();

	if ( totSt==1 )
		pmtname = "results/St"+to_string( stationsIds[0] )+"pmt"+to_string( pmtId );
 
  string doMonth = string(whichmonth);
	TH1F *receCh = new TH1F (); // Receive Ch from IoSdStation::HCharge

  unsigned int previusEvent = 0; // Avoiding read the same event
  unsigned int nrEventsRead = 0;
  unsigned int nrEvents = 0;
  bool found = false;

  TString tmpName;
  fstream file;

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

        if (event.Stations[i].Error==256) {
          receCh = event.Stations[i].HCharge(pmtId-1);
          tmpName.Form("Utc%ld.dat", event.UTCTime);
          tmpName = pmtname+tmpName;
          file.open(tmpName, ios_base::out);
          for ( int bin_i=1; bin_i<receCh->GetNbinsX()+1; bin_i++ ) {
            file << receCh->GetBinCenter( bin_i ) - 4 - receCh->GetBinCenter(0) << " " 
              << receCh->GetBinContent( bin_i ) << endl;
          }
          file.close();
        }
      }
    }
  }
 
	return 0;
}