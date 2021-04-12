#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#include <IoAuger.h>
#include <Ec.h>

#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TCanvas.h>


using namespace std;

float getmean( int arr[] ){
  float mean = 0.;
  for (int i=0; i<100; i++)
    mean += arr[i];
  return mean/100.;
}

float getrms( int arr[], float meanarr ){
  float rms = 0.;
  for (int i=0; i<100; i++)
    rms += (arr[i] - meanarr)*(arr[i] - meanarr);
  return rms/100.;
}


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
  const char* outputName = argv[2];
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
    
  IoSd *sdFile = NULL;
  AugerFile *adFile = NULL;
  TString nameStati = to_string( stationsIds[0] );

  TFile hfile("bl"+nameStati+".root","RECREATE","");

  TH1I eventId("eventId", "Events IDs for Station "+nameStati, 5356800, 1606780800, 1612137600);

  //TTree tree("T","");

  // For PMT 1
  TH1F pmt1lbmeanf("pmt1lbmeanf","Station "+nameStati+" Mean First 100 bins PMT 1 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt1lbmeanl("pmt1lbmeanl","Station "+nameStati+" Mean Last 100 bins PMT 1 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt1hbmeanf("pmt1hbmeanf","Station "+nameStati+" Mean First 100 PMT 1 High", 5356800, 1606780800, 1612137600);
  TH1F pmt1hbmeanl("pmt1hbmeanl","Station "+nameStati+" Mean Last 100 PMT 1 High", 5356800, 1606780800, 1612137600);

  TH1F pmt1lbrmsf("pmt1lbrmsf","Station "+nameStati+" RMS First 100 bins PMT 1 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt1lbrmsl("pmt1lbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 1 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt1hbrmsf("pmt1hbrmsf","Station "+nameStati+" RMS First 100 bins PMT 1 High", 5356800, 1606780800, 1612137600);
  TH1F pmt1hbrmsl("pmt1hbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 1 High", 5356800, 1606780800, 1612137600);

  // For PMT 2
  TH1F pmt2lbmeanf("pmt2lbmeanf","Station "+nameStati+" Mean First 100 bins PMT 2 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt2lbmeanl("pmt2lbmeanl","Station "+nameStati+" Mean Last 100 bins PMT 2 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt2hbmeanf("pmt2hbmeanf","Station "+nameStati+" Mean First 100 PMT 2 High", 5356800, 1606780800, 1612137600);
  TH1F pmt2hbmeanl("pmt2hbmeanl","Station "+nameStati+" Mean Last 100 PMT 2 High", 5356800, 1606780800, 1612137600);

  TH1F pmt2lbrmsf("pmt2lbrmsf","Station "+nameStati+" RMS First 100 bins PMT 2 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt2lbrmsl("pmt2lbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 2 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt2hbrmsf("pmt2hbrmsf","Station "+nameStati+" RMS First 100 bins PMT 2 High", 5356800, 1606780800, 1612137600);
  TH1F pmt2hbrmsl("pmt2hbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 2 High", 5356800, 1606780800, 1612137600);

  // For PMT 3
  TH1F pmt3lbmeanf("pmt3lbmeanf","Station "+nameStati+" Mean First 100 bins PMT 3 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt3lbmeanl("pmt3lbmeanl","Station "+nameStati+" Mean Last 100 bins PMT 3 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt3hbmeanf("pmt3hbmeanf","Station "+nameStati+" Mean First 100 PMT 3 High", 5356800, 1606780800, 1612137600);
  TH1F pmt3hbmeanl("pmt3hbmeanl","Station "+nameStati+" Mean Last 100 PMT 3 High", 5356800, 1606780800, 1612137600);

  TH1F pmt3lbrmsf("pmt3lbrmsf","Station "+nameStati+" RMS First 100 bins PMT 3 Low", 5356800, 1606780800, 1612137600);
  TH1F pmt3lbrmsl("pmt3lbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 3 Low", 5356800, 1606780800, 1612137600);

  TH1F pmt3hbrmsf("pmt3hbrmsf","Station "+nameStati+" RMS First 100 bins PMT 3 High", 5356800, 1606780800, 1612137600);
  TH1F pmt3hbrmsl("pmt3hbrmsl","Station "+nameStati+" RMS Last 100 bins PMT 3 High", 5356800, 1606780800, 1612137600);

  // For SPMT 
  TH1F spmtlbmeanf("spmtlbmeanf","Station "+nameStati+" Mean First 100 bins SPMT Low", 5356800, 1606780800, 1612137600);
  TH1F spmtlbmeanl("spmtlbmeanl","Station "+nameStati+" Mean Last 100 bins SPMT Low", 5356800, 1606780800, 1612137600);

  TH1F spmtlbrmsf("spmtlbrmsf","Station "+nameStati+" RMS First 100 bins SPMT Low", 5356800, 1606780800, 1612137600);
  TH1F spmtlbrmsl("spmtlbrmsl","Station "+nameStati+" RMS Last 100 bins SPMT Low", 5356800, 1606780800, 1612137600);

  // For SSDPMT
  TH1F ssdpmtlbmeanf("ssdpmtlbmeanf","Station "+nameStati+" Mean First 100 bins SSDPMT Low", 5356800, 1606780800, 1612137600);
  TH1F ssdpmtlbmeanl("ssdpmtlbmeanl","Station "+nameStati+" Mean Last 100 bins SSDPMT Low", 5356800, 1606780800, 1612137600);

  TH1F ssdpmthbmeanf("ssdpmthbmeanf","Station "+nameStati+" Mean First 100 SSDPMT High", 5356800, 1606780800, 1612137600);
  TH1F ssdpmthbmeanl("ssdpmthbmeanl","Station "+nameStati+" Mean Last 100 SSDPMT High", 5356800, 1606780800, 1612137600);

  TH1F ssdpmtlbrmsf("ssdpmtlbrmsf","Station "+nameStati+" RMS First 100 bins SSDPMT Low", 5356800, 1606780800, 1612137600);
  TH1F ssdpmtlbrmsl("ssdpmtlbrmsl","Station "+nameStati+" RMS Last 100 bins SSDPMT Low", 5356800, 1606780800, 1612137600);

  TH1F ssdpmthbrmsf("ssdpmthbrmsf","Station "+nameStati+" RMS First 100 bins SSDPMT High", 5356800, 1606780800, 1612137600);
  TH1F ssdpmthbrmsl("ssdpmthbrmsl","Station "+nameStati+" RMS Last 100 bins SSDPMT High", 5356800, 1606780800, 1612137600);


  int mpmt1lbmeanf[100];
  int mpmt1lbmeanl[100];

  int mpmt1hbmeanf[100];
  int mpmt1hbmeanl[100];

  int mpmt2lbmeanf[100];
  int mpmt2lbmeanl[100];

  int mpmt2hbmeanf[100];
  int mpmt2hbmeanl[100];

  int mpmt3lbmeanf[100];
  int mpmt3lbmeanl[100];

  int mpmt3hbmeanf[100];
  int mpmt3hbmeanl[100];

  int mspmtlbmeanf[100];
  int mspmtlbmeanl[100];

  int mssdpmtlbmeanf[100];
  int mssdpmtlbmeanl[100];

  int mssdpmthbmeanf[100];
  int mssdpmthbmeanl[100];


  int previusEvent = 0;
  int sameUtc = 0;
  int nsample = 0;

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
    for (unsigned int i = 0 ; i < event.Stations.size(); ++i){
      found = false;
      for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
             iter!= stationsIds.end(); ++iter){
        if (event.Stations[i].Id == *iter){
          found = true;
          }
        }
      
      if (!found)
        continue;
      if ( event.Stations[i].IsUUB && event.Id != previusEvent && sameUtc != event.utctime() ) {
        previusEvent = event.Id;
        sameUtc = event.utctime();

        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
             << endl;
        //cout << "# Error " << event.Stations[i].Error << endl;        

        IoSdEvent event(pos);
        
        if (event.Stations[i].Error==256) { //0+256
          nsample = event.Stations[i].UFadc->NSample;
          for (int k=0;k<100;k++){

            mpmt1lbmeanf[k] = event.Stations[i].UFadc->GetValue(0,1,k);
            mpmt1lbmeanl[k] = event.Stations[i].UFadc->GetValue(0,1,(nsample - k));

            mpmt1hbmeanf[k] = event.Stations[i].UFadc->GetValue(0,0,k);
            mpmt1hbmeanl[k] = event.Stations[i].UFadc->GetValue(0,0,(nsample - k));

            mpmt2lbmeanf[k] = event.Stations[i].UFadc->GetValue(1,1,k);
            mpmt2lbmeanl[k] = event.Stations[i].UFadc->GetValue(1,1,(nsample - k));

            mpmt2hbmeanf[k] = event.Stations[i].UFadc->GetValue(1,0,k);
            mpmt2hbmeanl[k] = event.Stations[i].UFadc->GetValue(1,0,(nsample - k));

            mpmt3lbmeanf[k] = event.Stations[i].UFadc->GetValue(2,1,k);
            mpmt3lbmeanl[k] = event.Stations[i].UFadc->GetValue(2,1,(nsample - k));

            mpmt3hbmeanf[k] = event.Stations[i].UFadc->GetValue(2,0,k);
            mpmt3hbmeanl[k] = event.Stations[i].UFadc->GetValue(2,0,(nsample - k));

            mspmtlbmeanf[k] = event.Stations[i].UFadc->GetValue(3,1,k);
            mspmtlbmeanl[k] = event.Stations[i].UFadc->GetValue(3,1,(nsample - k));

            mssdpmtlbmeanf[k] = event.Stations[i].UFadc->GetValue(4,1,k);
            mssdpmtlbmeanl[k] = event.Stations[i].UFadc->GetValue(4,1,(nsample - k));

            mssdpmthbmeanf[k] = event.Stations[i].UFadc->GetValue(4,0,k);
            mssdpmthbmeanl[k] = event.Stations[i].UFadc->GetValue(4,0,(nsample - k));
          }

          eventId.Fill(event.utctime(), event.Id);

          pmt1lbmeanf.Fill(event.utctime(), getmean(mpmt1lbmeanf));
          pmt1lbmeanl.Fill(event.utctime(), getmean(mpmt1lbmeanl));
          pmt1lbrmsf.Fill(event.utctime(), getrms(mpmt1lbmeanf, getmean(mpmt1lbmeanf)));
          pmt1lbrmsl.Fill(event.utctime(), getrms(mpmt1lbmeanl, getmean(mpmt1lbmeanl)));

          pmt1hbmeanf.Fill(event.utctime(), getmean(mpmt1hbmeanf));
          pmt1hbmeanl.Fill(event.utctime(), getmean(mpmt1hbmeanl));
          pmt1hbrmsf.Fill(event.utctime(), getrms(mpmt1hbmeanf, getmean(mpmt1hbmeanf)));
          pmt1hbrmsl.Fill(event.utctime(), getrms(mpmt1hbmeanl, getmean(mpmt1hbmeanl)));

          pmt2lbmeanf.Fill(event.utctime(), getmean(mpmt2lbmeanf));
          pmt2lbmeanl.Fill(event.utctime(), getmean(mpmt2lbmeanl));
          pmt2lbrmsf.Fill(event.utctime(), getrms(mpmt2lbmeanf, getmean(mpmt2lbmeanf)));
          pmt2lbrmsl.Fill(event.utctime(), getrms(mpmt2lbmeanl, getmean(mpmt2lbmeanl)));

          pmt2hbmeanf.Fill(event.utctime(), getmean(mpmt2hbmeanf));
          pmt2hbmeanl.Fill(event.utctime(), getmean(mpmt2hbmeanl));
          pmt2hbrmsf.Fill(event.utctime(), getrms(mpmt2hbmeanf, getmean(mpmt2hbmeanf)));
          pmt2hbrmsl.Fill(event.utctime(), getrms(mpmt2hbmeanl, getmean(mpmt2hbmeanl)));

          pmt3lbmeanf.Fill(event.utctime(), getmean(mpmt3lbmeanf));
          pmt3lbmeanl.Fill(event.utctime(), getmean(mpmt3lbmeanl));
          pmt3lbrmsf.Fill(event.utctime(), getrms(mpmt3lbmeanf, getmean(mpmt3lbmeanf)));
          pmt3lbrmsl.Fill(event.utctime(), getrms(mpmt3lbmeanl, getmean(mpmt3lbmeanl)));

          pmt3hbmeanf.Fill(event.utctime(), getmean(mpmt3hbmeanf));
          pmt3hbmeanl.Fill(event.utctime(), getmean(mpmt3hbmeanl));
          pmt3hbrmsf.Fill(event.utctime(), getrms(mpmt3hbmeanf, getmean(mpmt3hbmeanf)));
          pmt3hbrmsl.Fill(event.utctime(), getrms(mpmt3hbmeanl, getmean(mpmt3hbmeanl)));

          spmtlbmeanf.Fill(event.utctime(), getmean(mspmtlbmeanf));
          spmtlbmeanl.Fill(event.utctime(), getmean(mspmtlbmeanl));
          spmtlbrmsf.Fill(event.utctime(), getrms(mspmtlbmeanf, getmean(mspmtlbmeanf)));
          spmtlbrmsl.Fill(event.utctime(), getrms(mspmtlbmeanl, getmean(mspmtlbmeanl)));

          ssdpmtlbmeanf.Fill(event.utctime(), getmean(mssdpmtlbmeanf));
          ssdpmtlbmeanl.Fill(event.utctime(), getmean(mssdpmtlbmeanl));
          ssdpmtlbrmsf.Fill(event.utctime(), getrms(mssdpmtlbmeanf, getmean(mssdpmtlbmeanf)));
          ssdpmtlbrmsl.Fill(event.utctime(), getrms(mssdpmtlbmeanl, getmean(mssdpmtlbmeanl)));

          ssdpmthbmeanf.Fill(event.utctime(), getmean(mssdpmthbmeanf));
          ssdpmthbmeanl.Fill(event.utctime(), getmean(mssdpmthbmeanl));
          ssdpmthbrmsf.Fill(event.utctime(), getrms(mssdpmthbmeanf, getmean(mssdpmthbmeanf)));
          ssdpmthbrmsl.Fill(event.utctime(), getrms(mssdpmthbmeanl, getmean(mssdpmthbmeanl)));
        }
      }
    }
  }

  pmt1lbmeanf.GetXaxis()->SetTimeDisplay(1);
  pmt1lbmeanf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt1lbmeanf.GetXaxis()->SetLabelOffset(0.020);
  pmt1lbmeanf.GetXaxis()->SetLabelSize(0.03);
  pmt1lbmeanf.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) Low / FADC");
  pmt1lbmeanf.Draw("P");

  pmt1hbmeanl.GetXaxis()->SetTimeDisplay(1);
  pmt1hbmeanl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt1hbmeanl.GetXaxis()->SetLabelOffset(0.020);
  pmt1hbmeanl.GetXaxis()->SetLabelSize(0.03);
  pmt1hbmeanl.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) High / FADC");
  pmt1hbmeanl.Draw("P");

  pmt2lbmeanf.GetXaxis()->SetTimeDisplay(1);
  pmt2lbmeanf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt2lbmeanf.GetXaxis()->SetLabelOffset(0.020);
  pmt2lbmeanf.GetXaxis()->SetLabelSize(0.03);
  pmt2lbmeanf.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) Low / FADC");
  pmt2lbmeanf.Draw("P");

  pmt2hbmeanl.GetXaxis()->SetTimeDisplay(1);
  pmt2hbmeanl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt2hbmeanl.GetXaxis()->SetLabelOffset(0.020);
  pmt2hbmeanl.GetXaxis()->SetLabelSize(0.03);
  pmt2hbmeanl.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) High / FADC");
  pmt2hbmeanl.Draw("P");

  pmt3lbmeanf.GetXaxis()->SetTimeDisplay(1);
  pmt3lbmeanf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt3lbmeanf.GetXaxis()->SetLabelOffset(0.020);
  pmt3lbmeanf.GetXaxis()->SetLabelSize(0.03);
  pmt3lbmeanf.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) Low / FADC");
  pmt3lbmeanf.Draw("P");

  pmt3hbmeanl.GetXaxis()->SetTimeDisplay(1);
  pmt3hbmeanl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt3hbmeanl.GetXaxis()->SetLabelOffset(0.020);
  pmt3hbmeanl.GetXaxis()->SetLabelSize(0.03);
  pmt3hbmeanl.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) High / FADC");
  pmt3hbmeanl.Draw("P");

  spmtlbmeanf.GetXaxis()->SetTimeDisplay(1);
  spmtlbmeanf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  spmtlbmeanf.GetXaxis()->SetLabelOffset(0.020);
  spmtlbmeanf.GetXaxis()->SetLabelSize(0.03);
  spmtlbmeanf.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) / FADC");
  spmtlbmeanf.Draw("P");

  ssdpmtlbmeanf.GetXaxis()->SetTimeDisplay(1);
  ssdpmtlbmeanf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  ssdpmtlbmeanf.GetXaxis()->SetLabelOffset(0.020);
  ssdpmtlbmeanf.GetXaxis()->SetLabelSize(0.03);
  ssdpmtlbmeanf.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) Low / FADC");
  ssdpmtlbmeanf.Draw("P");

  ssdpmthbmeanl.GetXaxis()->SetTimeDisplay(1);
  ssdpmthbmeanl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  ssdpmthbmeanl.GetXaxis()->SetLabelOffset(0.020);
  ssdpmthbmeanl.GetXaxis()->SetLabelSize(0.03);
  ssdpmthbmeanl.GetYaxis()->SetTitle("Mean Baseline (First 100 bins) High / FADC");
  ssdpmthbmeanl.Draw("P");
  
  pmt1lbrmsf.GetXaxis()->SetTimeDisplay(1);
  pmt1lbrmsf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt1lbrmsf.GetXaxis()->SetLabelOffset(0.020);
  pmt1lbrmsf.GetXaxis()->SetLabelSize(0.03);
  pmt1lbrmsf.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) Low / FADC");
  pmt1lbrmsf.Draw("P");

  pmt1hbrmsl.GetXaxis()->SetTimeDisplay(1);
  pmt1hbrmsl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt1hbrmsl.GetXaxis()->SetLabelOffset(0.020);
  pmt1hbrmsl.GetXaxis()->SetLabelSize(0.03);
  pmt1hbrmsl.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) High / FADC");
  pmt1hbrmsl.Draw("P");

  pmt2lbrmsf.GetXaxis()->SetTimeDisplay(1);
  pmt2lbrmsf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt2lbrmsf.GetXaxis()->SetLabelOffset(0.020);
  pmt2lbrmsf.GetXaxis()->SetLabelSize(0.03);
  pmt2lbrmsf.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) Low / FADC");
  pmt2lbrmsf.Draw("P");

  pmt2hbrmsl.GetXaxis()->SetTimeDisplay(1);
  pmt2hbrmsl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt2hbrmsl.GetXaxis()->SetLabelOffset(0.020);
  pmt2hbrmsl.GetXaxis()->SetLabelSize(0.03);
  pmt2hbrmsl.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) High / FADC");
  pmt2hbrmsl.Draw("P");

  pmt3lbrmsf.GetXaxis()->SetTimeDisplay(1);
  pmt3lbrmsf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt3lbrmsf.GetXaxis()->SetLabelOffset(0.020);
  pmt3lbrmsf.GetXaxis()->SetLabelSize(0.03);
  pmt3lbrmsf.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) Low / FADC");
  pmt3lbrmsf.Draw("P");

  pmt3hbrmsl.GetXaxis()->SetTimeDisplay(1);
  pmt3hbrmsl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  pmt3hbrmsl.GetXaxis()->SetLabelOffset(0.020);
  pmt3hbrmsl.GetXaxis()->SetLabelSize(0.03);
  pmt3hbrmsl.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) High / FADC");
  pmt3hbrmsl.Draw("P");

  spmtlbrmsf.GetXaxis()->SetTimeDisplay(1);
  spmtlbrmsf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  spmtlbrmsf.GetXaxis()->SetLabelOffset(0.020);
  spmtlbrmsf.GetXaxis()->SetLabelSize(0.03);
  spmtlbrmsf.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) / FADC");
  spmtlbrmsf.Draw("P");

  ssdpmtlbrmsf.GetXaxis()->SetTimeDisplay(1);
  ssdpmtlbrmsf.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  ssdpmtlbrmsf.GetXaxis()->SetLabelOffset(0.020);
  ssdpmtlbrmsf.GetXaxis()->SetLabelSize(0.03);
  ssdpmtlbrmsf.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) Low / FADC");
  ssdpmtlbrmsf.Draw("P");

  ssdpmthbrmsl.GetXaxis()->SetTimeDisplay(1);
  ssdpmthbrmsl.GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
  ssdpmthbrmsl.GetXaxis()->SetLabelOffset(0.020);
  ssdpmthbrmsl.GetXaxis()->SetLabelSize(0.03);
  ssdpmthbrmsl.GetYaxis()->SetTitle("Baseline RMS (First 100 bins) High / FADC");
  ssdpmthbrmsl.Draw("P");

  //c1.Print("../plots/jan2021/bl"+nameStati+".pdf");

  hfile.Write();
  hfile.Close();

  return 0;
}
