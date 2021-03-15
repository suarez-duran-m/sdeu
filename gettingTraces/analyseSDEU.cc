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


int main (int argc, char *argv[]) {
  if ( argc < 4 ) {
    cout << endl << "Usage: " << argv[0] << " <stationsFile>  <output>  <files>" << endl
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
  
  /*cout << "======> I am  selecting events containing stations: "; 
  for (  vector<unsigned int>::const_iterator iter= stationsIds.begin();
         iter!= stationsIds.end(); ++iter)
    cout << *iter << " ";
  cout << " <====== " << endl;
  */
  
  IoSd *sdFile = NULL;
  AugerFile *adFile = NULL;

  TH1F pmt1l("pmt1l","Traces for PMT1 Low",2047,0,2048);
  TH1F pmt2l("pmt2l","Traces for PMT2 Low",2047,0,2048);
  TH1F pmt3l("pmt3l","Traces for PMT3 Low",2047,0,2048);
  TH1F spmtl("spmtl","Traces for SPMT Low",2047,0,2048);
  TH1F ssdpmtl("ssdpmtl","Traces for SSDPMT Low",2047,0,2048);

  TH1F pmt1h("pmt1h","Traces for PMT1 High",2047,0,2048);
  TH1F pmt2h("pmt2h","Traces for PMT2 High",2047,0,2048);
  TH1F pmt3h("pmt3h","Traces for PMT3 High",2047,0,2048);
  TH1F spmth("spmth","Traces for SPMT High",2047,0,2048);
  TH1F ssdpmth("ssdpmth","Traces for SSDPMT High",2047,0,2048);

  TCanvas c1("c1","test",0,0,3600,2400);
  c1.Divide(2,5);

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
      if ( event.Stations[i].IsUUB ) {
        cout << "# Event " << event.Id << " Station " << event.Stations[i].Id
             << endl;
        //cout << "# Error " << event.Stations[i].Error << endl;        

        IoSdEvent event(pos);
        if (event.RootClassName == "AugerEvent") {
          if (!adFile)
            adFile = new AugerFile(outputName, AugerFile::eWrite);
          adFile->Write(*(event.RawAugerEvent()), false);
        }     
        else { 
          if (!sdFile)
            sdFile = new IoSd(outputName, "w");
          sdFile->Write(event, "");
        }

        if (event.Stations[i].Error==256) { //0+256
          //cout << "# sizeTraces: " 
               //cout << event.Stations[i].UFadc->Traces.size() << endl;
          //for (int k=0;k<event.Stations[i].UFadc->NSample;k++) {
          for (int k=0;k<event.Stations[i].UFadc->NSample;k++) {
            //cout << k << endl;

            pmt1l.Fill( k, event.Stations[i].UFadc->GetValue(0,1,k) );
            pmt2l.Fill( k, event.Stations[i].UFadc->GetValue(1,1,k) ); 
            pmt3l.Fill( k, event.Stations[i].UFadc->GetValue(2,1,k) );
            spmtl.Fill( k, event.Stations[i].UFadc->GetValue(3,1,k) );
            ssdpmtl.Fill( k, event.Stations[i].UFadc->GetValue(4,1,k) );

            pmt1h.Fill( k, event.Stations[i].UFadc->GetValue(0,0,k) );
            pmt2h.Fill( k, event.Stations[i].UFadc->GetValue(1,0,k) );
            pmt3h.Fill( k, event.Stations[i].UFadc->GetValue(2,0,k) );
            spmth.Fill( k, event.Stations[i].UFadc->GetValue(3,1,k) );
            ssdpmth.Fill( k, event.Stations[i].UFadc->GetValue(4,1,k) );

            //for (int j=0;j<5;j++) // Run on type of PMT
              //cout << " " << event.Stations[i].UFadc->GetValue(j,0,k) << " " // 0 for high; 1 for low
              //     << event.Stations[i].UFadc->GetValue(j,1,k);
            //cout << endl;
          }
        }
            
        c1.cd(1);
        pmt1l.Draw();
        c1.cd(2);
        pmt1h.Draw();
        c1.cd(3);
        pmt2l.Draw();
        c1.cd(4);
        pmt2h.Draw();
        c1.cd(5);
        pmt3l.Draw();
        c1.cd(6);
        pmt3h.Draw();
        c1.cd(7);
        spmtl.Draw();
        c1.cd(8);
        spmth.Draw();
        c1.cd(9);
        ssdpmtl.Draw();
        c1.cd(10);
        ssdpmth.Draw();
        
        TString nameEvent = TString::UItoa(event.Id, 10);
        TString nameStati = TString::UItoa(event.Stations[i].Id, 10);

        //c1.Print("plots/jan2021/"+nameEvent+"-"+nameStati+".pdf");
        
        pmt1l.Reset();
        pmt2l.Reset();
        pmt3l.Reset();
        spmtl.Reset();
        ssdpmtl.Reset();

        pmt1h.Reset();
        pmt2h.Reset();
        pmt3h.Reset();
        spmth.Reset();
        ssdpmth.Reset();
       }
     }
  }
  return 0;
}
