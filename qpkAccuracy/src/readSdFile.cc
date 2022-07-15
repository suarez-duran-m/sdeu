#include "readSdFile.h"

#include "TFile.h"
#include "TH1.h"
#include "TString.h"
#include "TSystem.h"

#include "fitcharge.h"

using namespace std;

readSdFile::readSdFile():
  totalEvents(0),
  binsLeftRigh(0),
  tmpTestQpk(0.)
{}

readSdFile::readSdFile(int nfiles, char * fileNames[], int binLR) {
  sdFile = new AugerIoSd(nfiles, fileNames);
  totalEvents = sdFile->NumberOfEvents();
  binsLeftRigh = binLR;
}


void readSdFile::GetAndWriteChHisto(int StId, char *date) {
  TString outName(Form("results/chargeHistosSt%d_%s.root", StId, date));
  TFile outFile(outName,"RECREATE");

  // Store histograms
  TTree peakTree("peakHistos", "");
  // Store fitter results
  TTree *pmt1Tree = new TTree("pmt1","");
  TTree *pmt2Tree = new TTree("pmt2","");
  TTree *pmt3Tree = new TTree("pmt3","");

  chPmt1 = new TGraphErrors();
  chPmt2 = new TGraphErrors();
  chPmt3 = new TGraphErrors();
  TH1F *pkPmt1 = new TH1F();
  TH1F *pkPmt2 = new TH1F();
  TH1F *pkPmt3 = new TH1F();
  int evtId = 0;
  int utcTime = 0;
  
  for(int i=0; i<10; i++)
   fitVar[i] = 0.;

  SetTree(pmt1Tree);
  SetTree(pmt2Tree);
  SetTree(pmt3Tree);

  peakTree.Branch("peakPmt1","TH1F",&pkPmt1,32000,0);
  peakTree.Branch("peakPmt2","TH1F",&pkPmt2,32000,0);
  peakTree.Branch("peakPmt3","TH1F",&pkPmt3,32000,0);
  peakTree.Branch("eventId", &evtId,"evtId/I");
  
  EventPos pos;
  int nrEventsRead = 0;
  int previusEvent = 0;
  bool found = false;
  fitcharge *fitChHisto;

  for (pos = sdFile->FirstEvent(); pos<sdFile->LastEvent(); pos=sdFile->NextEvent()) {
    nrEventsRead++;
    if (nrEventsRead%1000 == 0)
       cout << "====> Read " << nrEventsRead << " out of " << totalEvents << endl;
    IoSdEvent event(pos);
    if ( event.Id == previusEvent )
      continue;
    previusEvent = event.Id;

    for ( unsigned int st_i = 0 ; st_i < event.Stations.size(); st_i++ ) {
      found = false;
      if ( !event.Stations[st_i].IsUUB && event.Stations[st_i].Error != 256 )
        continue;
      if ( StId == event.Stations[st_i].Id )
        found = true;    
      if ( !found ) 
        continue;
 
      evtId = event.Id;
      utcTime = event.UTCTime;
      IdAndTime[0] = event.Id;
      IdAndTime[1] = event.UTCTime;

      chPmt1 = FillingTree(event.Stations[st_i].HCharge(0),
          event.Stations[st_i].Histo->Offset[0+6], binsLeftRigh);
      pmt1Tree->Fill();

      chPmt2 = FillingTree(event.Stations[st_i].HCharge(1),
          event.Stations[st_i].Histo->Offset[1+6], binsLeftRigh);
      pmt2Tree->Fill();

      chPmt3 = FillingTree(event.Stations[st_i].HCharge(2),
          event.Stations[st_i].Histo->Offset[2+6], binsLeftRigh);
      pmt3Tree->Fill();
      
      pkPmt1 = event.Stations[st_i].HPeak(0);
      pkPmt2 = event.Stations[st_i].HPeak(1);
      pkPmt3 = event.Stations[st_i].HPeak(2);
      peakTree.Fill();
    }
  }
  outFile.Write();
  outFile.Close();
}


void readSdFile::SetTree(TTree *pmtTree) {  
  pmtTree->Branch("ChargeHisto","TGraphErrors", &chPmt1, 32000,0);
  pmtTree->Branch("eventId", &IdAndTime[0],"IdAndTime[0]/I");
  pmtTree->Branch("eventUTC", &IdAndTime[1],"IdAndTime[1]/I");
  pmtTree->Branch("Qpk", &fitVar[0],"fitQpk/D");
  pmtTree->Branch("QpkDer", &fitVar[1],"fitQpkDer/D");
  pmtTree->Branch("rangXmin", &fitVar[2],"fitRangXmin/D");
  pmtTree->Branch("rangXmax", &fitVar[3],"fitRangXmax/D");
  pmtTree->Branch("par0", &fitVar[4],"fitPar0/D");
  pmtTree->Branch("par1", &fitVar[5],"fitPar1/D");
  pmtTree->Branch("par2", &fitVar[6],"fitPar2/D");
  pmtTree->Branch("ndf", &fitVar[7],"fitNdf/D");
  pmtTree->Branch("chi2", &fitVar[8],"fitChi2/D");
  pmtTree->Branch("prob", &fitVar[9],"fitProb/D");
}

TGraphErrors *readSdFile::FillingTree(TH1F *hist, int Offset, int binsLeftRigh) {    
  TGraphErrors *retGraph;
  fitcharge *fitChHisto;
  fitChHisto = new fitcharge (hist, Offset, binsLeftRigh);
  fitVar[0] = fitChHisto->qpkPos;
  fitVar[1] = fitChHisto->qpkPosDeri;
  fitVar[2] = fitChHisto->rangXmin;
  fitVar[3] = fitChHisto->rangXmax;
  fitVar[4] = fitChHisto->par0;
  fitVar[5] = fitChHisto->par1;
  fitVar[6] = fitChHisto->par2;
  fitVar[7] = fitChHisto->ndf;
  fitVar[8] = fitChHisto->chi2;
  fitVar[9] = fitChHisto->prob;
  tmpTestQpk = fitChHisto->qpkPos;
  retGraph = (TGraphErrors*)fitChHisto->GetFitGraph()->Clone();
  fitChHisto->DeleteObjts();
  return retGraph;
}

const vector<vector<double>> readSdFile::GetQpk(int StId, char *date) {
  vector < vector < double > > retQpk(3);
  TString fileName(Form("results/chargeHistosSt%d_%s.root", StId, date));
  TSystem gsystem;
  if ( gsystem.AccessPathName(fileName) ) {
    cerr << "File " << fileName << " does not exits" << endl;    
    exit(0);
  }

  TFile *inFile = TFile::Open(fileName); 
  TTree *trCh = (TTree*)inFile->Get("chargeHistos");
  double qpkPmt1 = 0.;
  double qpkPmt2 = 0.;
  double qpkPmt3 = 0.;

  trCh->SetBranchAddress("QpkPmt1", &qpkPmt1);
  trCh->SetBranchAddress("QpkPmt2", &qpkPmt2);
  trCh->SetBranchAddress("QpkPmt3", &qpkPmt3);

  for ( int i=0; i<trCh->GetEntries(); i++ ) {
    trCh->GetEntry(i);
    retQpk[i].push_back( qpkPmt1 );
  }

  return retQpk;
}
