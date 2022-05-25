#include <DetectorGeometry.h>
#include <Traces.h>
#include <TraceType.h>

#include <TFile.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>

#include "include/readAdstFile.h"
#include "include/plotDiffDist.h"

using namespace std;

/********************/
/* Global variables */
/********************/


int main ( int argc, char** argv) {
  if ( argc < 3 ) { cout << endl
      << "Usage: " << argv[0] << " <st_for_analysis> <ifControl> <file_adstfiles>" << endl
        << "<st_for_analysis>: Station ID for analysis" << endl
        << "<ifControl>: True if you running for control porpoise -> output in results/control/" << endl
        << "<file_adstfiles>: File with list of adst files to read" 
        << endl;
    exit(0);
  }

  // Station for analysis
  int st2read = atoi(argv[1]);
  bool ifControl = argv[2];
  TString stName;
  stName.Form("%d", st2read);
  TString printPath = ifControl ? "results/control/" : "results/";

  // Dates to plot Total signal distributions
  const int gps1stDec2021 = ifControl ? 1196121618 : 1322352018;
  const int gps1stFeb2022 = 1327708818;  
  //vector < int > st2read = {846, 860, 1190, 1191, 1220, 1779};

  // Histograms for signals differences distribution
  double frstBin = -10.;
  double lstBin = 10.;
  int nBins = 200;

  TH1D dist12("dist12","",nBins, frstBin, lstBin);
  TH1D dist13("dist13","",nBins, frstBin, lstBin);
  TH1D dist23("dist23","",nBins, frstBin, lstBin);

  // Histograms for Muon signal decay, tau
  frstBin = 25.;
  lstBin = 125.;
  nBins = 50;
  TH1D dstMuonTau2021Pmt1("dstMuonTau2021Pmt1","", nBins, frstBin, lstBin);
  TH1D dstMuonTau2021Pmt2("dstMuonTau2021Pmt2","", nBins, frstBin, lstBin);
  TH1D dstMuonTau2021Pmt3("dstMuonTau2021Pmt3","", nBins, frstBin, lstBin);

  TH1D dstMuonTau2022Pmt1("dstMuonTau2022Pmt1","", nBins, frstBin, lstBin);
  TH1D dstMuonTau2022Pmt2("dstMuonTau2022Pmt2","", nBins, frstBin, lstBin);
  TH1D dstMuonTau2022Pmt3("dstMuonTau2022Pmt3","", nBins, frstBin, lstBin);

  // Histograms for Total signals
  frstBin = 0.;
  lstBin = 20.;
  nBins = 20;
  TH1D totSglDist2021("totSglDist2021","", nBins, frstBin, lstBin);
  TH1D totSglDist2022("totSglDist2022","", nBins, frstBin, lstBin);

  // Reading ADST files and extracting signals for each PMT
  for (int file_i=3; file_i<argc; file_i++) {
    readAdstFile dataFile;
    dataFile.setAdstFileName(argv[file_i]);
    //dataFile.checkNumberEvts();
    dataFile.readSignals(st2read);

    const vector < vector < double > > signalsPmt = dataFile.getSglPmt();
    const vector < vector < double > > muonTauPmt = dataFile.getMuonTauPmt();
    const vector < double > totalSignal = dataFile.getTotSgl(); 
    const vector < double > totalSignalErr = dataFile.getTotSglErr();
    const vector < double > signalGps = dataFile.getSglGpsSec();

    for(int i=0; i<signalsPmt[0].size(); i++) {
      dist12.Fill( (signalsPmt[0][i] - signalsPmt[1][i]) / totalSignalErr[i] );
      dist13.Fill( (signalsPmt[0][i] - signalsPmt[2][i]) / totalSignalErr[i] );
      dist23.Fill( (signalsPmt[1][i] - signalsPmt[2][i]) / totalSignalErr[i] );
    }

    for(int i=0; i<signalGps.size(); i++) {
      if ( signalGps[i] < gps1stDec2021 ) { 
        totSglDist2021.Fill( totalSignal[i] );
        dstMuonTau2021Pmt1.Fill( muonTauPmt[0][i] );
        dstMuonTau2021Pmt2.Fill( muonTauPmt[1][i] );
        dstMuonTau2021Pmt3.Fill( muonTauPmt[2][i] );
      }
      else {
        totSglDist2022.Fill( totalSignal[i] );
        dstMuonTau2022Pmt1.Fill( muonTauPmt[0][i] );
        dstMuonTau2022Pmt2.Fill( muonTauPmt[1][i] );
        dstMuonTau2022Pmt3.Fill( muonTauPmt[2][i] );
      }
    }
  }

  plotDiffDist pltDist(stName, printPath, &dist12, &dist13, &dist23, 
      &totSglDist2021, &totSglDist2022);
  
  pltDist.plotTausMuon(&dstMuonTau2021Pmt1, &dstMuonTau2021Pmt2, &dstMuonTau2021Pmt3,
      &dstMuonTau2022Pmt1, &dstMuonTau2022Pmt2, &dstMuonTau2022Pmt3);
  pltDist.writeRootFile();

  return 0;
}
