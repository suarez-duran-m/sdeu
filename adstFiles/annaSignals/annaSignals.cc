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

#include <fstream>

using namespace std;

/********************/
/* Global variables */
/********************/

void readPlotPerSt(unsigned int stid, TString stName, int argc, char** argv,
    const int gps1stFeb2021, TString printPath, vector<vector<double>> &beta,
    vector<vector<double>> &errBeta);


int main ( int argc, char** argv) {
  if ( argc < 3 ) { cout << endl
      << "Usage: " << argv[0] << " <file_St_for_analysis> <ifControl> <file_adstfiles>" << endl
        << "<file_St_for_analysis>: File with station IDs for analysis" << endl
        << "<ifControl>: 1 if you running for control porpoise -> output in results/control/ " << endl
        << "             0 otherwise" << endl  
        << "<file_adstfiles>: File with list of adst files to read" 
        << endl;
    exit(0);
  }

  // Station for analysis
  const char *filename = argv[1];
  ifstream fileWithIds(filename, ios::in);
  if (!fileWithIds.is_open()) {
    cout << "Could not open file: " << filename << endl;
    exit(0);
  }
  vector < double > st2read;
  vector < double > graphXerrors;
  vector < double > graphXbins;
  bool ifControl = atoi(argv[2]);
  unsigned int id = 0;
  unsigned int binx = 1;
  while ( fileWithIds.good() ) {
    id = 0;
    fileWithIds >> id;
    if (id) {
      st2read.push_back( id );
      graphXerrors.push_back(0.);
      graphXbins.push_back(binx);
      binx++;
    }
  }
  TString stName; 
  TString printPath = ifControl ? "results/control/" : "results/";

  // Dates to plot Total signal distributions
  const int gps1stDec2021 = ifControl ? 1196121618 : 1322352018;
  const int gps1stFeb2022 = 1327708818;

  // betTotSgnl: [0] for Bef, [1] for Aft
  vector < vector < double > > betTotSgnl(2);
  vector < vector < double > > errBetTotSgnl(2);

  // Extracting signals per Station
  for ( auto & st_i : st2read ) {
    stName.Form("%.f", st_i);
    readPlotPerSt(st_i, stName, argc, argv, gps1stDec2021, printPath,
        betTotSgnl, errBetTotSgnl);
  }

  TFile betasGraph ("results/betasGraph.root", "RECREATE");

  TH1D grpBef("grpBef","", st2read.size(), 0, st2read.size());
      //st2read.size(), &graphXbins[0], &betTotSgnl[0][0], &graphXerrors[0], &errBetTotSgnl[0][0]);
  TH1D grpAft("grpAft","", st2read.size(), 0, st2read.size());
      //st2read.size(), &graphXbins[0], &betTotSgnl[1][0], &graphXerrors[0], &errBetTotSgnl[1][0]);

 
  for ( int i=0; i<st2read.size(); i++ ) {
    stName.Form("%.f", st2read[i]);
    grpBef.SetBinContent(i+1, betTotSgnl[0][i]);
    grpBef.SetBinError(i+1, errBetTotSgnl[0][i]);
    grpBef.GetXaxis()->SetBinLabel(i+1, stName); 

    grpAft.SetBinContent(i+1, betTotSgnl[1][i]);
    grpAft.SetBinError(i+1, errBetTotSgnl[1][i]);
    grpAft.GetXaxis()->SetBinLabel(i+1, stName);
  } 

  TCanvas c("c", "c", 102, 76, 1600, 900);
  c.SetBorderMode(0);
  c.SetBorderSize(2);
  c.SetRightMargin(0.017);
  c.SetLeftMargin(0.074);
  c.SetTopMargin(0.014);
  c.SetBottomMargin(0.1);
  c.SetFrameBorderMode(0);
  c.SetFrameBorderMode(0); 
  c.cd();

  grpBef.SetStats(kFALSE);
  grpBef.GetXaxis()->SetTitle("St. Id");
  grpBef.GetYaxis()->SetTitle("[au]");
  grpBef.GetYaxis()->SetRangeUser(-5, 5);
  grpBef.SetLineColor(kBlue);
  grpBef.SetMarkerStyle(20);
  grpBef.Draw("E1");

  grpAft.SetLineColor(kRed);
  grpAft.SetMarkerColor(kRed);
  grpAft.SetMarkerStyle(24);
  grpAft.Draw("E1 same");

  TLegend legend(0.6, 0.5, 0.98, 0.95);
  legend.AddEntry(&grpBef, "Before","l");
  legend.AddEntry(&grpAft, "After","l");
  legend.SetTextSize(0.04);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.Draw();

  c.Print("results/betasGraph.pdf");

  grpBef.Write();
  grpAft.Write();
  c.Write();

  betasGraph.Write();
  betasGraph.Close();


  return 0;
}


void readPlotPerSt(unsigned int stid, TString stName, int argc, char** argv,
    const int gps1stDec2021,TString printPath, vector<vector<double>> &beta, 
    vector<vector<double>> &errBeta) {
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
    dataFile.readSignals(stid);

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

  vector<double> tmpBetas;
  tmpBetas = pltDist.getBetaFitTotSgnl();
  
  beta[0].push_back( tmpBetas[0] );
  beta[1].push_back( tmpBetas[1] );
  errBeta[0].push_back( tmpBetas[2] );
  errBeta[1].push_back( tmpBetas[3] );  

  pltDist.writeRootFile();
}
