#include <RecEvent.h>
#include <RecEventFile.h>
#include <DetectorGeometry.h>

#include <TFile.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>

TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 102, 76, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetRightMargin(0.04);
  canvas->SetLeftMargin(0.13);
  canvas->SetTopMargin(0.014);
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

using namespace std;

/********************/
/* Global variables */
/********************/

int main ( int argc, char** argv) {
  if ( argc < 2 ) {
    cout << endl
      << "Usage: " << argv[0] << " <file_adstfiles> <output> " << endl
      << "<file_adstfiles>: File with list of adst files to read" << endl
      << "<output>: name for adst_output file" << endl;
    exit(0);
  }
  // Vectors to store signals values
  vector < double > sig;
  vector < double > sigErr;
  // Reading ADST files
  for ( int adst_i = 1; adst_i<argc-1; adst_i++ ) {
    RecEventFile inputFile(argv[adst_i]);    
    RecEvent *theRecEvent = new RecEvent();
    inputFile.SetBuffers(&theRecEvent);
    cout << endl << endl << "Opened " << argv[adst_i] 
      << " with " << inputFile.GetNEvents() 
      << " events " << endl << endl;
    // Defining the Geometry to plot the Auger Map
    DetectorGeometry theGeometry;
    inputFile.ReadDetectorGeometry(theGeometry);
    // Reading events
    while ( inputFile.ReadNextEvent() == RecEventFile::eSuccess ) {
      const auto &sdEvent = theRecEvent->GetSDEvent();
      int nCand = sdEvent.GetNumberOfCandidates();
      for ( int i=0; i<nCand; i++ ) {
        sig.push_back( sdEvent.GetStation(i)->GetTotalSignal() );
        sigErr.push_back( sdEvent.GetStation(i)->GetTotalSignalError() );
      }
    }
  }
  // Plotting Signal distribution
  TH1D *histSig = new TH1D("histSig", "", 2000, 0, 2e3);
  for ( auto & sig_i : sig ) 
    histSig->Fill( sig_i );

  TCanvas *sigDistCanv = canvasStyle("sigDistCanv");
  sigDistCanv->cd();
  sigDistCanv->SetLogy(1);
  sigDistCanv->SetLogx(1);
  
  histSig->GetYaxis()->SetTitle("Counts [au]");
  histSig->GetXaxis()->SetTitle("Log10(S/VEM)");
  histSig->GetXaxis()->SetTitleOffset(1.3);
  histSig->Draw();
  sigDistCanv->Print("kk.pdf");
  
  return 1;
}
