#include <RecEvent.h>
#include <RecEventFile.h>
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

#include <vector>

using namespace std;

void setCanvasStyle(TCanvas& canvas);
void setTH1FontStyle(TH1D& histo);
TGraphErrors doHisto(vector<int> counts, vector<double> bins, TString stId, 
    TString typeOutlier);
void setTGraphFontStyle(TGraph& grph);
void plottingSt666(TGraph grph, TString typeHist, TString pmtId);


int main ( int argc, char** argv) {
  if ( argc < 2 ) { cout << endl
      << "Usage: " << argv[0] << " <file_adstfiles> <output> " << endl
      << "<file_adstfiles>: File with list of adst files to read" << endl
      << "<output>: name for adst_output file" << endl;
    exit(0);
  }

  TH1D distInteCQ("distInteCQ","", 250 , 0, 6);
  TH1D distInteCI("distInteCI","", 250 , 0, 6);

  vector < vector < double > > countsCI_st666(3);
  vector < vector < double > > countsCQ_st666(3);
  vector < double > st666_binsCI(150);
  vector < double > st666_binsCQ(600);
  vector < int > tot666Histos(3);
  for ( int i=0; i<3; i++ ) {
    countsCI_st666[i].resize(150);
    countsCQ_st666[i].resize(600);
    tot666Histos[i] = 0;
  }
 
  TCanvas cvnsOutlierCI4e4("cvnsOutlierCI4e4", "", 1.6e3, 9e2);
  cvnsOutlierCI4e4.Print("results/outliersCI4e4.pdf(");
  TCanvas cvnsOutlierCI1e4("cvnsOutlierCI1e4", "", 1.6e3, 9e2);
  cvnsOutlierCI1e4.Print("results/outliersCI1e4.pdf(");

  TCanvas cvnsOutlierCQ1e4("cvnsOutlierCQ1e4", "", 1.6e3, 9e2);
  cvnsOutlierCQ1e4.Print("results/outliersCQ1e4.pdf(");
  TCanvas cvnsOutlierCQ1p5e4("cvnsOutlierCQ1p5e4", "", 1.6e3, 9e2);
  cvnsOutlierCQ1p5e4.Print("results/outliersCQ1p5e4.pdf(");

  gErrorIgnoreLevel = kWarning;

  // Reading ADST files
  for ( int adst_i = 1; adst_i<argc-1; adst_i++ ) {
    RecEventFile inputFile(argv[adst_i]);
    RecEvent *theRecEvent = new RecEvent();
    inputFile.SetBuffers(&theRecEvent);

    // Reading events
    //
    double signal = 0.;
    while ( inputFile.ReadNextEvent() == RecEventFile::eSuccess ) {
      const auto &sdEvent = theRecEvent->GetSDEvent();
      const std::vector<SdRecStation>& stations = sdEvent.GetStationVector();
      // Reading all station for the event
      //
      for (unsigned int iS = 0; iS < stations.size(); ++iS) {
        const SdRecStation& recStation = stations[iS];
        if ( recStation.IsUUB() ) {
          const vector<Traces>& traces = recStation.GetPMTTraces();
          
          // Reading per PMT
          //
          for (vector<Traces>::const_iterator trIter = traces.begin();
              trIter != traces.end(); ++trIter) {
            const Traces& tr = *trIter;
            cout << "MSD eTotalTrace: " << eTotalTrace << endl;
            if (tr.GetType() != eTotalTrace)
              continue;
            const unsigned int pmtId = tr.GetPMTId();
            if ( pmtId > 3 )
              continue;
            if ( tr.IsVEMChargeFromHistogram() )
              cout << "MSD from Fit " << endl;
            else 
              cout << "MSD from online" << endl;

            // Correction by big-bins
            //
            double bigbins = 4.0;
            int integral = 0;
            const vector<Int_t>& valHistoCI = tr.GetCoinciPeakHistogram().GetValues();
            if ( !valHistoCI.size() )
              continue;

            if ( tr.GetCoinciCharge() > 2300 && tr.GetCoinciCharge() < 2600 )
              cout << "MSD 2500 " << theRecEvent->GetSDEvent().GetEventId() << endl;

            for ( int j = 0; j < 103; j++ ) {
              integral += valHistoCI[j];
              if ( recStation.GetId() == 666 )
                countsCI_st666[pmtId-1][j] += valHistoCI[j];              
            }
            for (int j = 0; j < 47; j++) {
              integral += valHistoCI[j + 103] / bigbins;
              if ( recStation.GetId() == 666 )
                countsCI_st666[pmtId-1][j+103] += valHistoCI[j+103] / bigbins;
            }            
            distInteCI.Fill ( log10(integral) );
            if ( recStation.GetId() == 666 )
              tot666Histos[pmtId-1]++;

            if ( integral > 4e4 ) {
              const vector< double >& bins = tr.GetCoinciPeakHistogram().GetBinning();
              if ( recStation.GetId() == 666 )
                for ( int i=0; i<bins.size(); i++ )
                st666_binsCI[i] = bins[i];

              TString stId = Form("%d", recStation.GetId());
              TString typeOutlier = "CI";
              cvnsOutlierCI4e4.cd();
              TGraphErrors graphHisto = doHisto( valHistoCI, bins, stId, typeOutlier );
              graphHisto.Draw("AL");
              TLegend lgnd(0.6, 0.6, 0.88, 0.88);
              lgnd.AddEntry(&graphHisto, Form("Evt: %d", 
                    theRecEvent->GetSDEvent().GetEventId()), "");
              lgnd.AddEntry(&graphHisto, "Station: "+stId, "");
              lgnd.AddEntry(&graphHisto, Form("PMT: %d", tr.GetPMTId()), "");
              lgnd.AddEntry(&graphHisto, Form("Entries: %d", integral), "");
              lgnd.SetBorderSize(0);
              lgnd.SetTextSize(0.055);
              lgnd.Draw();
              cvnsOutlierCI4e4.Print("results/outliersCI4e4.pdf");
            }
            if ( integral < 1e4 ) {
              const vector< double >& bins = tr.GetCoinciPeakHistogram().GetBinning(); 
              TString stId = Form("%d", recStation.GetId());
              TString typeOutlier = "CI";
              cvnsOutlierCI1e4.cd();
              TGraphErrors graphHisto = doHisto( valHistoCI, bins, stId, typeOutlier );
              graphHisto.Draw("AL");
              TLegend lgnd(0.6, 0.6, 0.88, 0.88);
              lgnd.AddEntry(&graphHisto, Form("Evt: %d", 
                    theRecEvent->GetSDEvent().GetEventId()), "");
              lgnd.AddEntry(&graphHisto, "Station: "+stId, "");
              lgnd.AddEntry(&graphHisto, Form("PMT: %d", tr.GetPMTId()), "");
              lgnd.AddEntry(&graphHisto, Form("Entries: %d", integral), "");
              lgnd.SetBorderSize(0);
              lgnd.SetTextSize(0.055);
              lgnd.Draw();
              cvnsOutlierCI1e4.Print("results/outliersCI1e4.pdf"); 
            }

            integral = 0;
            const vector<Int_t>& valHistoCQ = tr.GetCoinciChargeHistogram().GetValues();
            for ( int j = 0; j < 403; j++ ) {
              integral += valHistoCQ[j];
              if ( recStation.GetId() == 666 )
                countsCQ_st666[pmtId-1][j] += valHistoCQ[j];
            }
            for (int j = 0; j < 197; j++) {
              integral += valHistoCQ[j + 403] / bigbins;
              if ( recStation.GetId() == 666 )
                countsCQ_st666[pmtId-1][j+403] += valHistoCQ[j+403] / bigbins;
            }
            distInteCQ.Fill ( log10(integral) );

            if ( log10(integral) < 3.9 ) {
              const vector< double >& bins 
                = tr.GetCoinciChargeHistogram().GetBinning();
              TString stId = Form("%d", recStation.GetId());
              TString typeOutlier = "CQ";
              cvnsOutlierCQ1e4.cd();
              TGraphErrors graphHisto = doHisto( valHistoCQ, bins, stId, typeOutlier );
              graphHisto.Draw("AL");
              TLegend lgnd(0.6, 0.6, 0.88, 0.88);
              lgnd.AddEntry(&graphHisto, Form("Evt: %d", 
                    theRecEvent->GetSDEvent().GetEventId()), "");
              lgnd.AddEntry(&graphHisto, "Station: "+stId, "");
              lgnd.AddEntry(&graphHisto, Form("PMT: %d", tr.GetPMTId()), "");
              lgnd.AddEntry(&graphHisto, Form("Entries: %d", integral), "");
              lgnd.SetBorderSize(0);
              lgnd.SetTextSize(0.055);
              lgnd.Draw();
              cvnsOutlierCQ1e4.Print("results/outliersCQ1e4.pdf"); 
            }
            if ( log10(integral) > 3.9 && log10(integral) < 4.2 ) {
              const vector< double >& bins 
                = tr.GetCoinciChargeHistogram().GetBinning();
              if ( recStation.GetId() == 666 )
                for ( int i=0; i<bins.size(); i++ )
                  st666_binsCQ[i] = bins[i];

              TString stId = Form("%d", recStation.GetId());
              TString typeOutlier = "CQ";
              cvnsOutlierCQ1p5e4.cd();
              TGraphErrors graphHisto = doHisto( valHistoCQ, bins, stId, typeOutlier );
              graphHisto.Draw("AL");
              TLegend lgnd(0.6, 0.6, 0.88, 0.88);
              lgnd.AddEntry(&graphHisto, Form("Evt: %d", 
                    theRecEvent->GetSDEvent().GetEventId()), "");
              lgnd.AddEntry(&graphHisto, "Station: "+stId, "");
              lgnd.AddEntry(&graphHisto, Form("PMT: %d", tr.GetPMTId()), "");
              lgnd.AddEntry(&graphHisto, Form("Entries: %d", integral), "");
              lgnd.SetBorderSize(0);
              lgnd.SetTextSize(0.055);
              lgnd.Draw();
              cvnsOutlierCQ1p5e4.Print("results/outliersCQ1p5e4.pdf"); 
            }
          }
        }
      }
    }
  }

  gErrorIgnoreLevel = kPrint;
  cvnsOutlierCI4e4.Print("results/outliersCI4e4.pdf)"); 
  cvnsOutlierCI1e4.Print("results/outliersCI1e4.pdf)"); 
  cvnsOutlierCQ1e4.Print("results/outliersCQ1e4.pdf)");
  cvnsOutlierCQ1p5e4.Print("results/outliersCQ1p5e4.pdf)");

  TCanvas cvnsInteCQ("cvnsInteCQ", "", 1.6e3, 9e2);
  setCanvasStyle(cvnsInteCQ);
  cvnsInteCQ.cd();
  cvnsInteCQ.SetLogy();
  
  distInteCQ.SetStats(kFALSE);
  distInteCQ.Fit("gaus", "Q", "R", log10(9e3), log10(1.2e3));
  distInteCQ.GetXaxis()->SetTitle("Log10#left(#sumCQ_i#right) [au]");
  distInteCQ.GetYaxis()->SetTitle("Counts [au]");
  setTH1FontStyle(distInteCQ);
  distInteCQ.Draw();

  double cut = 4.1;

  TLine lineCQ(cut, 0.5, cut, 3e3);
  lineCQ.SetLineColor(kGreen+3);
  lineCQ.SetLineStyle(9);
  lineCQ.Draw();

  TLegend *lgnd = new TLegend(0.16, 0.7, 0.25, 0.95);
  lgnd->AddEntry(&distInteCQ, "Coincidence charge", "h");
  lgnd->AddEntry(&distInteCQ, Form("Entries: %.f", distInteCQ.GetEntries()), "h");
  lgnd->AddEntry(&lineCQ, Form("Cut at: %.f", pow(10,cut)), "l");
  lgnd->AddEntry(distInteCQ.GetFunction("gaus"), Form("#mu = %.f #pm %.f",
        pow(10, distInteCQ.GetFunction("gaus")->GetParameter(1)),
       pow(10, distInteCQ.GetFunction("gaus")->GetParError(1)) ), "l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.055);
  lgnd->Draw();

  cvnsInteCQ.Print("results/distInteCQ.pdf");

  TCanvas cvnsInteCI("cvnsInteCI", "", 1.6e3, 9e2);
  setCanvasStyle(cvnsInteCI);
  cvnsInteCI.cd();
  cvnsInteCI.SetLogy();

  distInteCI.Fit("gaus", "Q", "R", log10(9e3), log10(4e4));
  distInteCI.SetStats(kFALSE);
  distInteCI.GetXaxis()->SetTitle("Log10#left(#sumCI#right) [au]");
  distInteCI.GetYaxis()->SetTitle("Counts [au]");
  setTH1FontStyle(distInteCI);  
  distInteCI.Draw();

  cut = log10(9e3);
  TLine lineCI_1(cut, 0.5, cut, 3e3);
  lineCI_1.SetLineColor(kGreen+3);
  lineCI_1.SetLineStyle(9);
  lineCI_1.Draw();

  double cut2 = log10(4e4);
  TLine lineCI_2(cut2, 0.5, cut2, 3e3);
  lineCI_2.SetLineColor(kGreen+3);
  lineCI_2.SetLineStyle(9);
  lineCI_2.Draw();

  lgnd = new TLegend(0.16, 0.7, 0.3, 0.95);
  lgnd->AddEntry(&distInteCI, "Coincidence height pulse", "h");
  lgnd->AddEntry(&distInteCI, Form("Entries: %.f", distInteCI.GetEntries()), "h");
  lgnd->AddEntry(&lineCI_1, Form("Cut at: %.f", pow(10,cut)), "l");
  lgnd->AddEntry(&lineCI_2, Form("Cut at: %.f", pow(10,cut2)), "l");
  lgnd->AddEntry(distInteCI.GetFunction("gaus"), Form("#mu = %.f #pm %.f",
        pow(10, distInteCI.GetFunction("gaus")->GetParameter(1)),
        pow(10, distInteCI.GetFunction("gaus")->GetParError(1)) ), "l");
  lgnd->SetBorderSize(0);
  lgnd->SetTextSize(0.055);
  lgnd->Draw();

  cvnsInteCI.Print("results/distInteCI.pdf");

  double ave = 0.;
  for ( int j=0; j<3; j++ )
    for ( int i=0; i<150; i++ ) {
      countsCI_st666[j][i] /= tot666Histos[j];
      if ( j==0 )
        ave += countsCI_st666[j][i];
    }
  cout << "MSD " << ave << endl;

  TGraph histo666pmt1(st666_binsCI.size(), &st666_binsCI.front(), 
      &countsCI_st666[0].front());
  TGraph histo666pmt2(st666_binsCI.size(), &st666_binsCI.front(), 
      &countsCI_st666[1].front());
  TGraph histo666pmt3(st666_binsCI.size(), &st666_binsCI.front(), 
      &countsCI_st666[2].front());

  plottingSt666(histo666pmt1, "CI", "1");
  plottingSt666(histo666pmt2, "CI", "2");
  plottingSt666(histo666pmt3, "CI", "3");

  ave = 0;
  for ( int j=0; j<3; j++ )
    for ( int i=0; i<600; i++ ) {
      countsCQ_st666[j][i] /= tot666Histos[j];
      cout << "MSD " << j << " " << i << " " << countsCQ_st666[j][i] << endl;
      if ( j==0 )
        ave += countsCQ_st666[j][i];
    }
  cout << "MSD " << ave << endl;

  TGraph histoCQ666pmt1(st666_binsCQ.size(), &st666_binsCQ.front(), 
      &countsCQ_st666[0].front());
  TGraph histoCQ666pmt2(st666_binsCQ.size(), &st666_binsCQ.front(), 
      &countsCQ_st666[1].front());
  TGraph histoCQ666pmt3(st666_binsCQ.size(), &st666_binsCQ.front(), 
      &countsCQ_st666[2].front());

  plottingSt666(histoCQ666pmt1, "CQ", "1");
  plottingSt666(histoCQ666pmt2, "CQ", "2");
  plottingSt666(histoCQ666pmt3, "CQ", "3");


  return 0;
}

void setCanvasStyle(TCanvas& canvas) {
  canvas.SetTopMargin(0.03);
  canvas.SetBottomMargin(0.17);
  canvas.SetLeftMargin(0.11);
  canvas.SetRightMargin(0.07);
}

void setTH1FontStyle(TH1D& histo) {
  histo.GetXaxis()->SetTitleSize(0.07);
  histo.GetXaxis()->SetTitleOffset(1.2);
  histo.GetXaxis()->SetLabelSize(0.06);
  histo.GetYaxis()->SetTitleSize(0.08);
  histo.GetYaxis()->SetLabelSize(0.06);
  histo.GetYaxis()->SetTitleOffset(0.67);
}

TGraphErrors doHisto(vector<int> counts, vector<double> bins, TString stId,
    TString typeOutlier) {
  vector < double > val;
  vector < double > valEr;
  val.resize( bins.size() );
  valEr.resize( bins.size() );

  int binMax = (typeOutlier=="CI") ? 104 : 404;
  int binMax1 = (typeOutlier=="CI") ? 47 : 197;

  for ( int j = 0; j < binMax; j++ ) {
    val[j] = double(counts[j]);
    valEr[j] = sqrt( double(counts[j]) );
  }
  for (int j = 0; j < binMax1; j++) {
    val[j + binMax-1] = double(counts[j + binMax-1]) / 4.0;
    valEr[j + binMax-1] = sqrt( double(counts[j + binMax-1]) / 4.0 );
  }
  
  TGraphErrors graphHisto(bins.size(), &bins.front(), &val.front(), 
      0, &valEr.front());
  graphHisto.SetTitle("; "+typeOutlier+" [FADC]; Counts [au]"); 

  return graphHisto;
}

void setTGraphFontStyle(TGraph& grph) {
  grph.GetXaxis()->SetTitleSize(0.07);
  grph.GetXaxis()->SetTitleOffset(1.2);
  grph.GetXaxis()->SetLabelSize(0.06);
  grph.GetYaxis()->SetTitleSize(0.08);
  grph.GetYaxis()->SetLabelSize(0.06);
  grph.GetYaxis()->SetTitleOffset(0.67);
}


void plottingSt666(TGraph grph, TString typeHist, TString pmtId) {
  TCanvas cvnspmt("cvnspmt"+pmtId, "", 1.6e3, 9e2);
  cvnspmt.cd();
  setCanvasStyle(cvnspmt);
  cvnspmt.SetLogy();

  grph.SetTitle("; "+typeHist+" [FADC]; Average of counts [au]");
  setTGraphFontStyle(grph);
  grph.Draw("AL");

  TLegend lgnd(0.6, 0.7, 0.7, 0.88);
  lgnd.AddEntry(&grph, "St. 666, PMT "+pmtId, "h");
  if ( pmtId == "1" )
    lgnd.AddEntry(&grph, Form("Average of entries: %d", 5124), "h");
  
  lgnd.SetBorderSize(0);
  lgnd.SetTextSize(0.055);
  lgnd.Draw();

  cvnspmt.Print("results/cumulativeHisto"+typeHist+"St666Pmt"+pmtId+".pdf");
}
