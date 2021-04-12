/*
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
*/
 
void readBl() {
   // Create a histogram for the values we read.
   //auto myHist = new TH1F("h1","ntuple",100,-4,4);
   auto pmt0hbmean = new TH1F("pmt0hbmean","Mean Baseline for PMT0 High",5356800, 1606780800,1612137600);
   auto pmt0hbrms = new TH1F("pmt0hbrms","Baseline RMS for PMT0 High",5356800, 1606780800,1612137600);

    // Open the file containing the tree.
   auto blFile = TFile::Open("bl1211.root");
   if (!blFile || blFile->IsZombie()) {
      return;
   }

   auto filtMean = new TH1F("filtMean","Baseline RMS for PMT0 High",5356800, 1606780800,1612137600);
   auto filtRms = new TH1F("filtRms","Baseline RMS for PMT0 High",5356800, 1606780800,1612137600);


    // Open the file containing the tree.
   auto blFilter = TFile::Open("msdbl1211.root");
   if (!blFilter || blFile->IsZombie()) {
      return;
   }

   filtMean = (TH1F*)blFilter->Get("pmt0bh");
   filtRms = (TH1F*)blFilter->Get("pmt0bhrms");

   pmt0hbmean = (TH1F*)blFile->Get("pmt0bh");
   pmt0hbrms = (TH1F*)blFile->Get("pmt0bhrms");

   int nbinsMean = int(pmt0hbmean->GetMaximum() + 1);
   int nbinsRms = int(pmt0hbrms->GetMaximum() + 1);

   auto meanDist = new TH1F("meanDist", "Mean Distribution for PMT0 High", int(nbinsMean*100), 1, nbinsMean+1);   
   auto rmsDist = new TH1F("rmsDist", "RMS Distribution for PMT0 High", int(nbinsRms*100), 1, nbinsRms+1);

   for (int i=0; i<pmt0hbrms->GetNbinsX(); i++)
     rmsDist->Fill( pmt0hbrms->GetBinContent(i) );

   for (int i=0; i<pmt0hbmean->GetNbinsX(); i++)
     meanDist->Fill( pmt0hbmean->GetBinContent(i) ); 


   int binMax = int( meanDist->GetMean() + 4*meanDist->GetRMS() );
   int binMin = int( meanDist->GetMean() - 4*meanDist->GetRMS() );
   meanDist->Fit("gaus","","", binMin, binMax);
   TF1 *gMean = (TF1*) meanDist->GetListOfFunctions()->FindObject("gaus");
   float fitBlMm = gMean->GetParameter(1);
   float fitBlMs = gMean->GetParameter(2);

   binMax = int( rmsDist->GetMean() + 5*rmsDist->GetRMS() );
   binMin = int( rmsDist->GetMean() - 5*rmsDist->GetRMS() );
   rmsDist->Fit("gaus","","", binMin, binMax);
   TF1 *gRms = (TF1*) rmsDist->GetListOfFunctions()->FindObject("gaus");
   float fitBlRm = gRms->GetParameter(1);
   float fitBlRs = gRms->GetParameter(2);

   binMax = int( fitBlRm + 3.*fitBlRs );
   binMin = int( fitBlRm - 3.*fitBlRs );

   cout << fitBlRm - 3.*fitBlRs << " " 
     << fitBlRm + 3.*fitBlRs 
     << endl;
   
   auto c1 = new TCanvas("c1","Time Series for Baseline Mean and RMS",10,10,1200,700);
   c1->Divide(1,2);

   c1->cd(1);
   c1->GetPad(1)->SetGrid();

   pmt0hbmean->SetMarkerStyle(20);
   pmt0hbmean->SetMarkerColor(1);
   pmt0hbmean->GetYaxis()->SetRangeUser(200, 280);
   pmt0hbmean->Draw("P");
   
   filtMean->SetMarkerStyle(20);
   filtMean->SetMarkerColor(2);
   filtMean->Draw("SAME P");

   c1->cd(2);
   c1->GetPad(2)->SetGrid();
   pmt0hbrms->GetXaxis()->SetTimeDisplay(1);
   pmt0hbrms->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   pmt0hbrms->GetXaxis()->SetLabelOffset(0.020);
   pmt0hbrms->GetXaxis()->SetLabelSize(0.03);
   pmt0hbrms->GetYaxis()->SetTitle("Baseline RMS (First 100 bins) / FADC");
   pmt0hbrms->SetMarkerColor(1);
   pmt0hbrms->SetMarkerStyle(20);
   pmt0hbrms->Draw("P");

   //filtRms->GetXaxis()->SetTimeDisplay(1);
   //filtRms->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   
   filtRms->GetXaxis()->SetTimeDisplay(1);
   filtRms->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");   
   filtRms->SetMarkerColor(2);
   filtRms->SetMarkerStyle(20);
   filtRms->Draw("SAME P");

   c1->Print("timeMeanRMS1211.pdf");

   auto c2 = new TCanvas("c2","Baseline Mean and RMS Distribution",10,10,1200,700);
   c2->Divide(1, 2);

   c2->cd(1);
   c2->GetPad(1)->SetGrid();
   meanDist->GetYaxis()->SetTitle("Counts / au");
   meanDist->GetXaxis()->SetTitle("Baseline Mean (First 100 bins) / FADC");
   meanDist->Draw();
   
   c2->cd(2);
   c2->GetPad(2)->SetGrid();
   rmsDist->GetYaxis()->SetTitle("Counts / au");
   rmsDist->GetXaxis()->SetTitle("Baseline RMS / FADC");
   rmsDist->Draw();

   c2->Print("distMeanRms1211.pdf");

  
   auto c0 = new TCanvas("c0", "Filter Mean and RMS", 10, 10, 1200, 700);
   c0->Divide(1, 2);
   
   c0->cd(1);
   c0->GetPad(1)->SetGrid();
   filtMean->GetYaxis()->SetRangeUser(200, 280);
   filtMean->Draw("P");

   c0->cd(2);
   c0->GetPad(2)->SetGrid();
   filtRms->GetXaxis()->SetTimeDisplay(1);
   filtRms->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   filtRms->GetXaxis()->SetLabelOffset(0.020);
   filtRms->GetXaxis()->SetLabelSize(0.03);
   filtRms->GetYaxis()->SetTitle("Baseline RMS (First 100 bins) / FADC");
   filtRms->GetYaxis()->SetRangeUser(1, 3);
   filtRms->Draw("P");
   
   c0->Print("filterMeanRms1211.pdf");

   /*
   // Create a TTreeReader for the tree, for instance by passing the
   // TTree's name and the TDirectory / TFile it is in.
   TTreeReader myReader("ntuple", myFile);
 
   // The branch "px" contains floats; access them as myPx.
   //TTreeReaderValue<Float_t> myPx(myReader, "px");
   TTreeReaderValue<Float_t> myPx(myReader, "pmt0");
   // The branch "py" contains floats, too; access those as myPy.
   //TTreeReaderValue<Float_t> myPy(myReader, "py");
 
   // Loop over all entries of the TTree or TChain.
   while (myReader.Next()) {
      // Just access the data as if myPx and myPy were iterators (note the '*'
      // in front of them):
      pmt0->Fill(*myPx);
   }
 
   auto f = new TFile("ht.root");
   auto T = (TTree*)f->Get("T");
   auto c1 = new TCanvas("c1","test",10,10,600,400);
   T->Draw("hpx.GetRMS():hprof.GetMean()");
   c1->Print("htr3.png");

   pmt0->Draw();
   */
}
