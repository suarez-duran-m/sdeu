/*
#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
*/
 
void binsBl() {
   // Create a histogram for the values we read.
   //auto myHist = new TH1F("h1","ntuple",100,-4,4);
   auto pmt1hmeanf = new TH1F("pmt1hmeanf","Mean Baseline for First 100 bins PMT1 High",5356800, 1606780800,1612137600);
   auto pmt1hmeanl = new TH1F("pmt1hmeanl","Mean Baseline for First 100 bins PMT1 High",5356800, 1606780800,1612137600);
   auto pmt1hrmsf = new TH1F("pmt1hrmsf","Baseline RMS for Last 100 bins PMT1 High",5356800, 1606780800,1612137600);
   auto pmt1hrmsl = new TH1F("pmt1hrmsl","Baseline RMS for Las 100 bins PMT1 High",5356800, 1606780800,1612137600);
   auto meanDifftime = new TH1F("meanDifftime","Difference for Mean of First and Last 100 bins Baseline PMT1 High",5356800, 1606780800,1612137600);
   auto rmsDifftime = new TH1F("rmsmeanDifftime","Difference for RMS of First and Last 100 bins Baseline PMT1 High",5356800, 1606780800,1612137600);
   auto event = new TH1I ("event", "Events IDs for Station 1211", 615755, 61190079, 61805835);

    // Open the file containing the tree.
   auto blFile = TFile::Open("bl1211.root");
   if (!blFile || blFile->IsZombie()) {
      return;
   }

   pmt1hmeanf = (TH1F*)blFile->Get("pmt1lbmeanf");
   pmt1hmeanl = (TH1F*)blFile->Get("pmt1lbmeanl");
   pmt1hrmsf = (TH1F*)blFile->Get("pmt1hbrmsf");
   pmt1hrmsl = (TH1F*)blFile->Get("pmt1hbrmsl");
   event = (TH1I*)blFile->Get("eventId");


   auto diffMean= new TH1F("diffMean", "Diff for Mean PMT1 High", 10000, 0, 10);
   auto diffRms = new TH1F("diffRms", "Diff for RMS PMT1 High", 10000, -0, 100);

   for (int i=0; i<pmt1hmeanf->GetNbinsX(); i++)
     if (pmt1hmeanf->GetBinContent(i) > 0 ){
       diffMean->Fill( abs(pmt1hmeanf->GetBinContent(i)-pmt1hmeanl->GetBinContent(i)) );
       meanDifftime->Fill( pmt1hmeanf->GetXaxis()->GetBinCenter(i), abs(pmt1hmeanf->GetBinContent(i)-pmt1hmeanl->GetBinContent(i)) );
     }

   for (int i=0; i<pmt1hrmsf->GetNbinsX(); i++)
     if ( pmt1hrmsf->GetBinContent(i) > 0. ){
       diffRms->Fill( abs(pmt1hrmsf->GetBinContent(i)-pmt1hrmsl->GetBinContent(i)) );
       rmsDifftime->Fill( pmt1hrmsf->GetXaxis()->GetBinCenter(i), abs(pmt1hrmsf->GetBinContent(i)-pmt1hrmsl->GetBinContent(i)) );
       if ( pmt1hrmsf->GetBinContent(i)-pmt1hrmsl->GetBinContent(i) > 3000. )
         cout << pmt1hrmsf->GetBinContent(i)-pmt1hrmsl->GetBinContent(i) << " "
           << int(event->GetBinContent(i)) << endl;
     }
 
   auto c1 = new TCanvas("c1","Time Series for Baseline Mean and RMS",10,10,1200,700);
   c1->Divide(1,2);

   c1->cd(1);
   c1->GetPad(1)->SetGrid();

   pmt1hmeanf->SetMarkerStyle(20);
   pmt1hmeanf->SetMarkerColor(1);
   pmt1hmeanf->GetYaxis()->SetRangeUser(200, 280);
   pmt1hmeanf->Draw("P");
   
   pmt1hmeanl->SetMarkerStyle(20);
   pmt1hmeanl->SetMarkerColor(2);
   pmt1hmeanl->Draw("SAME P");

   c1->cd(2);
   c1->GetPad(2)->SetGrid();
   pmt1hrmsf->GetXaxis()->SetTimeDisplay(1);
   pmt1hrmsf->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   pmt1hrmsf->GetXaxis()->SetLabelOffset(0.020);
   pmt1hrmsf->GetXaxis()->SetLabelSize(0.03);
   pmt1hrmsf->GetYaxis()->SetTitle("Baseline RMS (First 100 bins) / FADC");
   pmt1hrmsf->SetMarkerColor(1);
   pmt1hrmsf->SetMarkerStyle(20);
   pmt1hrmsf->Draw("P");
   
   pmt1hrmsl->GetXaxis()->SetTimeDisplay(1);
   pmt1hrmsl->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");   
   pmt1hrmsl->SetMarkerColor(2);
   pmt1hrmsl->SetMarkerStyle(20);
   pmt1hrmsl->Draw("SAME P");

   c1->Print("timeBinsMeanRMS1211.pdf");

   auto c2 = new TCanvas("c2","Baseline Mean and RMS Distribution",10,10,1200,700);
   c2->Divide(1, 2);

   c2->cd(1);
   c2->GetPad(1)->SetGrid();
   diffMean->GetYaxis()->SetTitle("Counts / au");
   diffMean->GetXaxis()->SetTitle("Difference for the Mean of First and Last 100 bins of Baseline / FADC");
   diffMean->Draw();
   
   c2->cd(2);
   c2->GetPad(2)->SetGrid();

   int binMax = int( diffRms->GetMean() + 1.*diffRms->GetRMS() );
   int binMin = int( diffRms->GetMean() - 1.*diffRms->GetRMS() );
   diffRms->Fit("gaus","","", binMin, binMax);
   TF1 *gMean = (TF1*) diffRms->GetListOfFunctions()->FindObject("gaus");
   diffRms->GetYaxis()->SetTitle("Counts / au");
   diffRms->GetXaxis()->SetTitle("Difference for the Mean of First and Last 100 bins of Baseline RMS / FADC");
   diffRms->Draw();

   c2->Print("diffBinsMeanRms1211.pdf");

   auto c3 = new TCanvas("c3","Temporal Difference for Mean and RMS",10,10,1200,700);
   c3->Divide(1, 3);
   c3->cd(1);
   c3->GetPad(1)->SetGrid();
   meanDifftime->GetXaxis()->SetTimeDisplay(1);
   meanDifftime->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   meanDifftime->GetXaxis()->SetLabelOffset(0.020);
   meanDifftime->GetXaxis()->SetLabelSize(0.03);
   meanDifftime->SetMarkerStyle(20);
   meanDifftime->Draw("hist P");

   c3->cd(2);
   c3->GetPad(2)->SetGrid();
   rmsDifftime->GetXaxis()->SetTimeDisplay(1);
   rmsDifftime->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   rmsDifftime->GetXaxis()->SetLabelOffset(0.020);
   rmsDifftime->GetXaxis()->SetLabelSize(0.03);
   rmsDifftime->SetMarkerStyle(20);
   rmsDifftime->Draw("hist P");

   c3->cd(3);
   c3->GetPad(3)->SetGrid();
   event->GetXaxis()->SetTimeDisplay(1);
   event->GetXaxis()->SetTimeFormat("#splitline{%m/%d/%y}{%H:%M:%S}%F1970-01-01 05:00:00");
   event->GetXaxis()->SetLabelOffset(0.020);
   event->GetXaxis()->SetLabelSize(0.03);
   event->SetMarkerStyle(20);
   event->Draw("hist P");
}
