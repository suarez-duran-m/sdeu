TCanvas *canvasStyle(TString name)
{
  TCanvas *canvas = new TCanvas(name, name, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.11); 
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.02); 
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

void histoStyle(TH1F *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void histoStyle(TGraphErrors *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}

void fillingChisPk( TString bname, TString st, int pmt, TH1F *hist, TH1F *hprob)
{
  TString monthUub[] = {"aug", "sep", "oct", "nov"};
  TString pmtId;
  pmtId.Form("%d", pmt);

  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *peakInfo;
  
  double tmpchi2 = 0.;
  double tmpndf = 0.;
  double tmpprob = 0.;
  double tmpchindf = 0.;
  int tmpcnt = 0;

  for ( int month=0; month<4; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    peakInfo = (TTree*)f->Get("PeakData");

    tmpchi2 = 0.;
    tmpndf = 0.;
    tmpprob = 0.;

    peakInfo->SetBranchAddress("chi2", &tmpchi2);
    peakInfo->SetBranchAddress("ndf", &tmpndf);
    peakInfo->SetBranchAddress("prob", &tmpprob);

    for( int etry=0; etry<peakInfo->GetEntries(); etry++)
    {
      peakInfo->GetEntry(etry);
      tmpchindf = tmpchi2 / tmpndf;

      if ( tmpchindf < 12 )
      {
        hist->Fill( tmpchindf );
        hprob->Fill( tmpprob );
      }
      else
      {
        tmpcnt++;
        hist->Fill(12);
      }
    }
  }
  f->Delete();
  cerr << "bad " << pmt << " " << tmpcnt << endl;
}


void fillingChisCh( TString bname, TString st, int pmt, TH1F *hist, TH1F *hprob)
{
  TString monthUub[] = {"aug", "sep", "oct", "nov"};
  TString pmtId;
  pmtId.Form("%d", pmt);

  TString fname = bname + pmtId+st+"Mth";

  TFile *f;
  TTree *chargeInfo;
  
  double tmpchi2 = 0.;
  double tmpndf = 0.;
  double tmpprob = 0.;
  double tmpchindf = 0.;

  for ( int month=0; month<3; month++ )
  {
    f = TFile::Open(fname+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpchi2 = 0.;
    tmpndf = 0.;
    tmpprob = 0.;

    chargeInfo->SetBranchAddress("chi2", &tmpchi2);
    chargeInfo->SetBranchAddress("ndf", &tmpndf);
    chargeInfo->SetBranchAddress("prob", &tmpprob);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      tmpchindf = tmpchi2 / tmpndf;

      if ( tmpchindf < 6 )
      {
        hist->Fill( tmpchindf );
        hprob->Fill( tmpprob );
      }
      else
        hist->Fill(6);
    }
  }
  f->Delete();
}

// =================================
// *** *** *** MAIN CODE *** *** ***

void plottingFitsQltyPkCh(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "ubAoPPMT";

  TPaveStats *ptstats;

  // ====================================
  // *** *** *** Reading Chis *** *** ***

  // ==============================
  // *** *** Doing for Peak *** ***
  TH1F *hChiPkpmt1 = new TH1F("hChiPkpmt1", "", 200, 0, 20);
  TH1F *hChiPkpmt2 = new TH1F("hChiPkpmt2", "", 200, 0, 20);
  TH1F *hChiPkpmt3 = new TH1F("hChiPkpmt3", "", 200, 0, 20);

  TH1F *hProbPkpmt1 = new TH1F("hProbPkpmt1", "", 1e2, 0, 1);
  TH1F *hProbPkpmt2 = new TH1F("hProbPkpmt2", "", 1e2, 0, 1);
  TH1F *hProbPkpmt3 = new TH1F("hProbPkpmt3", "", 1e2, 0, 1);

  fillingChisPk( basename, statId, 1, hChiPkpmt1, hProbPkpmt1 );
  fillingChisPk( basename, statId, 2, hChiPkpmt2, hProbPkpmt2 );
  fillingChisPk( basename, statId, 3, hChiPkpmt3, hProbPkpmt3 );

  // =========================
  // *** Plotting for Chis ***

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  statId = "";
  statId.Form("%d", st);

  hChiPkpmt1->GetXaxis()->SetTitle("#chi^{2} / ndf [au]");
  hChiPkpmt1->GetYaxis()->SetTitle("Counts [au]");
  hChiPkpmt1->GetYaxis()->SetRangeUser(0, 380);// 270);
  hChiPkpmt1->GetXaxis()->SetRangeUser(0, 5); //10);
  hChiPkpmt1->SetLineColor(kBlue);
  hChiPkpmt1->SetLineWidth(2);
  histoStyle(hChiPkpmt1);
  hChiPkpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hChiPkpmt1->SetName("Station "+statId+" PMT1");
  hChiPkpmt1->GetListOfFunctions()->Add(ptstats);

  hChiPkpmt2->SetLineColor(kGreen+2);
  hChiPkpmt2->SetLineWidth(2);
  hChiPkpmt2->Draw("sames");
  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hChiPkpmt2->SetName("Station "+statId+" PMT2");
  hChiPkpmt2->GetListOfFunctions()->Add(ptstats);

  hChiPkpmt3->SetLineColor(kRed+1);
  hChiPkpmt3->SetLineWidth(2);
  hChiPkpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hChiPkpmt3->SetName("Station "+statId+" PMT3");
  hChiPkpmt3->GetListOfFunctions()->Add(ptstats);

  c1->Print("../../plots/ubChisDistPmtsPk"+statId+".pdf");


  // =========================
  // *** Plotting for Prob ***

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  statId = "";
  statId.Form("%d", st);

  c2->SetLogy();
  hProbPkpmt1->GetXaxis()->SetTitle("Prob. [au]");
  hProbPkpmt1->GetYaxis()->SetTitle("Counts [au]");
  hProbPkpmt1->SetLineColor(kBlue);
  hProbPkpmt1->SetFillColor(kBlue);
  hProbPkpmt1->SetFillStyle(3001);
  histoStyle(hProbPkpmt1);
  hProbPkpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hProbPkpmt1->SetName("Station "+statId+" PMT1");
  hProbPkpmt1->GetListOfFunctions()->Add(ptstats);

  hProbPkpmt2->SetLineColor(kGreen+2);
  hProbPkpmt2->SetFillColor(kGreen+2);
  hProbPkpmt2->SetFillStyle(3001);
  hProbPkpmt2->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hProbPkpmt2->SetName("Station "+statId+" PMT2");
  hProbPkpmt2->GetListOfFunctions()->Add(ptstats);

  hProbPkpmt3->SetLineColor(kRed+1);
  hProbPkpmt3->SetFillColor(kRed+1);
  hProbPkpmt3->SetFillStyle(3001);
  hProbPkpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hProbPkpmt3->SetName("Station "+statId+" PMT3");
  hProbPkpmt3->GetListOfFunctions()->Add(ptstats);

  c2->Print("../../plots/ubProbDistPmtsPk"+statId+".pdf");


  // ================================
  // *** *** Doing for Charge *** ***

  TH1F *hChiChpmt1 = new TH1F("hChiChpmt1", "", 100, 0, 10);
  TH1F *hChiChpmt2 = new TH1F("hChiChpmt2", "", 100, 0, 10);
  TH1F *hChiChpmt3 = new TH1F("hChiChpmt3", "", 100, 0, 10);

  TH1F *hProbChpmt1 = new TH1F("hProbChpmt1", "", 1e2, 0, 1);
  TH1F *hProbChpmt2 = new TH1F("hProbChpmt2", "", 1e2, 0, 1);
  TH1F *hProbChpmt3 = new TH1F("hProbChpmt3", "", 1e2, 0, 1);


  statId = "";
  statId.Form("St%d", st);

  fillingChisCh( basename, statId, 1, hChiChpmt1, hProbChpmt1 );
  fillingChisCh( basename, statId, 2, hChiChpmt2, hProbChpmt2 );
  fillingChisCh( basename, statId, 3, hChiChpmt3, hProbChpmt3 );

  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  statId = "";
  statId.Form("%d", st);

  hChiChpmt1->GetXaxis()->SetTitle("#chi^{2} / ndf [au]");
  hChiChpmt1->GetYaxis()->SetTitle("Counts [au]");
  hChiChpmt1->GetYaxis()->SetRangeUser(0, 550); //540);
  hChiChpmt1->GetXaxis()->SetRangeUser(0, 2.5);
  hChiChpmt1->SetLineColor(kBlue);
  hChiChpmt1->SetLineWidth(2);
  histoStyle(hChiChpmt1);
  hChiChpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hChiChpmt1->SetName("Station "+statId+" PMT1");
  hChiChpmt1->GetListOfFunctions()->Add(ptstats);

  hChiChpmt2->SetLineColor(kGreen+2);
  hChiChpmt2->SetLineWidth(2);
  hChiChpmt2->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hChiChpmt2->SetName("Station "+statId+" PMT2");
  hChiChpmt2->GetListOfFunctions()->Add(ptstats);

  hChiChpmt3->SetLineColor(kRed+1);
  hChiChpmt3->SetLineWidth(2);
  hChiChpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hChiChpmt3->SetName("Station "+statId+" PMT3");
  hChiChpmt3->GetListOfFunctions()->Add(ptstats);

  c3->Print("../../plots/ubChisDistPmtsCh"+statId+".pdf");

  // =========================
  // *** Plotting for Prob ***

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  statId = "";
  statId.Form("%d", st);

  hProbChpmt1->GetXaxis()->SetTitle("Prob. [au]");
  hProbChpmt1->GetYaxis()->SetTitle("Counts [au]");
  hProbChpmt1->GetYaxis()->SetRangeUser(0, 110);
  hProbChpmt1->SetLineColor(kBlue);
  hProbChpmt1->SetFillColor(kBlue);
  hProbChpmt1->SetFillStyle(3001);
  histoStyle(hProbChpmt1);
  hProbChpmt1->Draw();

  ptstats = new TPaveStats(0.73, 0.77, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hProbChpmt1->SetName("Station "+statId+" PMT1");
  hProbChpmt1->GetListOfFunctions()->Add(ptstats);

  hProbChpmt2->SetLineColor(kGreen+2);
  hProbChpmt2->SetFillColor(kGreen+2);
  hProbChpmt2->SetFillStyle(3001);
  hProbChpmt2->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.55, 0.96, 0.75,"brNDC");
  ptstats->SetTextColor(kGreen+2);
  hProbChpmt2->SetName("Station "+statId+" PMT2");
  hProbChpmt2->GetListOfFunctions()->Add(ptstats);

  hProbChpmt3->SetLineColor(kRed+1);
  hProbChpmt3->SetFillColor(kRed+1);
  hProbChpmt3->SetFillStyle(3001);
  hProbChpmt3->Draw("sames");

  ptstats = new TPaveStats(0.73, 0.33, 0.96, 0.53,"brNDC");
  ptstats->SetTextColor(kRed+1);
  hProbChpmt3->SetName("Station "+statId+" PMT3");
  hProbChpmt3->GetListOfFunctions()->Add(ptstats);

  c4->Print("../../plots/ubProbDistPmtsCh"+statId+".pdf");
}
