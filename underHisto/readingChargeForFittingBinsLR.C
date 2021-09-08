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
  hist->GetYaxis()->SetTitleOffset(0.9);
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


void histoStyle(TGraph *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


vector < double > fillingCh( TString bname, TString st, TString strlrb, int pmt, TString whichInfo)
{
  TString monthUub[] = {"dec", "jan", "feb", "mar", "apr", "may", "jun", "jul"};
  TString pmtId;
  pmtId.Form("%d", pmt);
  TString fname = bname + pmtId+st+"lrb";

  TFile *f;
  TTree *chargeInfo;
  double tmpVals = 0.;
  vector < double > returnVals;
  int tmp = 7;

  for ( int month=0; month<tmp+1; month++ )
  {
    f = TFile::Open(fname+strlrb+monthUub[month]+".root");
    chargeInfo = (TTree*)f->Get("ChargeData");

    tmpVals = 0.;
    chargeInfo->SetBranchAddress(whichInfo, &tmpVals);

    for( int etry=0; etry<chargeInfo->GetEntries(); etry++)
    {
      chargeInfo->GetEntry(etry);
      if ( tmpVals > 0 )
        returnVals.push_back( tmpVals );
      //if ( whichInfo=="chargeValDer" && tmpVals > 1400 )
        //cout << monthUub[month] << endl;
    }
  }
 
  chargeInfo->Delete();
  f->Delete();
  return returnVals;
}

double getmean( vector<double> arr )
{
  double mean = 0.;
  int nb = arr.size();
  int ngoodb = 0;
    for (int i=0; i<nb; i++)
      if ( arr[i] > 100 )
      {
        mean += arr[i];
        ngoodb++;
      }
  return mean/ngoodb;
}

double getrms( vector<double> arr, double meanarr )
{
  double rms = 0.;
  int nb = arr.size();
  int ngoodb = 0;
  for (int i=0; i<nb; i++)
    if ( arr[i] > 100 )
    {
      rms += (arr[i] - meanarr)*(arr[i] - meanarr);
      ngoodb++;
    }
  return sqrt(rms/ngoodb)/meanarr;
}

// =================================
// *** *** *** MAIN CODE *** *** ***

void readingChargeForFittingBinsLR(int st)
{
  TString statId;
  statId.Form("St%d", st);
  TString basename = "uubChPkPMT";
  TString strlrb;

  TPaveStats *ptstats;
  TLine *line;
  TLegend *leg;

  bool gettime = true;
  TString whinfo = "chargeVal";

  vector < double > chargePmt1;
  vector < double > chargePmt2;
  vector < double > chargePmt3;
  vector < double > chargeDerPmt1;
  vector < double > chargeDerPmt2;
  vector < double > chargeDerPmt3;

  int nLrb = 10;
  int nPoints = nLrb;
  double tmp = 0.;
  double xlrb [ nLrb ];
  double yvemFitPmt1 [ nLrb ];
  double yvemFitPmt2 [ nLrb ];
  double yvemFitPmt3 [ nLrb ];
  double yvemDer;

  TH1F *vemDistPmt1Fit40 = new TH1F ("vemDistPmt1Fit40", "", 2000, 0, 2000);
  TH1F *vemDistPmt2Fit40 = new TH1F ("vemDistPmt2Fit40", "", 2000, 0, 2000);
  TH1F *vemDistPmt3Fit40 = new TH1F ("vemDistPmt3Fit40", "", 2000, 0, 2000);

  for ( int lrb=5; lrb<=50; lrb+=5 )
  {
    strlrb.Form("%d", lrb);
    whinfo = "chargeVal";
    chargePmt1 = fillingCh( basename, statId, strlrb, 1, whinfo);
    chargePmt2 = fillingCh( basename, statId, strlrb, 2, whinfo);
    chargePmt3 = fillingCh( basename, statId, strlrb, 3, whinfo);
    xlrb[int(lrb/5)-1] = lrb;
    tmp = getmean( chargePmt1 );
    yvemFitPmt1[int(lrb/5)-1] = getrms( chargePmt1, tmp );
    tmp = getmean( chargePmt2 );
    yvemFitPmt2[int(lrb/5)-1] = getrms( chargePmt2, tmp );
    tmp = getmean( chargePmt3 );
    yvemFitPmt3[int(lrb/5)-1] = getrms( chargePmt3, tmp );
    if ( lrb==40 )
    {
      for ( int i=0; i<chargePmt2.size(); i++ )
        vemDistPmt1Fit40->Fill( chargePmt2[i] );
      //for ( int i=0; i<chargePmt2.size(); i++ )
        //vemDistPmt2Fit40->Fill( chargePmt2[i] );
      for ( int i=0; i<chargePmt3.size(); i++ )
         vemDistPmt3Fit40->Fill( chargePmt3[i] );
    }
    if ( lrb==25 )
      for ( int i=0; i<chargePmt2.size(); i++ )
        vemDistPmt2Fit40->Fill( chargePmt2[i] );
  }
  whinfo = "chargeValDer";
  chargeDerPmt1 = fillingCh( basename, statId, "40", 1, whinfo);
  tmp = getmean( chargeDerPmt1 );
  yvemDer = getrms( chargeDerPmt1, tmp );

  TGraph *chFitPmt1 = new TGraph(nPoints, xlrb, yvemFitPmt1);
  TGraph *chFitPmt2 = new TGraph(nPoints, xlrb, yvemFitPmt2);
  TGraph *chFitPmt3 = new TGraph(nPoints, xlrb, yvemFitPmt3);

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();
  chFitPmt1->SetTitle("");
  chFitPmt1->GetXaxis()->SetTitle("Bins for fit (towards left, towards right) [au]");
  chFitPmt1->GetYaxis()->SetTitle("RMS/Mean for VEM-Charge [FADC]");
  chFitPmt1->GetYaxis()->SetRangeUser(0, .08);
  chFitPmt1->GetXaxis()->SetRangeUser(18, 52);
  chFitPmt1->SetMarkerStyle(8);
  chFitPmt1->SetMarkerSize(2);
  chFitPmt1->SetMarkerColor(kAzure+10);
  histoStyle(chFitPmt1);
  chFitPmt1->Draw("AP same");

  chFitPmt2->SetMarkerStyle(8);
  chFitPmt2->SetMarkerSize(2);
  chFitPmt2->SetMarkerColor(kGreen+3);
  chFitPmt2->Draw("P same");

  chFitPmt3->SetMarkerStyle(8);
  chFitPmt3->SetMarkerSize(2);
  chFitPmt3->SetMarkerColor(kOrange+10);
  chFitPmt3->Draw("P same");

  leg = new TLegend(0.75,0.7,0.95,0.95);
  leg->AddEntry(chFitPmt1, "PMT 1", "p");
  leg->AddEntry(chFitPmt2, "PMT 2", "p");
  leg->AddEntry(chFitPmt3, "PMT 3", "p");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  //c1->Print("../plots/uubChRmsFitBinsLrSt863zoom.pdf");
  c1->Print("../plots/uubChRmsFitBinsLrSt863.pdf");


  TCanvas *c2 = canvasStyle("c2"); 
  c2->cd();

  vemDistPmt1Fit40->SetTitle("");
  vemDistPmt1Fit40->GetYaxis()->SetRangeUser(0, 62);
  vemDistPmt1Fit40->GetYaxis()->SetTitle("Counts [au]");
  vemDistPmt1Fit40->GetXaxis()->SetRangeUser(1.1e3, 1.7e3);
  vemDistPmt1Fit40->GetXaxis()->SetTitle("VEM Charge [FADC]");
  histoStyle(vemDistPmt1Fit40);
  vemDistPmt1Fit40->SetLineWidth(2);
  vemDistPmt1Fit40->SetLineColor(kAzure+10);
  vemDistPmt1Fit40->Draw();

  vemDistPmt2Fit40->SetLineColor(kGreen+4);
  vemDistPmt2Fit40->SetLineWidth(2);
  vemDistPmt2Fit40->Draw("same");

  vemDistPmt3Fit40->SetLineColor(kOrange+10);
  vemDistPmt3Fit40->SetLineWidth(2);
  vemDistPmt3Fit40->Draw("same");

  //cout << vemDistPmt2Fit40->GetMean()
}
