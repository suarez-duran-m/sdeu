TCanvas *canvasStyle(TString name) {
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

void histoStyle(TH1F *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void histoStyle(TGraphErrors *hist) {
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


TH1F *getSmooth(TH1F &hist, double xb[]) {
	unsigned int nb = 600;
  double yi = 0.;

  TString tmpname;
  tmpname.Form("%d",rand());

	TH1F *hstSmooth = new TH1F(tmpname, "", nb, xb);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b>6 && b<nb-7 )
    {
      yi = hist.GetBinContent(b+1 - 7) 
        + hist.GetBinContent(b+1 - 6) 
        + hist.GetBinContent(b+1 - 5) 
        + hist.GetBinContent(b+1 - 4) 
        + hist.GetBinContent(b+1 - 3) 
        + hist.GetBinContent(b+1 - 2)
        + hist.GetBinContent(b+1 - 1) 
        + hist.GetBinContent(b+1) 
        + hist.GetBinContent(b+1 + 1) 
        + hist.GetBinContent(b+1 + 2)
        + hist.GetBinContent(b+1 + 3)
        + hist.GetBinContent(b+1 + 4)
        + hist.GetBinContent(b+1 + 5)
        + hist.GetBinContent(b+1 + 6)
        + hist.GetBinContent(b+1 + 7);
      yi = yi/15.;
    }
    else
      yi = hist.GetBinContent(b+1);
    
    hstSmooth->SetBinContent(b+1, yi);
  }
  return hstSmooth;
}


TH1F *histDerivative(TH1F &hist, double xb[]) {// Central differences
  int nbins = 600;
  int h = 0;
  TString tmpname;
  tmpname.Form("%d",rand());
  TH1F *derihist = new TH1F(tmpname, "", nbins, xb);
  double der = 0.;
  for ( int kk=1; kk<nbins-1; kk++ )
  {
    h = ( hist.GetBinCenter( kk+1 ) - hist.GetBinCenter( kk ) );
    der = ( hist.GetBinContent( kk+1 ) -  hist.GetBinContent( kk-1 ) )/ (2.*h);
    derihist->SetBinContent( kk, der );

  }
  return derihist;
}


vector < double > getFitRange( TH1F &h ) {
  vector < double > minmax;
  int minRng = 0; // Min for fitting
	int maxRng = 0; // Max for fitting

  double binMin = 0.;
  double binMax = 0.;
  double rawbinMax = 0.;
  double rawbinMin = 0.;
  int rightleftBins = 35;

  for ( int kk=300; kk>80; kk-- ) // From 300 FADC backward
    if ( h.GetBinContent(kk) > 0 ) {
      if ( h.GetBinContent(kk) < fabs(h.GetBinContent(kk+1)) ) {
        binMax = h.GetBinCenter(kk); // tmp FADC for VEM
        rawbinMax = kk; // Bin for tmp VEM
      }
      else {
        binMax = h.GetBinCenter(kk+1);
        rawbinMax = kk; // Bin for tmp VEM
      }
      break;
    }
  
  minRng = h.GetBinCenter(rawbinMax-rightleftBins);
  maxRng = h.GetBinCenter(rawbinMax+rightleftBins);

  minmax.push_back( minRng );
  minmax.push_back( maxRng );
  minmax.push_back( rawbinMax-rightleftBins );
  minmax.push_back( rawbinMax );
  minmax.push_back( binMax ); 

  cout << "MSD: " << binMax << " " << rawbinMax << " " << minRng << " " << maxRng << endl;

  return minmax;
}

void getResiduals( TH1F *hist, TF1 *func,
    double rangMin, double rangMax,
    vector < double > &x, vector < double > &y, vector < double > &err ) {

  int nbins = hist->GetXaxis()->GetNbins();
  double tmp = 0.;
  for ( int kbin=1; kbin<nbins; kbin++ )
    if ( hist->GetBinCenter(kbin) >= rangMin && hist->GetBinCenter(kbin) <= rangMax )
    {
      x.push_back( hist->GetBinCenter(kbin) );
      tmp = func->Eval( hist->GetBinCenter(kbin) ) - hist->GetBinContent(kbin);
      y.push_back( tmp / sqrt( hist->GetBinContent(kbin) ) );
      err.push_back( 
          sqrt( pow(sqrt( hist->GetBinContent(kbin) ),2)
            + pow(sqrt( sqrt(func->Eval( hist->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( hist->GetBinContent(kbin) )
          );
    }
}


// ==================================
// *** ***  *** MAIN CODE *** *** *** 

void readChargeForFitting() {
  TString cdasmeName = "ubChPkPMT1St863lrb35aug.root";
  TFile *cdasmeF = TFile::Open(cdasmeName);
  TTree *cdasmeInfo = (TTree*)cdasmeF->Get("ChargeData");

  TH1F *charge = new TH1F();
  TGraphErrors *graphFitted = new TGraphErrors();
  int evt;
  double chChi2;
  int chNdf;

  cdasmeInfo->SetBranchAddress("chargeForFit", &charge);
  cdasmeInfo->SetBranchAddress("graph", &graphFitted);
  cdasmeInfo->SetBranchAddress("eventId", &evt);
  cdasmeInfo->SetBranchAddress("chi2", &chChi2);
  cdasmeInfo->SetBranchAddress("ndf", &chNdf);

  cdasmeInfo->GetEntry(5);
  cout << evt << endl;

  int nXbins = charge->GetXaxis()->GetNbins();
  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins; b++ )
    xfadc[b] = charge->GetBinCenter(b);

  TH1F *chargeSmooth = getSmooth(*charge, xfadc);
  TH1F *chargeSmooDer = histDerivative(*chargeSmooth, xfadc);
  TH1F *chargeSmooDerSmth = getSmooth(*chargeSmooDer, xfadc);
  chargeSmooDerSmth = getSmooth(*chargeSmooDerSmth, xfadc);

  TLine *line;
  TLegend *leg;
  TPaveStats *ptstats;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  charge->SetStats(0);
  charge->SetTitle("");
  charge->SetLineColor(kBlue);
  charge->GetXaxis()->SetRangeUser(0, 6e2);
  charge->GetYaxis()->SetRangeUser(0, 1550);
  charge->SetLineWidth(1);
  charge->GetXaxis()->SetTitle("[FADC]");
  charge->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(charge);
  charge->Draw();

  chargeSmooth->SetLineColor(kBlack); //Orange+9);
  chargeSmooth->SetLineWidth(1);
  chargeSmooth->Draw("same");

  TString strEvt;
  strEvt.Form("%d", evt); 

  leg = new TLegend(0.62,0.65,0.95,0.96);
  leg->SetHeader("#splitline{Charge histogrma Station 863}{(Event "+strEvt+")}");
  leg->AddEntry(charge,"Charge histogram","l");
  leg->AddEntry(chargeSmooth,"H_{S} histogram","l");
  leg->Draw();  
  
  //c1->Print("../../plots/ubChargeHisto863.pdf");
  

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  chargeSmooDer->SetStats(0);
  chargeSmooDer->SetLineColor(kBlack);
  chargeSmooDer->SetLineWidth(1);
  chargeSmooDer->GetXaxis()->SetRangeUser(0, 6e2);
  chargeSmooDer->GetXaxis()->SetTitle("[FADC]");
  chargeSmooDer->GetYaxis()->SetRangeUser(-25, 12);
  chargeSmooDer->GetYaxis()->SetTitle("[au]");
  histoStyle(chargeSmooDer);
  chargeSmooDer->Draw();

  chargeSmooDerSmth->SetStats(0);
  chargeSmooDerSmth->SetLineColor(kOrange+9);
  chargeSmooDerSmth->SetLineWidth(1);
  chargeSmooDerSmth->Draw("sames");

  line = new TLine(0, 0, 3e3, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(120, -25, 120, 12);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->SetLineColor(kGreen+3);
  line->Draw();

  leg = new TLegend(0.6,0.3,0.95,0.6);
  leg->AddEntry(chargeSmooDer,"H_{DS}","l");
  leg->AddEntry(chargeSmooDerSmth,"H_{SDS}","l");
  leg->AddEntry(chargeSmooDer, "First estimated", "");
  leg->AddEntry(chargeSmooDer, "VEM-hump (green line)", "");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();
  //c2->Print("../plots/chargeDerHisto863.pdf");

  // ===============
  // *** Fitting ***

	double rangXmin = 0; // Min for fitting
	double rangXmax = 0; // Max for fitting
  vector < double > rangeValues;  

  rangeValues = getFitRange(*chargeSmooDerSmth);
  rangXmin = rangeValues[0];
  rangXmax = rangeValues[1];

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( unsigned int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( charge->GetBinContent(b) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( charge->GetBinCenter(b) );
  }
	
	TGraphErrors* chToFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2;

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;
  double tmp = 0.;
  int bigRsd = 0;
  double reduceFactor = 0.05;
  double chi2Ndf = 500.; 
  double bestXmin = 0.;
  double bestXmax = 0.;

  poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  chToFit->Fit("poly2","QR");
                                                                                          
  xResid.clear();
  yResid.clear();
  errResid.clear();
                                                                                          
  // Calculating residuals
  for ( int kbn=1; kbn<nXbins; kbn++ )
    if ( charge->GetBinCenter(kbn) >= rangXmin && charge->GetBinCenter(kbn) <= rangXmax )
    {
      xResid.push_back( charge->GetBinCenter(kbn) );
      tmp = poly2->Eval( charge->GetBinCenter(kbn) ) - charge->GetBinContent(kbn);
      yResid.push_back( tmp / sqrt( charge->GetBinContent(kbn) ) );
      errResid.push_back( sqrt(
            pow(sqrt( charge->GetBinContent(kbn) ),2) 
            + pow(sqrt( sqrt(poly2->Eval( charge->GetBinCenter(kbn) ) ) ),2)
            ) / sqrt( charge->GetBinContent(kbn) ) );
    }

  TGraphErrors *residPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  double vemPk = -1.*poly2->GetParameter(1)/(2.*poly2->GetParameter(0));
 
  TCanvas *c3 = canvasStyle("c3");
  c3->cd();

  chToFit->SetTitle("");
  chToFit->GetXaxis()->SetRangeUser(0, 6e2);
  chToFit->GetYaxis()->SetRangeUser(0, 1550);
  chToFit->SetLineColor(kGray+2);
  chToFit->SetLineWidth(1);
  chToFit->GetXaxis()->SetTitle("[FADC]");
  chToFit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(chToFit);
  chToFit->Draw();

  poly2->SetLineWidth(4);
  poly2->SetLineColor(kBlue);
  poly2->Draw("same");

  TString strVem;
  strVem.Form("%.2f", vemPk);
  leg = new TLegend(0.42,0.7,0.88,0.92);
  TLegendEntry* l1 = leg->AddEntry(poly2,"VEM Fitting-Poly2 (35-bin): "+strVem,"");
  l1->SetTextColor(kBlue);
  strVem.Form("%.2f", rangeValues[4]);
  l1 = leg->AddEntry(poly2,"VEM-hump.: "+strVem,"");
  l1->SetTextColor(kGreen+3);
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->Draw();

  line = new TLine(vemPk, 0, vemPk, 1550);
  line->SetLineColor(kBlue);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(rangeValues[4], 0, rangeValues[4], 1550);
  line->SetLineColor(kGreen+3);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();

  //c3->Print("../plots/chargeFitPoly2863.pdf");

  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  residPoly2->SetTitle("");
  residPoly2->GetXaxis()->SetRangeUser(rangXmin-5, rangXmax+5);
  residPoly2->GetYaxis()->SetRangeUser(-6, 4);
  residPoly2->SetLineColor(kBlack);
  residPoly2->SetLineWidth(1);
  residPoly2->GetXaxis()->SetTitle("[FADC]");
  residPoly2->GetYaxis()->SetTitle("Residuals");
  residPoly2->SetMarkerSize(3);
  residPoly2->SetMarkerColor(kBlack);
  residPoly2->SetMarkerStyle(23);
  histoStyle(residPoly2);
  residPoly2->GetYaxis()->SetTitleOffset(0.8);
  residPoly2->Draw("APL");

  TString chindf;
  TString strVemPk;
  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
  strVemPk.Form("%.2f", vemPk);
  TString strChi;
  TString strNdf; 
  strChi.Form("%.2f", poly2->GetChisquare());
    strNdf.Form("%d", poly2->GetNDF());

  leg = new TLegend(0.5,0.19,0.9,0.40);
  l1 = leg->AddEntry(residPoly2,"#frac{#chi^{2}}{ndf} = #frac{"+strChi+"}{"+strNdf+"}","");
  //l1 = leg->AddEntry(residPoly2,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{VEM-Charge: "+strVemPk+"}","");
  l1->SetTextColor(kBlue);
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  leg->Draw();

  line = new TLine(rangXmin-5, 0, rangXmax+5, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(vemPk, -6, vemPk, 4);
  line->SetLineColor(kBlue);
  line->SetLineStyle(1);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(rangXmin-5, 1, rangXmax+5, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-5, -1, rangXmax+5, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  //c4->Print("../plots/chargeFitResidualsPoly2863.pdf");
  
}
