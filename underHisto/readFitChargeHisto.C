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

double fitFunctionCh(double *x0, double *par) 
{
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
}


TH1F *getSmooth(TH1F &hist, double xb[])
{
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


TH1F *histDerivative(TH1F &hist, double xb[]) // Central differences
{
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


vector < double > getFitRange( TH1F &h )
{
  vector < double > minmax;
  int minRng = 0; // Min for fitting
	int maxRng = 0; // Max for fitting

  double binMin = 0.;
  double binMax = 0.;
  double rawbinMax = 0.;
  double rawbinMin = 0.;
  double vemDer = 0.;

  for ( int kk=275; kk>125; kk-- ) // from ~2000 FADC backward
    if ( h.GetBinContent(kk) < 0 )
    {
      if ( binMax < fabs( h.GetBinContent( kk ) ) )
      {
        binMax = fabs( h.GetBinContent(kk));
        maxRng = h.GetBinCenter(kk); // Change of concavity
      }
    }
    else
    {
      binMax = h.GetBinCenter(kk); // tmp FADC for VEM
      rawbinMax = kk; // Bin for tmp VEM
      vemDer = h.GetBinCenter(kk) + 4.;
      // +4 because the zero is between this pixel and
      // the next one.
      break;
    }

  binMin = 0;
  for ( int kk=rawbinMax; kk>0; kk-- )
    if ( h.GetBinContent(kk) > 0 )
    {
      if ( binMin < h.GetBinContent(kk) )
      {
        binMin = h.GetBinContent(kk); 
        minRng = h.GetBinCenter(kk); // Change of concavity
      }
    }
    else
    {
      binMin =  h.GetBinCenter(kk);
      minRng = binMin;
      rawbinMin = kk;
      break;
    }

  minmax.push_back( minRng );
  minmax.push_back( maxRng );
  minmax.push_back( binMax );
  minmax.push_back( rawbinMin );
  minmax.push_back( vemDer );

  return minmax;
}


double chargeMaxPk(TF1* f) 
{
	TCanvas* c = new TCanvas("c","c",800,600);
  TGraph* der = (TGraph*)f->DrawDerivative("goff");
  double * x = der->GetX();
  double * yd1 = der->GetY();
  int n = der->GetN();
  double d0x = 0;
  double d0y = 1000;
  for (int j = 1; j < n; ++j)
  {
    const double yd2 = f->Derivative2(x[j]);
    if (yd2 <0 && d0y > fabs(yd1[j]))
    {
      d0y = fabs(yd1[j]);
      d0x = x[j];
    }
  }
  delete der;
  delete c;
  return d0x;
}


void getResiduals( TH1F *hist, TF1 *func,
    double rangMin, double rangMax,
    vector < double > &x, vector < double > &y, vector < double > &err )
{
  int nbins = hist->GetXaxis()->GetNbins();
  double tmp = 0.;
  for ( int kbin=1; kbin<nbins; kbin++ )
    if ( hist->GetBinCenter(kbin) >= rangMin && hist->GetBinCenter(kbin) <= rangMax )
    {
      x.push_back( hist->GetBinCenter(kbin) );
      tmp = func->Eval( hist->GetBinCenter(kbin) ) - hist->GetBinContent(kbin);
      y.push_back( tmp / sqrt( hist->GetBinContent(kbin) ) );
      err.push_back( sqrt( pow(sqrt( hist->GetBinContent(kbin) ),2)
            + pow(sqrt( sqrt(func->Eval( hist->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( hist->GetBinContent(kbin) ) );
    }
}


// ==================================
// *** ***  *** MAIN CODE *** *** *** 

void readFitChargeHisto()
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitChargeHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  //in.open(Form("%skkcharge.dat",dir.Data()));
  in.open(Form("%schargeHist61219267.dat",dir.Data()));

  const int nbins = 600;

  double tmpXb;
  double tmpXfadc;
  double tmpYcnt;

  double xbin[nbins+1];
  double xfadc[nbins+1];
  double ycnt[nbins+1];

  int bincnt = 0;

  while ( bincnt < nbins+1 ) 
  {
    in >> tmpXb >> tmpXfadc >> tmpYcnt;
    xbin[bincnt] = tmpXb;
    xfadc[bincnt] = tmpXfadc+4; // Correct for binCenter
    ycnt[bincnt] = tmpYcnt;
    bincnt++;
  }
  in.close();

  TCanvas *c1 = canvasStyle("c1"); 
  TCanvas *c2 = canvasStyle("c2"); 

  TH1F *charge = new TH1F("charge", "", nbins, xfadc);
 
  for ( int kk=0; kk<nbins; kk++ )
    charge->SetBinContent(kk, ycnt[kk]);

  TH1F *chargeDer = histDerivative(*charge, xfadc);

  TH1F *chargeSmooth = getSmooth(*charge, xfadc);
  TH1F *chargeSmooDer = histDerivative(*chargeSmooth, xfadc);
  TH1F *chargeSmooDerSmth = getSmooth(*chargeSmooDer, xfadc);
  chargeSmooDerSmth = getSmooth(*chargeSmooDerSmth, xfadc);

  TLine *line;
  TLegend *leg;

  c1->cd();
  charge->SetStats(0);
  charge->SetLineColor(kBlue);
  charge->SetLineWidth(1);
  charge->GetXaxis()->SetTitle("[FADC]");
  charge->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(charge);
  charge->Draw();
  
  chargeSmooth->SetLineColor(kOrange+10);
  chargeSmooth->SetLineWidth(1);
  chargeSmooth->Draw("same"); 
  
  leg = new TLegend(0.62,0.65,0.95,0.96);
  leg->SetHeader("#splitline{Charge histogram Station 863}{(Event 61219267)}");
  leg->AddEntry(charge,"Charge histogram","f");
  leg->AddEntry(chargeSmooth,"Smooth charge histogram ","f");
  leg->Draw();   
  c1->Print("../plots/chargeHisto863.pdf");

  c2->cd();
  chargeDer->SetStats(0);
  chargeDer->SetLineColor(kBlack);
  chargeDer->SetLineWidth(1);
  chargeDer->GetXaxis()->SetTitle("[FADC]");
  chargeDer->GetYaxis()->SetTitle("[au]");
  chargeDer->GetYaxis()->SetRangeUser(-3,4);
  chargeDer->GetXaxis()->SetRangeUser(0,3000);
  histoStyle(chargeDer);
  chargeDer->Draw();

  chargeSmooDer->SetStats(0);
  chargeSmooDer->SetLineColor(kOrange+9);
  chargeSmooDer->SetLineWidth(1);
  chargeSmooDer->Draw("sames");

  chargeSmooDerSmth->SetStats(0);
  chargeSmooDerSmth->SetLineColor(kGreen+3);
  chargeSmooDerSmth->SetLineWidth(2);
  chargeSmooDerSmth->Draw("sames");

  line = new TLine(0, 0, 3000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.53,0.78,0.96,0.97);
  leg->SetHeader("From charge histogram","C");
  leg->AddEntry(chargeDer,"First derivative Charge-H.","l");
  leg->AddEntry(chargeSmooDer,"First derivative Smooth-Charge-H.","l");
  leg->AddEntry(chargeSmooDerSmth,"First derivative Smooth-Charge-H.","l");

  leg->SetTextSize(0.04);
  leg->Draw();
  c2->Print("../plots/chargeDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c3 = canvasStyle("c3");
  TCanvas *c4 = canvasStyle("c4");

	int rangXmin = 0; // Min for fitting
	int rangXmax = 0; // Max for fitting
  double chargeVal = 0.;
	int nXbins = charge->GetXaxis()->GetNbins(); // Number of bins for fitting

  vector < double > fitRange = getFitRange(*chargeSmooDerSmth);
  rangXmin = fitRange[0];
  rangXmax = fitRange[1];
  
	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( charge->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( charge->GetBinCenter(b+1) );
	}

	TGraphErrors* chToFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2;

  // ===========================
  // *** Computing Residuals ***

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;
  double tmp = 0.;
  int bigRsd = 0;
  double reduceFactor = 0.05;
  double chi2Ndf = 500.; 
  double bestXmin = 0.;
  double bestXmax = 0.;

  for ( int nl=0; nl<10; nl++ ) 
  {
    bigRsd = 0;
    tmp = 0;

    poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
    chToFit->Fit("poly2","QR");

    xResid.clear();
    yResid.clear();
    errResid.clear();

    for ( int kbin=1; kbin<nXbins; kbin++ )
      if ( charge->GetBinCenter(kbin) >= rangXmin && charge->GetBinCenter(kbin) <= rangXmax )
      {
        xResid.push_back( charge->GetBinCenter(kbin) );
        tmp = poly2->Eval( charge->GetBinCenter(kbin) ) - charge->GetBinContent(kbin);
        yResid.push_back( tmp / sqrt( charge->GetBinContent(kbin) ) );
        errResid.push_back( sqrt(
              pow(sqrt( charge->GetBinContent(kbin) ),2) 
              + pow(sqrt( sqrt(poly2->Eval( charge->GetBinCenter(kbin) ) ) ),2)
              ) / sqrt( charge->GetBinContent(kbin) ) );
      }

    for ( int rsd=0; rsd<yResid.size(); rsd++ )
      if ( fabs ( yResid[rsd] ) > 1 )
        bigRsd++;

    if ( chi2Ndf > poly2->GetChisquare() / poly2->GetNDF() )
    {
      chi2Ndf = poly2->GetChisquare() / poly2->GetNDF();
      bestXmin = rangXmin;
      bestXmax = rangXmax;
    }
    if ( bigRsd < 0.3*yResid.size() )
      break;
    else
    {
      rangXmin = (1.+reduceFactor)*fitRange[0];
      rangXmax = (1.-reduceFactor)*fitRange[1];
      reduceFactor += 0.01;
    }
  }

  if ( chi2Ndf < poly2->GetChisquare() / poly2->GetNDF() )
  {
    rangXmin = bestXmin;
    rangXmax = bestXmax;
    poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
    chToFit->Fit("poly2","QR");

    xResid.clear();
    yResid.clear();
    errResid.clear();
    tmp = 0;

    for ( int kbin=1; kbin<nXbins; kbin++ )
      if ( charge->GetBinCenter(kbin) >= rangXmin && charge->GetBinCenter(kbin) <= rangXmax )
      {
        xResid.push_back( charge->GetBinCenter(kbin) );
        tmp = poly2->Eval( charge->GetBinCenter(kbin) ) - charge->GetBinContent(kbin);
        yResid.push_back( tmp / sqrt( charge->GetBinContent(kbin) ) );
        errResid.push_back( sqrt(
              pow(sqrt( charge->GetBinContent(kbin) ),2) 
              + pow(sqrt( sqrt(poly2->Eval( charge->GetBinCenter(kbin) ) ) ),2)
              ) / sqrt( charge->GetBinContent(kbin) ) );
      }
  }
  chargeVal = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0));

  cerr << rangXmin << " " << rangXmax << " " << chargeVal << endl;
  cerr << poly2->GetChisquare() << " " 
    << poly2->GetNDF() << " " 
    << poly2->GetChisquare() / poly2->GetNDF()
    << endl;

  TGraphErrors* residGraph = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  c3->cd();
  TPaveStats *ptstats;
  //gStyle->SetOptStat();
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  ptstats->SetFillColor(0);
  chToFit->GetListOfFunctions()->Add(ptstats);
  chToFit->SetTitle("");
  chToFit->GetXaxis()->SetRangeUser(100, 3000);
  chToFit->GetYaxis()->SetRangeUser(100, 1300);
  chToFit->SetLineColor(kBlue);
  chToFit->SetLineWidth(1);
  chToFit->GetXaxis()->SetTitle("[FADC]");
  chToFit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(chToFit);
  chToFit->Draw();

  poly2->SetLineWidth(4);
  poly2->Draw("same");

  line = new TLine(chargeVal, 1e2, chargeVal, 1300);
  line->SetLineStyle(4);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  line->Draw();

  line = new TLine(fitRange[4], 1e2, fitRange[4], 1300);
  line->SetLineStyle(4);
  line->SetLineColor(kGreen+3);
  line->SetLineWidth(3);
  line->Draw();
  c3->Print("../plots/chargeFittedHisto863.pdf");


  c4->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(rangXmin-50, rangXmax+50);
  residGraph->GetYaxis()->SetRangeUser(-9, 4.);
  residGraph->SetLineColor(kBlack);
  residGraph->SetLineWidth(1);
  residGraph->SetMarkerSize(1.5);
  residGraph->SetMarkerStyle(20);
  residGraph->GetXaxis()->SetTitle("[FADC]");
  residGraph->GetYaxis()->SetTitle("Residuals");
  residGraph->SetMarkerSize(3);
  residGraph->SetMarkerColor(kBlack);
  residGraph->SetMarkerStyle(23);
  histoStyle(residGraph);
  residGraph->GetYaxis()->SetTitleOffset(0.8);
  residGraph->Draw("APL");

  TString chindf;
  TString valcharge;
  TString valchargeDer;
  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
  valcharge.Form("%.2f", chargeVal);
  valchargeDer.Form("%.2f", fitRange[4]);


  leg = new TLegend(0.39,0.16,0.66,0.50);
  //leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valcharge+"}","");
  leg->AddEntry(residGraph,"#frac{#chi^{2}}{ndf} = "+chindf,"");
  leg->AddEntry(residGraph,"Peak val. from fit: "+valcharge,"");
  leg->AddEntry(residGraph,"Peak val. from BXL-method: "+valchargeDer,"");
  leg->SetLineWidth(0);
  leg->SetTextSize(0.06);
  leg->Draw();

  line = new TLine(rangXmin-50, 0, rangXmax+50, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chargeVal, -9, chargeVal, 4.);
  line->SetLineStyle(4);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);
  line->Draw();

  line = new TLine(fitRange[4], -9, fitRange[4], 4.);
  line->SetLineStyle(4);
  line->SetLineColor(kGreen+3);
  line->SetLineWidth(3);
  line->Draw();

  line = new TLine(rangXmin-50, 1, rangXmax+50, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-50, -1, rangXmax+50, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c4->Print("../plots/chargeFitResiduals863.pdf");

/*
  // =============================
  // *** *** Fitting Poly2 *** ***

  TCanvas *c8 = canvasStyle("c8");
  TCanvas *c9 = canvasStyle("c9");

  TGraphErrors* poly2Fit = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  rangXmin = binMax*0.8;
  rangXmax = binMax*1.3;

  TF1 *poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2Fit->Fit("poly2","QR");

  chargeVal = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0));

  // =====================================
  // *** Computing residuals for Poly2 ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0;

  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( charge->GetBinCenter(kbin) >= rangXmin && charge->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( charge->GetBinCenter(kbin) );
      tmp = poly2->Eval( charge->GetBinCenter(kbin) ) - charge->GetBinContent(kbin); 
      yResid.push_back( tmp / sqrt( charge->GetBinContent(kbin) ) );
      errResid.push_back( sqrt( 
            pow(sqrt( charge->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2->Eval( charge->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( charge->GetBinContent(kbin) ) );
    }

  TGraphErrors* residGraphPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *poly2Resid = new TH1F("poly2Resid", "", 11, -5, 5);
  for ( int val=0; val<yResid.size(); val++ )
    poly2Resid->Fill( yResid[val] );


  c8->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2Fit->GetListOfFunctions()->Add(ptstats);
  poly2Fit->SetTitle("");
  poly2Fit->GetXaxis()->SetRangeUser(100, 3000);
  poly2Fit->GetYaxis()->SetRangeUser(0, 1600);
  poly2Fit->SetLineColor(kBlue);
  poly2Fit->SetLineWidth(1);
  poly2Fit->GetXaxis()->SetTitle("[FADC]");
  poly2Fit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2Fit);
  poly2Fit->Draw();

  line = new TLine(chargeVal, 0, chargeVal, 1600);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c8->Print("../plots/chargeFitPoly2863.pdf");

  c9->cd();
  residGraphPoly2->SetTitle("");
  residGraphPoly2->GetYaxis()->SetRangeUser(-7, 4.2);
  residGraphPoly2->SetLineColor(kBlue);
  residGraphPoly2->SetLineWidth(1);
  residGraphPoly2->SetMarkerSize(1.5);
  residGraphPoly2->SetMarkerStyle(20);
  residGraphPoly2->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2->GetYaxis()->SetTitle("Residuals [au]");
  histoStyle(residGraphPoly2);
  residGraphPoly2->GetYaxis()->SetTitleOffset(0.8);
  residGraphPoly2->Draw("APL");

  chindf;
  valcharge;
  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
  valcharge.Form("%.2f", chargeVal);

  leg = new TLegend(0.58,0.19,0.96,0.40);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valcharge+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(rangXmin-50, 0, rangXmax+53, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chargeVal, -120, chargeVal, 160);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(rangXmin-50, 1, rangXmax+53, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-50, -1, rangXmax+53, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c9->Print("../plots/chargeFitResidualsPoly2863.pdf");

  // ============================================
  // *** *** Aplying for Smooth Histogram *** ***

  TCanvas *c10 = canvasStyle("c10");
  TCanvas *c11 = canvasStyle("c11");

  xbins.clear();
  ycnts.clear();
  yerrs.clear();
	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( chargeFFT->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( chargeFFT->GetBinCenter(b+1) );
	}

  binMax = 0;
  binMin = 0;
  for ( int kk=312; kk>125; kk-- ) // from 2960 FADC backward
    if ( chargeFFTDer->GetBinContent( kk ) < 0 )
    {
      if ( binMax < fabs(chargeFFTDer->GetBinContent( kk ) ) )
      {
        binMax = fabs(chargeFFTDer->GetBinContent(kk));
        rangXmax = chargeFFTDer->GetBinCenter(kk);
      }
    }
    else
    {
      binMax = chargeFFTDer->GetBinCenter(kk);
      break;
    }

  binMin = 0;
  tmpneg = 0;
  for ( int kk=25; kk<binMax; kk++ ) // 200 FADC after 0 FADC
  {
    if ( chargeFFTDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( chargeFFTDer->GetBinContent( kk ) < 0 )
      if ( binMin < fabs( chargeFFTDer->GetBinContent( kk ) ) )
      {
        rangXmin = chargeFFTDer->GetBinCenter(kk);
        binMin = fabs( chargeFFTDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  }
 
  chFit = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76);
  chFit->Fit("fitFcn","QR");

  expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);

  chargeVal = chargeMaxPk(fitFcn);

  expon->SetParameter(0, fitFcn->GetParameter(3));
  expon->SetParameter(1, fitFcn->GetParameter(4));

  lognorm->SetParameter(0, fitFcn->GetParameter(0));
  lognorm->SetParameter(1, fitFcn->GetParameter(1));
  lognorm->SetParameter(2, fitFcn->GetParameter(2));

  // ===========================================
  // *** Computing residuals for FFT LogNorm ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0;

  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( chargeFFT->GetBinCenter(kbin) >= rangXmin && chargeFFT->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( chargeFFT->GetBinCenter(kbin) );
      tmp = fitFcn->Eval( chargeFFT->GetBinCenter(kbin) ) - chargeFFT->GetBinContent(kbin); 
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( chargeFFT->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(fitFcn->Eval( chargeFFT->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraphLogNormSmooth = new TGraphErrors( xResid.size(), &xResid.front(), &yResid.front(), 0, &errResid.front() );

  TH1F *logNormResidSmooth = new TH1F("logNormResidSmooth", "", 100, -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    logNormResidSmooth->Fill( yResid[val] );

  c10->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  chFit->GetListOfFunctions()->Add(ptstats);
  chFit->SetTitle("");
  chFit->SetLineColor(kBlue);
  chFit->SetLineWidth(1);
  chFit->GetXaxis()->SetTitle("[FADC]");
  chFit->GetYaxis()->SetTitle("Counts [au]");
  chFit->GetXaxis()->SetRangeUser(100, 3000);
  chFit->GetYaxis()->SetRangeUser(0, 1600);
  histoStyle(chFit);
  chFit->Draw();

  expon->SetLineColor(kGreen+2);
  expon->Draw("same");
  lognorm->SetLineColor(kMagenta+2);
  lognorm->Draw("same");

  leg = new TLegend(0.25,0.67,0.61,0.97);
  leg->AddEntry(charge,"Charge histogram","f");
  leg->AddEntry(fitFcn,"Fit Exp+Lognormal","f");
  leg->AddEntry(expon,"Fit Exp","f");
  leg->AddEntry(lognorm,"Fit Lognormal","f");
  leg->Draw();

  line = new TLine(chargeVal, 0, chargeVal, 1600);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c10->Print("../plots/chargeFFTFitLogNorm863.pdf");

  c11->cd();
  residGraphLogNormSmooth->SetTitle("");
  residGraphLogNormSmooth->SetLineColor(kBlue);
  residGraphLogNormSmooth->SetLineWidth(1);
  residGraphLogNormSmooth->GetXaxis()->SetTitle("[FADC]");
  residGraphLogNormSmooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residGraphLogNormSmooth->GetYaxis()->SetRangeUser(-120,160);
  residGraphLogNormSmooth->SetMarkerSize(1.5);
  residGraphLogNormSmooth->SetMarkerStyle(20);
  histoStyle(residGraphLogNormSmooth);
  residGraphLogNormSmooth->Draw("APL");

  chindf;
  valcharge;
  chindf.Form("%.2f", fitFcn->GetChisquare() / fitFcn->GetNDF() );
  valcharge.Form("%.2f", chargeVal);

  leg = new TLegend(0.58,0.75,0.96,0.97);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valcharge+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(rangXmin-50, 0, rangXmax+50, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chargeVal, -120, chargeVal, 160);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c11->Print("../plots/chargeFFTLogNormResiduals863.pdf");

  // =====================
  // *** For TTF Poly2 ***
  TCanvas *c12 = canvasStyle("c12");
  TCanvas *c13 = canvasStyle("c13");

  xbins.clear();
  ycnts.clear();
  yerrs.clear();
	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( chargeFFT->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( chargeFFT->GetBinCenter(b+1) );
	}

  rangXmin = binMax*0.8;
  rangXmax = binMax*1.3;

  TGraphErrors* poly2FitSmooth = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2Smooth = new TF1("poly2Smooth","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2FitSmooth->Fit("poly2Smooth","QR");

  chargeVal = -poly2Smooth->GetParameter(1) / (2.*poly2Smooth->GetParameter(0));

  // =====================================
  // *** Computing residuals for Poly2 ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0;

  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( chargeFFT->GetBinCenter(kbin) >= rangXmin && chargeFFT->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( chargeFFT->GetBinCenter(kbin) );
      tmp = poly2Smooth->Eval( chargeFFT->GetBinCenter(kbin) ) - chargeFFT->GetBinContent(kbin); 
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( chargeFFT->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2Smooth->Eval( chargeFFT->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraphPoly2Smooth = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *residPoly2Smooth = new TH1F("residPoly2Smooth", "", 101, -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    residPoly2Smooth->Fill( yResid[val] );

  c12->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2FitSmooth->GetListOfFunctions()->Add(ptstats);
  poly2FitSmooth->SetTitle("");
  poly2FitSmooth->SetLineColor(kBlue);
  poly2FitSmooth->SetLineWidth(1);
  poly2FitSmooth->GetXaxis()->SetTitle("[FADC]");
  poly2FitSmooth->GetYaxis()->SetTitle("Counts [au]");
  poly2FitSmooth->GetXaxis()->SetRangeUser(100, 3000);
  poly2FitSmooth->GetYaxis()->SetRangeUser(0, 1600);
  histoStyle(poly2FitSmooth);
  poly2FitSmooth->Draw();

  line = new TLine(chargeVal, 0, chargeVal, 1600);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c12->Print("../plots/chargeFFTFitPoly2863.pdf");

  c13->cd();
  residGraphPoly2Smooth->SetTitle("");
  residGraphPoly2Smooth->SetLineColor(kBlue);
  residGraphPoly2Smooth->SetLineWidth(1);
  residGraphPoly2Smooth->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2Smooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residGraphPoly2Smooth->GetYaxis()->SetRangeUser(-120,160);
  residGraphPoly2Smooth->SetMarkerSize(1.5);
  residGraphPoly2Smooth->SetMarkerStyle(20);
  histoStyle(residGraphPoly2Smooth);
  residGraphPoly2Smooth->Draw("APL");

  chindf;
  valcharge;
  chindf.Form("%.2f", poly2Smooth->GetChisquare() / poly2Smooth->GetNDF() );
  valcharge.Form("%.2f", chargeVal);

  leg = new TLegend(0.58,0.75,0.96,0.97);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valcharge+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(rangXmin-50, 0, rangXmax+50, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chargeVal, -120, chargeVal, 160);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c13->Print("../plots/chargeFitFFTResidualsPoly2863.pdf");


  // ============================================================
  // *** *** *** Plotting for residuals distributions *** *** ***

  TCanvas *c14 = canvasStyle("c14");
  TCanvas *c15 = canvasStyle("c15");

  c14->cd();
  logNormResid->SetStats(1);
  logNormResid->SetTitle("");
  logNormResid->SetLineColor(kBlack);
  logNormResid->SetLineWidth(1);
  logNormResid->SetFillColor(kBlack);
  logNormResid->SetFillStyle(3001);
  logNormResid->GetXaxis()->SetTitle("Residuals [au]");
  logNormResid->GetYaxis()->SetTitle("Counts [au]");
  logNormResid->GetYaxis()->SetRangeUser(0, 82);
  logNormResid->GetXaxis()->SetRangeUser(-4, 4);
  histoStyle(logNormResid);
  logNormResid->Draw();

  poly2Resid->SetTitle("");
  poly2Resid->SetLineColor(kBlue);
  poly2Resid->SetFillColor(kBlue);
  poly2Resid->SetFillStyle(3001); 
  poly2Resid->SetLineWidth(1);
  poly2Resid->Draw("sames");

  ptstats = new TPaveStats(0.13, 0.7, 0.36, 0.97,"brNDC");
  ptstats->SetTextColor(kBlack);
  logNormResid->SetName("Resid. LogNorm");
  logNormResid->GetListOfFunctions()->Add(ptstats);

  ptstats = new TPaveStats(0.72, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  poly2Resid->SetName("Resid. Poly2");
  poly2Resid->GetListOfFunctions()->Add(ptstats);

  c14->Print("../plots/chargeResidualsDist863.pdf");

  c15->cd();
  logNormResidSmooth->SetStats(1);
  logNormResidSmooth->SetTitle("");
  logNormResidSmooth->SetLineColor(kBlack);
  logNormResidSmooth->SetLineWidth(1);
  logNormResidSmooth->SetFillColor(kBlack);
  logNormResidSmooth->SetFillStyle(3001);
  logNormResidSmooth->GetXaxis()->SetTitle("y_{fit} - y_{data} [au]");
  logNormResidSmooth->GetYaxis()->SetTitle("Counts [au]");
  logNormResidSmooth->GetYaxis()->SetRangeUser(0, 25);
  logNormResidSmooth->GetXaxis()->SetRangeUser(-25, 25);
  histoStyle(logNormResidSmooth);
  logNormResidSmooth->Draw();

  residPoly2Smooth->SetTitle("");
  residPoly2Smooth->SetLineColor(kBlue);
  residPoly2Smooth->SetFillColor(kBlue);
  residPoly2Smooth->SetFillStyle(3001); 
  residPoly2Smooth->SetLineWidth(1);
  residPoly2Smooth->Draw("sames");

  ptstats = new TPaveStats(0.13, 0.7, 0.36, 0.97,"brNDC");
  ptstats->SetTextColor(kBlack);
  logNormResidSmooth->SetName("Resid. FFT LogNorm");
  logNormResidSmooth->GetListOfFunctions()->Add(ptstats);

  ptstats = new TPaveStats(0.72, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  residPoly2Smooth->SetName("Resid. FFT Poly2");
  residPoly2Smooth->GetListOfFunctions()->Add(ptstats);

  c15->Print("../plots/chargeFFTResidualsDist863.pdf");
*/
}
