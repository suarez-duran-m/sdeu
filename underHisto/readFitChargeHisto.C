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

double fitFunctionCh(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
}


TH1F *getSmooth(TH1F &hist, double xb[])
{
	unsigned int nb = 600;
  double yi = 0.;

	TH1F *hstSmooth = new TH1F("hstSmooth", "", nb, xb);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b > 1 && b<nb-2 )
    {
      yi = hist.GetBinContent(b+1 - 2) 
        + 2*hist.GetBinContent(b+1 - 1) 
        + 3*hist.GetBinContent(b+1) 
        + 2*hist.GetBinContent(b+1 + 1) 
        + hist.GetBinContent(b+1 + 2);
      yi = yi/9.;
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
  TH1F *derihist = new TH1F("derihist", "", nbins, xb);
  double der = 0.;
  for ( int kk=1; kk<nbins-1; kk++ )
  {
    h = ( hist.GetBinCenter( kk+1 ) - hist.GetBinCenter( kk ) );
    der = ( hist.GetBinContent( kk+1 ) -  hist.GetBinContent( kk-1 ) )/ (2.*h);
    derihist->SetBinContent( kk, der );

  }
  return derihist;
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


// ==================================
// *** ***  *** MAIN CODE *** *** *** 

void readFitChargeHisto()
{

  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitChargeHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  //in.open(Form("%skkcharge.dat",dir.Data()));
  in.open(Form("%schargeHist.dat",dir.Data()));

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
  TCanvas *c3 = canvasStyle("c3"); 
  TCanvas *c4 = canvasStyle("c4");

  TH1F *charge = new TH1F("charge", "", nbins, xfadc);
 
  for ( int kk=0; kk<nbins; kk++ )
    charge->SetBinContent(kk, ycnt[kk]);

  TH1F *chargeDer = histDerivative(*charge, xfadc);

  TH1F *chargeSmooth = getSmooth(*charge, xfadc);
  TH1F *chargeSmooDer = histDerivative(*chargeSmooth, xfadc);
  chargeSmooDer->Smooth(500);

  TH1 *test = 0;
  TVirtualFFT::SetTransform(0);
  test = charge->FFT(test, "PH");

  int n = test->GetXaxis()->GetNbins();
  Double_t *re_full = new Double_t[n];
  Double_t *im_full = new Double_t[n];

  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  fft->GetPointsComplex(re_full,im_full);

  int first = 0.03*n;
  for ( int k=first; k<test->GetXaxis()->GetNbins(); k++ )
  {
    re_full[k] = 0;
    im_full[k] = 0;
  }

  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_full,im_full);
  fft_back->Transform();

  TH1 *hbC = 0;
  hbC = TH1::TransformHisto(fft_back,hbC,"Re");

  double xc[601];
  for (int j = 0; j < 402; j++)
    xc[j] = 8.*j;
  for (int j = 0; j < 201; j++)
    xc[400 + j] = 400*8. + 4.*8.*j;

  TH1F *chargeFFT = new TH1F("chargeFFT","", 600, xc);
  for (int j = 0; j < 402; j++)
    chargeFFT->SetBinContent(j + 1, hbC->GetBinContent(j)/600.);
  for (int j = 0; j < 200; j++)
    chargeFFT->SetBinContent(j + 1 + 400, hbC->GetBinContent(j+400)/600.);

  TH1 *chargeFFTDer = histDerivative(*chargeFFT, xfadc);


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

  chargeFFT->SetLineColor(kMagenta+1);
  chargeFFT->SetLineWidth(1);
  chargeFFT->Draw("same");
  
  leg = new TLegend(0.62,0.65,0.95,0.96);
  //leg->SetHeader("#splitline{Charge chargeFFTgrma Station 1740}{(Event 61203949)}");
  leg->SetHeader("#splitline{Charge chargeFFTgrma Station 863}{(Event 61219267)}");
  leg->AddEntry(charge,"Charge chargeFFTgram","f");
  leg->AddEntry(chargeSmooth,"Smooth charge chargeFFTgram","f");
  leg->AddEntry(chargeFFT,"Smooth FFT charge histogram","f");
  leg->Draw();  
  
  c1->Print("../plots/chargeHisto863.pdf");

  c2->cd();
  charge->SetStats(0);
  charge->SetLineColor(kBlue);
  charge->SetLineWidth(1);
  charge->GetXaxis()->SetTitle("[FADC]");
  charge->GetYaxis()->SetTitle("Counts [au]");
  charge->GetXaxis()->SetRangeUser(1300, 1700);
  histoStyle(charge);
  charge->Draw();

  chargeSmooth->SetLineColor(kOrange+10);
  chargeSmooth->SetLineWidth(1);
  chargeSmooth->Draw("same");

  chargeFFT->SetLineColor(kMagenta+1);
  chargeFFT->SetLineWidth(1);
  chargeFFT->Draw("same");
                                                                               
  leg = new TLegend(0.62,0.65,0.95,0.96);
  leg->SetHeader("#splitline{Charge histogrma Station 863}{(Event 61219267)}");
  leg->AddEntry(charge,"Charge histogram","f");
  leg->AddEntry(chargeSmooth,"Smooth charge histogram","f");
  leg->AddEntry(chargeFFT,"Smooth FFT charge histogram","f");
  leg->Draw();  
 
  c2->Print("../plots/chargeHistoZoom863.pdf");

  c3->cd();
  chargeDer->SetStats(0);
  chargeDer->SetLineColor(kBlue);
  chargeDer->SetLineWidth(1);
  chargeDer->GetXaxis()->SetTitle("[FADC]");
  chargeDer->GetYaxis()->SetTitle("[au]");
  chargeDer->GetYaxis()->SetRangeUser(-10,10);
  chargeDer->GetXaxis()->SetRangeUser(0,3000);
  histoStyle(chargeDer);
  chargeDer->Draw();

  line = new TLine(0, 0, 3000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.65,0.7,0.95,0.87);
  leg->SetHeader("From charge histogram","C");
  leg->AddEntry(chargeSmooDer,"First derivative","f");
  leg->Draw();
  c3->Print("../plots/chargeDerHisto863.pdf");


  c4->cd(); 
  chargeSmooDer->SetLineColor(kOrange+10);
  chargeSmooDer->SetLineWidth(1);
  chargeSmooDer->SetStats(0);
  chargeSmooDer->GetYaxis()->SetTitle("[au]");
  chargeSmooDer->GetXaxis()->SetTitle("[FADC]");
  chargeSmooDer->GetYaxis()->SetRangeUser(-10,10);
  chargeSmooDer->GetXaxis()->SetRangeUser(0,3000);
  histoStyle(chargeSmooDer);
  chargeSmooDer->Draw();

  chargeFFTDer->SetLineColor(kMagenta+1);
  chargeFFTDer->SetLineWidth(1);
  chargeFFTDer->Draw("hist same");

  line = new TLine(0, 0, 3000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(1);
  line->Draw();
  
  leg = new TLegend(0.55,0.7,0.95,0.92);
  leg->SetHeader("From Smooth charge histogram","C");
  leg->AddEntry(chargeSmooDer,"First derivative smooth","f");
  leg->AddEntry(chargeFFTDer,"First derivative FFT","f");
  leg->Draw();
  c4->Print("../plots/chargeSmoothDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c5 = canvasStyle("c5");
  TCanvas *c6 = canvasStyle("c6");
  TCanvas *c7 = canvasStyle("c7");

	int rangXmin = 0; // Min for fitting
	int rangXmax = 0; // Max for fitting
	int nXbins = charge->GetXaxis()->GetNbins(); // Number of bins for fitting

  double binMax = 0;
  double binMin = 0;
  for ( int kk=275; kk>125; kk-- ) // from 2496 FADC backward
    if ( chargeSmooDer->GetBinContent( kk ) < 0 )
    {
      if ( binMax < fabs(chargeSmooDer->GetBinContent( kk ) ) )
      {
        binMax = fabs(chargeSmooDer->GetBinContent(kk));
        rangXmax = chargeSmooDer->GetBinCenter(kk);
      }
    }
    else
    {
      binMax = chargeSmooDer->GetBinCenter(kk);
      break;
    }

  rangXmax *= 1.2;


  binMin = 0;
  int tmpneg = 0;
  for ( int kk=25; kk<binMax; kk++ ) // 200 FADC after 0 FADC
  {
    if ( chargeSmooDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( chargeSmooDer->GetBinContent( kk ) < 0 )
      if ( binMin < fabs( chargeSmooDer->GetBinContent( kk ) ) )
      {
        rangXmin = chargeSmooDer->GetBinCenter(kk); // 20 FADC fordward EM-Peak
        binMin = fabs( chargeSmooDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  }

  //rangXmin *= 2.5;
  
	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( charge->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( charge->GetBinCenter(b+1) );
	}


	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76);
  chFit->Fit("fitFcn","QR");

  cerr << rangXmin << " " << rangXmax << " " << binMax << endl;
  cerr << fitFcn->GetChisquare() << " " << fitFcn->GetNDF() << endl;

  TF1 *expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  TF1 *lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);

  double chargeVal = 0.;
  chargeVal = chargeMaxPk(fitFcn);

  expon->SetParameter(0, fitFcn->GetParameter(3));
  expon->SetParameter(1, fitFcn->GetParameter(4));

  lognorm->SetParameter(0, fitFcn->GetParameter(0));
  lognorm->SetParameter(1, fitFcn->GetParameter(1));
  lognorm->SetParameter(2, fitFcn->GetParameter(2));

  // ===========================
  // *** Computing Residuals ***

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;
  double tmp = 0.;
  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( charge->GetBinCenter(kbin) >= rangXmin && charge->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( charge->GetBinCenter(kbin) );
      tmp = fitFcn->Eval( charge->GetBinCenter(kbin) ) - charge->GetBinContent(kbin);
      yResid.push_back( tmp );
      errResid.push_back( sqrt(
            pow(sqrt( charge->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(fitFcn->Eval( charge->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraph = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *logNormResid = new TH1F("logNormResid", "", 601/20., -300, 300);
  for ( int val=0; val<yResid.size(); val++ )
    logNormResid->Fill( yResid[val] );

  int nPoints = chFit->GetN();
  int nbins2 = ( chFit->GetPointX(nPoints - 1) - chFit->GetPointX(0) ) / 8.;
  TH1F *tmp2 = new TH1F("tmp2", "", nbins2, chFit->GetPointX(0), chFit->GetPointX(nPoints-1));

  double x, y;
  for(int i=0; i < nPoints; i++ )
  {
    chFit->GetPoint(i, x, y);
    tmp2->SetBinContent(i, y);
  }
 
  c5->cd();
  TPaveStats *ptstats;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  chFit->GetListOfFunctions()->Add(ptstats);
  chFit->SetTitle("");
  chFit->GetXaxis()->SetRangeUser(100, 3000);
  chFit->GetYaxis()->SetRangeUser(0, 1600);
  chFit->SetLineColor(kBlue);
  chFit->SetLineWidth(1);
  chFit->GetXaxis()->SetTitle("[FADC]");
  chFit->GetYaxis()->SetTitle("Counts [au]");
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

  c5->Print("../plots/chargeFittedHisto863.pdf");

  /*
  c6->cd();
  TRatioPlot *rpmean;
  TF1 *fitmean21 = chFit->GetFunction("fitFcn");
  tmp2->Fit(fitmean21,"QR");
  tmp2->GetXaxis()->SetRangeUser(100, 2200);
  rpmean = new TRatioPlot(tmp2);
  rpmean->SetGraphDrawOpt("APL*");
  rpmean->Draw();
  rpmean->GetLowerRefYaxis()->SetRangeUser(-5, 5);
  c6->Update();
  */

  c7->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(100, 3000);
  residGraph->GetYaxis()->SetRangeUser(-120, 160);
  residGraph->SetLineColor(kBlue);
  residGraph->SetLineWidth(1);
  residGraph->SetMarkerSize(1.5);
  residGraph->SetMarkerStyle(20);
  residGraph->GetXaxis()->SetTitle("[FADC * 8.33 ns]");
  residGraph->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residGraph);
  residGraph->Draw("APL");

  TString chindf;
  TString valcharge;
  chindf.Form("%.2f", fitFcn->GetChisquare() / fitFcn->GetNDF() );
  valcharge.Form("%.2f", chargeVal);

  leg = new TLegend(0.58,0.75,0.96,0.97);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valcharge+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(400, 0, 1800, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(chargeVal, -120, chargeVal, 160);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c7->Print("../plots/chargeFitResiduals863.pdf");


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
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( charge->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2->Eval( charge->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraphPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *poly2Resid = new TH1F("poly2Resid", "", 601/20., -300, 300);
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
  residGraphPoly2->GetYaxis()->SetRangeUser(-120, 160);
  residGraphPoly2->SetLineColor(kBlue);
  residGraphPoly2->SetLineWidth(1);
  residGraphPoly2->SetMarkerSize(1.5);
  residGraphPoly2->SetMarkerStyle(20);
  residGraphPoly2->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residGraphPoly2);
  residGraphPoly2->Draw("APL");

  chindf;
  valcharge;
  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
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
  logNormResid->GetXaxis()->SetTitle("y_{fit} - y_{data} [au]");
  logNormResid->GetYaxis()->SetTitle("Counts [au]");
  logNormResid->GetYaxis()->SetRangeUser(0, 54);
  logNormResid->GetXaxis()->SetRangeUser(-100, 100);
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

}
