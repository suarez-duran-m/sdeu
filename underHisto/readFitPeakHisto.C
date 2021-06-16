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

double fitFunctionPk(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
}


TH1F *getSmooth(TH1F &hist, double xb[])
{
	unsigned int nb = 150;
  double yi = 0.;

	TH1F *hstSmooth = new TH1F("hstSmooth", "", nb, xb);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b > 1 && b<148 )
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
  int nbins = 150;
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


double peakMaxPk(TF1* f) 
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

void readFitPeakHisto()
{

  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitPeakHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  //in.open(Form("%skkpeak.dat",dir.Data()));
  in.open(Form("%speakHist.dat",dir.Data()));
  
  srand (time(NULL)); 
  const int nbins = 150;

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
    xfadc[bincnt] = tmpXfadc+2; // Correct for binCenter
    ycnt[bincnt] = tmpYcnt;
    bincnt++;
  }
  in.close();

  TCanvas *c1 = canvasStyle("c1"); 
  TCanvas *c2 = canvasStyle("c2"); 
  TCanvas *c3 = canvasStyle("c3");

  TH1F *peak = new TH1F("peak", "", nbins, xfadc);
  
  for ( int kk=0; kk<nbins; kk++ )
    peak->SetBinContent(kk, ycnt[kk]);

  TH1F *peakDer = histDerivative(*peak, xfadc);

  TH1F *peakSmooth = getSmooth(*peak, xfadc);
  TH1F *peakSmooDer = histDerivative(*peakSmooth, xfadc);

  TH1 *test = 0;
  TVirtualFFT::SetTransform(0);
  test = peak->FFT(test, "PH");

  int n = test->GetXaxis()->GetNbins();
  Double_t *re_full = new Double_t[n];
  Double_t *im_full = new Double_t[n];

  TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
  fft->GetPointsComplex(re_full,im_full);

  int first = 0.06*n;
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
  hbC->SetTitle("The backward transform result");

  double xc[151];
  for (int j = 0; j < 102; j++)
    xc[j] = 4.*j - 12;
  for (int j = 0; j < 51; j++)
    xc[100 + j] = 100*4. + 4.*4.*j;

  TH1F *peakFFT = new TH1F("peakFFT","", 150, xc);
  for (int j = 0; j < 102; j++)
    peakFFT->SetBinContent(j + 1, hbC->GetBinContent(j)/150.);
  for (int j = 0; j < 50; j++)
    peakFFT->SetBinContent(j + 1 + 100, hbC->GetBinContent(j+100)/150.);

  TH1 *peakFFTDer = histDerivative(*peakFFT, xfadc);

  TLine *line;
  TLegend *leg;

  c1->cd();
  peak->SetStats(0);
  peak->SetLineColor(kBlue);
  peak->GetXaxis()->SetRangeUser(-100, 1000);
  peak->SetLineWidth(1);
  peak->GetXaxis()->SetTitle("[FADC]");
  peak->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(peak);
  peak->Draw();

  peakSmooth->SetLineColor(kOrange+10);
  peakSmooth->SetLineWidth(1);
  peakSmooth->Draw("same");

  peakFFT->SetLineColor(kMagenta+1);
  peakFFT->SetLineWidth(1);
  peakFFT->Draw("same");

  leg = new TLegend(0.62,0.65,0.95,0.96);
  //leg->SetHeader("#splitline{Peak histogrma Station 1740}{(Event 61509855)}");
  leg->SetHeader("#splitline{Peak histogrma Station 863}{(Event 61219267)}");
  leg->AddEntry(peak,"Peak histogram","f");
  leg->AddEntry(peakSmooth,"Smooth peak histogram","f");
  leg->AddEntry(peakFFTDer,"First derivative FFT","f");
  leg->Draw();  
  c1->Print("../plots/peakHisto863.pdf");

  c2->cd();
  peakDer->SetStats(0);
  peakDer->SetLineColor(kBlue);
  peakDer->SetLineWidth(1);
  peakDer->GetXaxis()->SetRangeUser(-10, 1000);
  peakDer->GetXaxis()->SetTitle("[FADC]");
  peakDer->GetYaxis()->SetRangeUser(-50, 50);
  peakDer->GetYaxis()->SetTitle("[au]");
  histoStyle(peakDer);
  peakDer->Draw();

  line = new TLine(0, 0, 1000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.65,0.7,0.95,0.87);
  leg->SetHeader("From peak histogram","C");
  leg->AddEntry(peakSmooDer,"First derivative","f");
  leg->Draw();
  c2->Print("../plots/peakDerHisto863.pdf");

  c3->cd();  
  peakSmooDer->SetLineColor(kOrange+10);
  peakSmooDer->SetLineWidth(1);
  peakSmooDer->SetStats(0);
  peakSmooDer->GetXaxis()->SetRangeUser(-10, 1000);
  peakSmooDer->GetXaxis()->SetTitle("[FADC]");
  peakSmooDer->GetYaxis()->SetRangeUser(-50, 50);
  peakSmooDer->GetYaxis()->SetTitle("[au]");
  histoStyle(peakSmooDer);
  peakSmooDer->Draw();

  peakFFTDer->SetLineColor(kMagenta+1);
  peakFFTDer->SetLineWidth(2);
  peakFFTDer->Draw("hist same");

  line = new TLine(0, 0, 1000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.5,0.64,0.95,0.9);
  leg->SetHeader("From Smooth peak histogram","C");
  leg->AddEntry(peakSmooDer,"First derivative smooth","f");
  leg->AddEntry(peakFFTDer,"First derivative FFT","f");
  leg->Draw();
  c3->Print("../plots/peakSmoothDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c4 = canvasStyle("c4");
  //TCanvas *c5 = canvasStyle("c5");
  TCanvas *c6 = canvasStyle("c6");

	int rangXmin = 0; // Min for fitting
	int rangXmax = 0; // Max for fitting
	int nXbins = peak->GetXaxis()->GetNbins(); // Number of bins for fitting

  double binMin = 0.;
  double binMax = 0.;

  for ( int kk=77; kk>27; kk-- ) // from 300 FADC backward
    if ( peakFFTDer->GetBinContent(kk) < 0 )
    {
      if ( binMax < fabs(peakFFTDer->GetBinContent( kk ) ) )
      {
        binMax = fabs(peakFFTDer->GetBinContent(kk));
        rangXmax = peakFFTDer->GetBinCenter(kk);
      }
    }
    else
    {
      binMax = peakFFTDer->GetBinCenter(kk);
      break;
    }

  rangXmax *= 1.5;

  binMin = 0;
  int tmpneg = 0;
  for ( int kk=7; kk<130; kk++ ) // 20 FADC after 0 FADC
  {
    if ( peakFFTDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( peakFFTDer->GetBinContent( kk ) < 0 )
      if ( binMin < fabs( peakFFTDer->GetBinContent( kk ) ) )
      {
        rangXmin = peakFFTDer->GetBinCenter(kk);
        binMin = fabs( peakFFTDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  } 
  //rangXmin += 10;

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( peak->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peak->GetBinCenter(b+1) );
	}

	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  cerr << "MSD " << rangXmin << " " << rangXmax << endl;
  TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(12.7, binMax, -3., 5., -266.2); //Set  init. fit par.
  chFit->Fit("fitFcn","RQ");

  TF1 *expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  TF1 *lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);

  double peakVal = 0.;
  peakVal = peakMaxPk(fitFcn);

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
    if ( peak->GetBinCenter(kbin) >= rangXmin && peak->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( peak->GetBinCenter(kbin) );
      tmp = fitFcn->Eval( peak->GetBinCenter(kbin) ) - peak->GetBinContent(kbin);
      yResid.push_back( tmp / sqrt( peak->GetBinContent(kbin) ) );
      errResid.push_back( sqrt( 
            pow(sqrt( peak->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(fitFcn->Eval( peak->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( peak->GetBinContent(kbin) ) );
    }

  TGraphErrors* residGraph = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *logNormResid = new TH1F("logNormResid", "", 11, -5, 5);
  for ( int val=0; val<yResid.size(); val++ )
    logNormResid->Fill( yResid[val] );
  
  int nPoints = chFit->GetN();
  int nbins2 = ( chFit->GetPointX(nPoints - 1) - chFit->GetPointX(0) ) / 4.;
  TH1F *tmp2 = new TH1F("tmp2", "", nbins2, chFit->GetPointX(0), chFit->GetPointX(nPoints-1));
  
  double x, y;
  
  for(int i=0; i < nPoints; i++ )
  {
    chFit->GetPoint(i, x, y);
    tmp2->SetBinContent(i, y);
  }  

  c4->cd();
  TPaveStats *ptstats;
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  chFit->GetListOfFunctions()->Add(ptstats);
  chFit->SetTitle("");
  chFit->GetXaxis()->SetRangeUser(10, 500);
  chFit->GetYaxis()->SetRangeUser(0, 4000);
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

  leg = new TLegend(0.63,0.34,0.96,0.65);
  leg->AddEntry(peak,"Peak histogram","f");
  leg->AddEntry(fitFcn,"Fit Exp+Lognormal","f");
  leg->AddEntry(expon,"Fit Exp","f");
  leg->AddEntry(lognorm,"Fit Lognormal","f");
  leg->Draw();

  line = new TLine(peakVal, 0, peakVal, 4000);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c4->Print("../plots/peakFittedHisto863.pdf");

  //TRatioPlot *rpmean;
  //c5->cd();
  //TF1 *fitmean21 = chFit->GetFunction("fitFcn");
  //tmp2->Fit(fitmean21,"QR");
  //tmp2->GetXaxis()->SetRangeUser(50, 400);
  //rpmean = new TRatioPlot(tmp2);
  //rpmean->SetGraphDrawOpt("APL*");
  //rpmean->Draw();
  //rpmean->GetLowerRefYaxis()->SetRangeUser(-5, 5);
  //c5->Update();

  c6->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(50, 500);
  residGraph->GetYaxis()->SetRangeUser(-8, 4);
  residGraph->SetLineColor(kBlue);
  residGraph->SetLineWidth(1);
  residGraph->GetXaxis()->SetTitle("[FADC]");
  residGraph->GetYaxis()->SetTitle("Residuals [au]");
  residGraph->SetMarkerStyle(20);
  residGraph->SetMarkerSize(1.5);
  histoStyle(residGraph);
  residGraph->GetYaxis()->SetTitleOffset(0.8);
  residGraph->Draw("APL");

  TString chindf;
  TString valpeak;
  chindf.Form("%.2f", fitFcn->GetChisquare() / fitFcn->GetNDF() );
  valpeak.Form("%.2f", peakVal);

  leg = new TLegend(0.61,0.19,0.96,0.40);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(50, 0, 350, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(peakVal, -8, peakVal, 4);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(50, 1, 350, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(50, -1, 350, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c6->Print("../plots/peakFitResiduals863.pdf");


  // =============================
  // *** *** Fitting Poly2 *** ***

  TCanvas *c7 = canvasStyle("c7");
  TCanvas *c8 = canvasStyle("c8");

  rangXmax = 1.3*binMax;
  rangXmin = 0.7*binMax;

  TGraphErrors* poly2Fit = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2Fit->Fit("poly2","QR");
  peakVal = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0));

  // =====================================
  // *** Computing residuals for Poly2 ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0;

  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( peak->GetBinCenter(kbin) >= rangXmin && peak->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( peak->GetBinCenter(kbin) );
      tmp = poly2->Eval( peak->GetBinCenter(kbin) ) - peak->GetBinContent(kbin); 
      yResid.push_back( tmp / sqrt( peak->GetBinContent(kbin) ) );
      errResid.push_back( sqrt( 
            pow(sqrt( peak->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2->Eval( peak->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( peak->GetBinContent(kbin) ) );
    }

  TGraphErrors* residGraphPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *poly2Resid = new TH1F("poly2Resid", "", 11, -5, 5);
  for ( int val=0; val<yResid.size(); val++ )
    poly2Resid->Fill( yResid[val] );


  c7->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2Fit->GetListOfFunctions()->Add(ptstats);
  poly2Fit->SetTitle("");
  poly2Fit->GetXaxis()->SetRangeUser(10, 500);
  poly2Fit->GetYaxis()->SetRangeUser(0, 4000);
  poly2Fit->SetLineColor(kBlue);
  poly2Fit->SetLineWidth(1);
  poly2Fit->GetXaxis()->SetTitle("[FADC]");
  poly2Fit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2Fit);
  poly2Fit->Draw();

  line = new TLine(peakVal, 0, peakVal, 4000);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c7->Print("../plots/peakFitPoly2863.pdf");

  c8->cd();
  residGraphPoly2->SetTitle("");
  residGraphPoly2->GetXaxis()->SetRangeUser(50, 500);
  residGraphPoly2->GetYaxis()->SetRangeUser(-8, 4);
  residGraphPoly2->SetLineColor(kBlue);
  residGraphPoly2->SetLineWidth(1);
  residGraphPoly2->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2->GetYaxis()->SetTitle("Residuals [au]");
  residGraphPoly2->SetMarkerSize(1.5);
  residGraphPoly2->SetMarkerStyle(20);
  histoStyle(residGraphPoly2);
  residGraphPoly2->GetYaxis()->SetTitleOffset(0.8);
  residGraphPoly2->Draw("APL");

  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
  valpeak.Form("%.2f", peakVal);

  leg = new TLegend(0.14,0.19,0.50,0.40);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","f");
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(rangXmin-5, 0, rangXmax+5, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(peakVal, -8, peakVal, 4);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(rangXmin-10, 1, rangXmax+5, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-10, -1, rangXmax+5, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c8->Print("../plots/peakFitResidualsPoly2863.pdf");

  // ====================================================
  // *** *** *** Aplying for Smooth Histogram *** *** *** 

  // ===============================
  // *** *** Doing for Poly2 *** ***
  TCanvas *c9 = canvasStyle("c9");
  TCanvas *c10 = canvasStyle("c10");

  xbins.clear();
  ycnts.clear();
  yerrs.clear();
	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( peakFFT->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peakFFT->GetBinCenter(b+1) );
	}

  TGraphErrors* poly2FitSmooth = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2Smooth = new TF1("poly2Smooth","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2FitSmooth->Fit("poly2Smooth","QR");

  peakVal = -poly2Smooth->GetParameter(1) / (2.*poly2Smooth->GetParameter(0));

  // =====================================
  // *** Computing residuals for Poly2 ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0;

  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( peakFFT->GetBinCenter(kbin) >= rangXmin && peakFFT->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( peakFFT->GetBinCenter(kbin) );
      tmp = poly2Smooth->Eval( peakFFT->GetBinCenter(kbin) ) - peakFFT->GetBinContent(kbin);
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( peakFFT->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2Smooth->Eval( peakFFT->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraphPoly2Smooth = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *residPoly2Smooth = new TH1F("residPoly2Smooth", "", 33, -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    residPoly2Smooth->Fill( yResid[val] );

  c9->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2FitSmooth->GetListOfFunctions()->Add(ptstats);
  poly2FitSmooth->SetTitle("");
  poly2FitSmooth->GetXaxis()->SetRangeUser(10, 500);
  poly2FitSmooth->GetYaxis()->SetRangeUser(0, 4000);
  poly2FitSmooth->SetLineColor(kBlue);
  poly2FitSmooth->SetLineWidth(1);
  poly2FitSmooth->GetXaxis()->SetTitle("[FADC]");
  poly2FitSmooth->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2FitSmooth);
  poly2FitSmooth->Draw();

  line = new TLine(peakVal, 0, peakVal, 4000);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c9->Print("../plots/peakSmoothFitPoly2863.pdf");

  c10->cd();
  residGraphPoly2Smooth->SetTitle("");
  residGraphPoly2Smooth->GetXaxis()->SetRangeUser(50, 500);
  residGraphPoly2Smooth->GetYaxis()->SetRangeUser(-320,200);
  residGraphPoly2Smooth->SetLineColor(kBlue);
  residGraphPoly2Smooth->SetLineWidth(1);
  residGraphPoly2Smooth->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2Smooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residGraphPoly2Smooth->SetMarkerStyle(20);
  residGraphPoly2Smooth->SetMarkerSize(1.5);
  histoStyle(residGraphPoly2Smooth);
  residGraphPoly2Smooth->Draw("APL");

  chindf.Form("%.2f", poly2Smooth->GetChisquare() / poly2Smooth->GetNDF() );
  valpeak.Form("%.2f", peakVal);


  leg = new TLegend(0.6,0.19,0.96,0.47);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(rangXmin-5, 0, rangXmax+5, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();


  line = new TLine(peakVal, -320, peakVal, 200);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c10->Print("../plots/peakFitSmoothResidualsPoly2863.pdf");


  // ====================================================
  // *** *** Fitting Log-Norm to Smooth histogram *** ***

  TCanvas *c11 = canvasStyle("c11");
  TCanvas *c12 = canvasStyle("c12");

  rangXmax = 1.8*binMax;
  binMax = 0.;

  for ( int kk=77; kk>27; kk-- ) // from 300 FADC backward
    if ( peakFFTDer->GetBinContent(kk) < 0 )
    {
      if ( binMax < fabs(peakFFTDer->GetBinContent( kk ) ) )
      {
        binMax = fabs(peakFFTDer->GetBinContent(kk));
        rangXmax = peakFFTDer->GetBinCenter(kk);
      }
    }
    else
    {
      binMax = peakFFTDer->GetBinCenter(kk);
      break;
    }

  rangXmax *= 1.5;

  binMin = 0;
  tmpneg = 0;
  for ( int kk=7; kk<130; kk++ ) // 20 FADC after 0 FADC
  {
    if ( peakFFTDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( peakFFTDer->GetBinContent( kk ) < 0 )
      if ( binMin < fabs( peakFFTDer->GetBinContent( kk ) ) )
      {
        rangXmin = peakFFTDer->GetBinCenter(kk); // 20 FADC fordward EM-Peak
        binMin = fabs( peakFFTDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  }

  //rangXmin += 10;

  xbins.clear();
  ycnts.clear();
  yerrs.clear();

	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( peakFFT->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peakFFT->GetBinCenter(b+1) );
	}

	chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(12.7, binMax, -3., 5., -266.2); //Set  init. fit par.
  chFit->Fit("fitFcn","QR");
  
  peakVal = peakMaxPk(fitFcn);

  expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);
  
  expon->SetParameter(0, fitFcn->GetParameter(3));
  expon->SetParameter(1, fitFcn->GetParameter(4));
  
  lognorm->SetParameter(0, fitFcn->GetParameter(0));
  lognorm->SetParameter(1, fitFcn->GetParameter(1));
  lognorm->SetParameter(2, fitFcn->GetParameter(2));

  // ===================================================
  // *** *** Doing residuals for Log-Norm Smooth *** ***
  xResid.clear();
  yResid.clear();
  errResid.clear();
  tmp = 0.;
  for ( int kbin=1; kbin<nXbins; kbin++ )
    if ( peakFFT->GetBinCenter(kbin) >= rangXmin && peakFFT->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( peakFFT->GetBinCenter(kbin) );
      tmp = fitFcn->Eval( peakFFT->GetBinCenter(kbin) ) - peakFFT->GetBinContent(kbin);
      yResid.push_back( tmp ); /// sqrt( peak->GetBinContent(kbin) ) );
      errResid.push_back( sqrt( 
            pow(sqrt( peakFFT->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(fitFcn->Eval( peakFFT->GetBinCenter(kbin) ) ) ),2)
            ) );
    }
    
  TGraphErrors *residGraphLogNormSmooth = new TGraphErrors( xResid.size(), &xResid.front(), &yResid.front(), 0, &errResid.front() );

  TH1F *residLogNormSmooth = new TH1F("residLogNormSmooth", "", 33, -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    residLogNormSmooth->Fill( yResid[val] );
      
  c11->cd();

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  chFit->GetListOfFunctions()->Add(ptstats);
  
  chFit->SetTitle("");
  chFit->GetXaxis()->SetRangeUser(10, 500);
  chFit->GetYaxis()->SetRangeUser(0, 4000);
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
                                               
  leg = new TLegend(0.63,0.34,0.96,0.65);
  leg->AddEntry(peak,"Peak histogram","f");
  leg->AddEntry(fitFcn,"Fit Exp+Lognormal","f");
  leg->AddEntry(expon,"Fit Exp","f");
  leg->AddEntry(lognorm,"Fit Lognormal","f");
  leg->Draw();

  line = new TLine(peakVal, 0, peakVal, 4000);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c11->Print("../plots/peakSmoothFittedHisto863.pdf");

  c12->cd();

  residGraphLogNormSmooth->SetTitle("");
  residGraphLogNormSmooth->GetXaxis()->SetRangeUser(50, 500);
  residGraphLogNormSmooth->GetYaxis()->SetRangeUser(-320, 200);
  residGraphLogNormSmooth->SetLineColor(kBlue);
  residGraphLogNormSmooth->SetLineWidth(1);
  residGraphLogNormSmooth->GetXaxis()->SetTitle("[FADC]");
  residGraphLogNormSmooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residGraphLogNormSmooth->SetMarkerStyle(20);
  residGraphLogNormSmooth->SetMarkerSize(1.5);
  histoStyle(residGraphLogNormSmooth);
  residGraphLogNormSmooth->Draw("APL");

  chindf.Form("%.2f", fitFcn->GetChisquare() / fitFcn->GetNDF() );
  valpeak.Form("%.2f", peakVal);

  leg = new TLegend(0.6,0.19,0.96,0.47);
  leg->AddEntry(residGraphLogNormSmooth,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","f"); 
  leg->SetTextSize(0.065);
  leg->Draw();

  line = new TLine(50, 0, 400, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(peakVal, -320, peakVal, 200);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c12->Print("../plots/peakSmoothFitResiduals863.pdf");

  // ========================================================
  // *** *** *** Plotting Residuals distributions *** *** ***

  TCanvas *c13 = canvasStyle("c13");
  TCanvas *c14 = canvasStyle("c14");

  c13->cd();
  logNormResid->SetStats(1);
  logNormResid->SetTitle("");
  logNormResid->SetLineColor(kBlack);
  logNormResid->SetLineWidth(1);
  logNormResid->SetFillColor(kBlack);
  logNormResid->SetFillStyle(3001);
  logNormResid->GetXaxis()->SetTitle("Residuals [au]");
  logNormResid->GetYaxis()->SetTitle("Counts [au]");
  logNormResid->GetYaxis()->SetRangeUser(0, 15.5);
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

  ptstats = new TPaveStats(0.73, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  poly2Resid->SetName("Resid. Poly2");
  poly2Resid->GetListOfFunctions()->Add(ptstats);

  c13->Print("../plots/peakResidualsDist863.pdf");

  c14->cd();
  residLogNormSmooth->SetStats(1);
  residLogNormSmooth->SetTitle("");
  residLogNormSmooth->SetLineColor(kBlack);
  residLogNormSmooth->SetLineWidth(1);
  residLogNormSmooth->SetFillColor(kBlack);
  residLogNormSmooth->SetFillStyle(3001);
  residLogNormSmooth->GetXaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residLogNormSmooth->GetYaxis()->SetTitle("Counts [au]");
  residLogNormSmooth->GetYaxis()->SetRangeUser(0, 11.2);
  residLogNormSmooth->GetXaxis()->SetRangeUser(-90, 90);
  histoStyle(residLogNormSmooth);
  residLogNormSmooth->Draw();

  residPoly2Smooth->SetTitle("");
  residPoly2Smooth->SetLineColor(kBlue);
  residPoly2Smooth->SetFillColor(kBlue);
  residPoly2Smooth->SetFillStyle(3001); 
  residPoly2Smooth->SetLineWidth(1);
  residPoly2Smooth->Draw("sames");

  ptstats = new TPaveStats(0.13, 0.7, 0.36, 0.97,"brNDC");
  ptstats->SetTextColor(kBlack);
  residLogNormSmooth->SetName("Resid. FFT LogNorm");
  residLogNormSmooth->GetListOfFunctions()->Add(ptstats);

  ptstats = new TPaveStats(0.73, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  residPoly2Smooth->SetName("Resid. FFT Poly2");
  residPoly2Smooth->GetListOfFunctions()->Add(ptstats);

  c14->Print("../plots/peakSmoothResidualsDist863.pdf");

}
