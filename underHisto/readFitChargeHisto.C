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

// ==================================
// *** ***  *** MAIN CODE *** *** *** 

void readFitChargeHisto()
{

  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitChargeHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
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
  TH1F *charge2Der = histDerivative(*chargeDer, xfadc);

  TH1F *chargeSmooth = getSmooth(*charge, xfadc);
  TH1F *chargeSmooDer = histDerivative(*chargeSmooth, xfadc);
  TH1F *chargeSmoo2Der = histDerivative(*chargeSmooDer, xfadc);

  //TH1F *test = new TH1F(*charge);
  //test->Smooth(50);
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
  //Let's look at the output
  hbC = TH1::TransformHisto(fft_back,hbC,"Re");
  hbC->SetTitle("The backward transform result");

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

  TLatex *latex; 

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
  leg->AddEntry(chargeSmooDer,"First derivative","f");
  leg->AddEntry(chargeFFTDer,"First derivative FFT","f");
  leg->Draw();
  c4->Print("../plots/chargeSmoothDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c5 = canvasStyle("c5");
  TCanvas *c6 = canvasStyle("c6");
  TCanvas *c7 = canvasStyle("c7");

	int emPkb = 0; // Bin for EM peak
	int emPkc = 0; // Counts for EM peak
	int rangXmin = 0; // Min for fitting
	int rangXmax = 0; // Max for fitting
	int nXbins = charge->GetXaxis()->GetNbins(); // Number of bins for fitting

  int binMax = 0;
  int binMin = 0;
  for ( int kk=370; kk>10; kk-- ) // from 2960 FADC backward
    if ( chargeFFTDer->GetBinContent( kk ) > 0 )
    {
      binMax = charge->GetBinCenter(kk);
      break;
    }

  rangXmax = binMax + 0.5*(3000-binMax);
  
  for ( int kk=25; kk<binMax; kk++ ) // 200 FADC after 0 FADC
    if ( chargeFFTDer->GetBinContent( kk ) < 0 )
    {
      binMin = charge->GetBinCenter(kk);
      break;
    }
  
  rangXmin = binMin*(1.5);
  
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

  TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76);
  chFit->Fit("fitFcn","QR");

  TF1 *expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  TF1 *lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);

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

  leg = new TLegend(0.2,0.67,0.56,0.97);
  leg->AddEntry(charge,"Charge histogram","f");
  leg->AddEntry(fitFcn,"Fit Exp+Lognormal","f");
  leg->AddEntry(expon,"Fit Exp","f");
  leg->AddEntry(lognorm,"Fit Lognormal","f");
  leg->Draw();
  
  c5->Print("../plots/chargeFittedHisto863.pdf");

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

  c7->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(100, 3000);
  residGraph->GetYaxis()->SetRangeUser(-120, 160);
  residGraph->SetLineColor(kBlue);
  residGraph->SetLineWidth(1);
  residGraph->GetXaxis()->SetTitle("[FADC * 8.33 ns]");
  residGraph->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residGraph);
  residGraph->Draw("APL*");

  line = new TLine(400, 0, 2200, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c7->Print("../plots/chargeFitResiduals863.pdf");


  // =============================
  // *** *** Fitting Poly2 *** ***

  TCanvas *c8 = canvasStyle("c8");
  TCanvas *c9 = canvasStyle("c9");
  TCanvas *c10 = canvasStyle("c10");
  TCanvas *c11 = canvasStyle("c11");
  TCanvas *c12 = canvasStyle("c12");

  TGraphErrors* poly2Fit = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  rangXmin = binMax*0.8;
  rangXmax = binMax*1.3;
  cerr << binMax << endl;
  cerr << "MSD " << rangXmin << " " << rangXmax << endl;

  TF1 *poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2Fit->Fit("poly2","QR");

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

  TGraphErrors* residPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *hResid = new TH1F("hResid", "", 601/20., -300, 300);
  for ( int val=0; val<yResid.size(); val++ )
    hResid->Fill( yResid[val] );


  c8->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2Fit->GetListOfFunctions()->Add(ptstats);
  poly2Fit->SetTitle("");
  poly2Fit->GetXaxis()->SetRangeUser(100, 3000);
  poly2Fit->SetLineColor(kBlue);
  poly2Fit->SetLineWidth(1);
  poly2Fit->GetXaxis()->SetTitle("[FADC]");
  poly2Fit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2Fit);
  poly2Fit->Draw();

  c8->Print("../plots/chargeFitPoly2863.pdf");

  c9->cd();
  residPoly2->SetTitle("");
  residPoly2->GetYaxis()->SetRangeUser(-120, 160);
  residPoly2->SetLineColor(kBlue);
  residPoly2->SetLineWidth(1);
  residPoly2->GetXaxis()->SetTitle("[FADC]");
  residPoly2->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residPoly2);
  residPoly2->Draw("APL*");

  line = new TLine(rangXmin-50, 0, rangXmax+50, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c9->Print("../plots/chargeFitResidualsPoly2863.pdf");

  // ====================================
  // *** Aplying for Smooth Histogram ***

  xbins.clear();
  ycnts.clear();
  yerrs.clear();
	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( chargeFFT->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( chargeFFT->GetBinCenter(b+1) );
	}

  TGraphErrors* poly2FitSmooth = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2Smooth = new TF1("poly2Smooth","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2FitSmooth->Fit("poly2Smooth","QR");

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
      tmp = poly2->Eval( chargeFFT->GetBinCenter(kbin) ) - chargeFFT->GetBinContent(kbin); 
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( chargeFFT->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2Smooth->Eval( chargeFFT->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TH1F *hResidSmooth = new TH1F("hResidSmooth", "", 201/10., -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    hResidSmooth->Fill( yResid[val] );

  TGraphErrors* residPoly2Smooth = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  c10->cd();
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
  histoStyle(poly2FitSmooth);
  poly2FitSmooth->Draw();
  c10->Print("../plots/chargeFFTFitPoly2863.pdf");

  c11->cd();
  residPoly2Smooth->SetTitle("");
  residPoly2Smooth->SetLineColor(kBlue);
  residPoly2Smooth->SetLineWidth(1);
  residPoly2Smooth->GetXaxis()->SetTitle("[FADC]");
  residPoly2Smooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residPoly2Smooth->GetYaxis()->SetRangeUser(-120,160);
  histoStyle(residPoly2Smooth);
  residPoly2Smooth->Draw("APL*");

  line = new TLine(rangXmin-50, 0, rangXmax+50, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c11->Print("../plots/chargeFitFFTResidualsPoly2863.pdf");

  c12->cd();
  hResid->SetStats(1);
  hResid->SetTitle("");
  hResid->SetLineColor(kBlack);
  hResid->SetLineWidth(1);
  hResid->SetFillColor(kBlack);
  hResid->SetFillStyle(3001);
  hResid->GetXaxis()->SetTitle("y_{fit} - y_{data} [au]");
  hResid->GetYaxis()->SetTitle("Counts [au]");
  hResid->GetYaxis()->SetRangeUser(0, 38);
  hResid->GetXaxis()->SetRangeUser(-100, 100);
  histoStyle(hResid);
  hResid->Draw();

  hResidSmooth->SetTitle("");
  hResidSmooth->SetLineColor(kBlue);
  hResidSmooth->SetFillColor(kBlue);
  hResidSmooth->SetFillStyle(3001); 
  hResidSmooth->SetLineWidth(1);
  hResidSmooth->Draw("sames");

  ptstats = new TPaveStats(0.13, 0.7, 0.36, 0.97,"brNDC");
  ptstats->SetTextColor(kBlack);
  hResid->SetName("Resid. Peak histogram");
  hResid->GetListOfFunctions()->Add(ptstats);

  ptstats = new TPaveStats(0.7, 0.67, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hResidSmooth->SetName("Resid. FFT Peak histogram");
  hResidSmooth->GetListOfFunctions()->Add(ptstats);

  c12->Print("../plots/chargeResidualsDist863.pdf");
}
