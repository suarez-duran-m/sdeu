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

void readFitPeakHisto()
{

  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitPeakHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%speakHist.dat",dir.Data()));

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
  TH1F *peak2Der = histDerivative(*peakDer, xfadc);

  TH1F *peakSmooth = getSmooth(*peak, xfadc);
  TH1F *peakSmooDer = histDerivative(*peakSmooth, xfadc);
  TH1F *peakSmoo2Der = histDerivative(*peakSmooDer, xfadc);

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

  peakSmooth->GetXaxis()->SetRangeUser(-10, 1000);
  peakSmooth->SetLineColor(kOrange+10);
  peakSmooth->SetLineWidth(1);
  peakSmooth->GetXaxis()->SetTitle("[FADC]");
  peakSmooth->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(peakSmooth);
  peakSmooth->Draw("same");

  leg = new TLegend(0.62,0.65,0.95,0.96);
  leg->SetHeader("#splitline{Peak histogrma Station 863}{(Event 61219267)}");
  leg->AddEntry(peak,"Peak histogram","f");
  leg->AddEntry(peakSmooth,"Smooth peak histogram","f");
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
  peakSmooDer->Draw("SAME");

  line = new TLine(0, 0, 1000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.5,0.64,0.95,0.9);
  leg->SetHeader("From Smooth peak histogram","C");
  leg->AddEntry(peakSmooDer,"First derivative","f");
  leg->Draw();
  c3->Print("../plots/peakSmoothDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c4 = canvasStyle("c4");
  TCanvas *c5 = canvasStyle("c5");
  TCanvas *c6 = canvasStyle("c6");

	int emPkb = 0; // Bin for EM peak
	int emPkc = 0; // Counts for EM peak
	int rangXmin = 0; // Min for fitting
	int rangXmax = 0; // Max for fitting
	int nXbins = peak->GetXaxis()->GetNbins(); // Number of bins for fitting

  int binMax = 0;
  int binMin = 0;
  for ( int kk=100; kk>10; kk-- ) // from 600 FADC backward
    if ( peakSmooDer->GetBinContent( kk ) > 0 )
    {
      binMax = peakSmooth->GetBinCenter(kk);
      break;
    }

  rangXmax = 2.*binMax;  

  for ( int kk=7; kk<130; kk++ ) // 20 FADC after 0 FADC
    if ( peakSmooDer->GetBinContent( kk ) < 0 )
    {
      binMin = peakSmooth->GetBinCenter(kk+5); // 20 FADC fordward EM-Peak
      break;
    }

  rangXmin = 1.3*binMin;
  
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

  TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(12.7, 150., -3., 5., -266.2); //Set  init. fit par.
  chFit->Fit("fitFcn","QR");

  TF1 *expon = new TF1("expon", "exp( [0] - [1]/x )", rangXmin, rangXmax);
  TF1 *lognorm = new TF1("lognorm", "(exp([0])/x) * exp( -0.5*pow( (( log([1]) - log(x) )*[2]),2 ) )", rangXmin, rangXmax);

  TGraph* der = (TGraph*)fitFcn->DrawDerivative("goff");

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
      yResid.push_back( tmp ); /// sqrt( peak->GetBinContent(kbin) ) );
      errResid.push_back( sqrt( 
            pow(sqrt( peak->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(fitFcn->Eval( peak->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residGraph = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

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
  c4->Print("../plots/peakFittedHisto863.pdf");

  TRatioPlot *rpmean;
  c5->cd();
  TF1 *fitmean21 = chFit->GetFunction("fitFcn");
  tmp2->Fit(fitmean21,"QR");
  tmp2->GetXaxis()->SetRangeUser(50, 400);
  rpmean = new TRatioPlot(tmp2);
  rpmean->SetGraphDrawOpt("APL*");
  rpmean->Draw();
  rpmean->GetLowerRefYaxis()->SetRangeUser(-5, 5);
  c5->Update();

  c6->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(50, 500);
  residGraph->SetLineColor(kBlue);
  residGraph->SetLineWidth(1);
  residGraph->GetXaxis()->SetTitle("[FADC]");
  residGraph->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residGraph);
  residGraph->Draw("APL*");

  line = new TLine(50, 0, 400, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c6->Print("../plots/peakFitResiduals863.pdf");


  // =============================
  // *** *** Fitting Poly2 *** ***

  TCanvas *c7 = canvasStyle("c7");
  TCanvas *c8 = canvasStyle("c8");
  TCanvas *c9 = canvasStyle("c9");
  TCanvas *c10 = canvasStyle("c10");
  TCanvas *c11 = canvasStyle("c11");

  rangXmax = 1.3*binMax;
  rangXmin = 0.8*binMax;

  TGraphErrors* poly2Fit = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  poly2Fit->Fit("poly2","QR");

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
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( peak->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2->Eval( peak->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TGraphErrors* residPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *hResid = new TH1F("hResid", "", 601/20., -300, 300);
  for ( int val=0; val<yResid.size(); val++ )
    hResid->Fill( yResid[val] );


  c7->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2Fit->GetListOfFunctions()->Add(ptstats);
  poly2Fit->SetTitle("");
  poly2Fit->GetXaxis()->SetRangeUser(50, 500);
  poly2Fit->SetLineColor(kBlue);
  poly2Fit->SetLineWidth(1);
  poly2Fit->GetXaxis()->SetTitle("[FADC]");
  poly2Fit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2Fit);
  poly2Fit->Draw();
  c7->Print("../plots/peakFitPoly2863.pdf");

  c8->cd();
  residPoly2->SetTitle("");
  residPoly2->GetXaxis()->SetRangeUser(50, 500);
  residPoly2->SetLineColor(kBlue);
  residPoly2->SetLineWidth(1);
  residPoly2->GetXaxis()->SetTitle("[FADC]");
  residPoly2->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  histoStyle(residPoly2);
  residPoly2->Draw("APL*");

  line = new TLine(115, 0, 220, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c8->Print("../plots/peakFitResidualsPoly2863.pdf");

  // ====================================
  // *** Aplying for Smooth Histogram ***

  xbins.clear();
  ycnts.clear();
  yerrs.clear();
	for( int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( peakSmooth->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peakSmooth->GetBinCenter(b+1) );
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
    if ( peakSmooth->GetBinCenter(kbin) >= rangXmin && peakSmooth->GetBinCenter(kbin) <= rangXmax )
    {
      xResid.push_back( peakSmooth->GetBinCenter(kbin) );
      tmp = poly2->Eval( peakSmooth->GetBinCenter(kbin) ) - peakSmooth->GetBinContent(kbin); 
      yResid.push_back( tmp );
      errResid.push_back( sqrt( 
            pow(sqrt( peakSmooth->GetBinContent(kbin) ),2) 
            + pow(sqrt( sqrt(poly2Smooth->Eval( peakSmooth->GetBinCenter(kbin) ) ) ),2)
            ) );
    }

  TH1F *hResidSmooth = new TH1F("hResidSmooth", "", 201/10., -100, 100);
  for ( int val=0; val<yResid.size(); val++ )
    hResidSmooth->Fill( yResid[val] );

  TGraphErrors* residPoly2Smooth = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  c9->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  poly2FitSmooth->GetListOfFunctions()->Add(ptstats);
  poly2FitSmooth->SetTitle("");
  poly2FitSmooth->GetXaxis()->SetRangeUser(50, 500);
  poly2FitSmooth->SetLineColor(kBlue);
  poly2FitSmooth->SetLineWidth(1);
  poly2FitSmooth->GetXaxis()->SetTitle("[FADC]");
  poly2FitSmooth->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(poly2FitSmooth);
  poly2FitSmooth->Draw();
  c9->Print("../plots/peakSmoothFitPoly2863.pdf");

  c10->cd();
  residPoly2Smooth->SetTitle("");
  residPoly2Smooth->GetXaxis()->SetRangeUser(50, 500);
  residPoly2Smooth->SetLineColor(kBlue);
  residPoly2Smooth->SetLineWidth(1);
  residPoly2Smooth->GetXaxis()->SetTitle("[FADC]");
  residPoly2Smooth->GetYaxis()->SetTitle("y_{fit} - y_{data} [au]");
  residPoly2Smooth->GetYaxis()->SetRangeUser(-330,230);
  histoStyle(residPoly2Smooth);
  residPoly2Smooth->Draw("APL*");

  line = new TLine(115, 0, 220, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  c10->Print("../plots/peakFitSmoothResidualsPoly2863.pdf");

  c11->cd();
  hResid->SetStats(1);
  hResid->SetTitle("");
  hResid->SetLineColor(kBlack);
  hResid->SetLineWidth(1);
  hResid->SetFillColor(kBlack);
  hResid->SetFillStyle(3001);
  hResid->GetXaxis()->SetTitle("y_{fit} - y_{data} [au]");
  hResid->GetYaxis()->SetTitle("Counts [au]");
  hResid->GetYaxis()->SetRangeUser(0, 6.5);
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

  ptstats = new TPaveStats(0.73, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kBlue);
  hResidSmooth->SetName("Resid. Smooth Peak histogram");
  hResidSmooth->GetListOfFunctions()->Add(ptstats);

  c11->Print("../plots/peakResidualsDist863.pdf");
}
