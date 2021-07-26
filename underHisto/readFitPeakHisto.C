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

  TString tmpname;
  tmpname.Form("%d",rand());

	TH1F *hstSmooth = new TH1F(tmpname, "", nb, xb);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b > 2 && b<147 )
    {
      yi = hist.GetBinContent(b+1 - 3) 
        + hist.GetBinContent(b+1 - 2)
        + hist.GetBinContent(b+1 - 1) 
        + hist.GetBinContent(b+1) 
        + hist.GetBinContent(b+1 + 1) 
        + hist.GetBinContent(b+1 + 2)
        + hist.GetBinContent(b+1 + 3);
      yi = yi/7.;
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


vector < double > getFitRange( TH1F &h )
{
  vector < double > minmax;
  int tmpMin = 0; // Min for fitting
	int tmpMax = 0; // Max for fitting

  double tmpBinMin = 0.;
  double tmpBinMax = 0.;
  double rawbinMax = 0.;
  double rawbinMin = 0.;

  for ( int kk=77; kk>27; kk-- ) // from 300 FADC backward
    if ( h.GetBinContent(kk) < 0 )
    {
      if ( tmpBinMax < fabs( h.GetBinContent( kk ) ) )
      {
        tmpBinMax = fabs( h.GetBinContent(kk));
        tmpMax = h.GetBinCenter(kk); // Change of concavity
      }
    }
    else
    {
      tmpBinMax = h.GetBinCenter(kk); // tmp FADC for VEM
      rawbinMax = kk; // Bin for tmp VEM
      break;
    }

  tmpBinMin = 0;
  int nroot = 0;
  for ( int kk=rawbinMax; kk>0; kk-- ) // 20 FADC after 0 FADC
  {
    if ( tmpBinMin < fabs(h.GetBinContent( kk ) ) )
    {
      tmpBinMin = fabs( h.GetBinContent( kk ) ); 
      tmpMin = h.GetBinCenter(kk-1); // Change of concavity
    }
    if ( h.GetBinContent( kk ) < 0 && nroot==0 )
    {
      nroot = 1;
      rawbinMin = h.GetBinCenter(kk); // Bin for local minimum
    }
    else if ( h.GetBinContent( kk ) > 0 && nroot==1 )
      break;
  }

  minmax.push_back( tmpMin );
  minmax.push_back( tmpMax );
  minmax.push_back( tmpBinMax );
  minmax.push_back( rawbinMin );

  return minmax;
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
      err.push_back( 
          sqrt( pow(sqrt( hist->GetBinContent(kbin) ),2)
            + pow(sqrt( sqrt(func->Eval( hist->GetBinCenter(kbin) ) ) ),2)
            ) / sqrt( hist->GetBinContent(kbin) )
          );
    }
}


// ==================================
// *** ***  *** MAIN CODE *** *** *** 

void readFitPeakHisto()
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("readFitPeakHisto.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%ssample61219267.dat",dir.Data()));
  //in.open(Form("%speakHist61219267.dat",dir.Data()));
  
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

  TH1F *peak = new TH1F("peak", "", nbins, xfadc);
  
  for ( int kk=0; kk<nbins; kk++ )
    peak->SetBinContent(kk, ycnt[kk]);

  TH1F *peakDer = histDerivative(*peak, xfadc); // Peak raw Derivative
  TH1F *peakSmooth = getSmooth(*peak, xfadc); // Smooth Peak
  TH1F *peakSmooDer = histDerivative(*peakSmooth, xfadc); // Peak Smooth Derivative
  TH1F *peakSmooDerSmth = getSmooth(*peakSmooDer, xfadc); // Smooth Peak-Smooth-Derivative

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

  peakSmooth->SetLineColor(kOrange+9);
  peakSmooth->SetLineWidth(1);
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

  peakSmooDer->SetStats(0);
  peakSmooDer->SetLineColor(kOrange+9);
  peakSmooDer->SetLineWidth(1);
  peakSmooDer->Draw("sames");

  peakSmooDerSmth->SetStats(0);
  peakSmooDerSmth->SetLineColor(kGreen+3);
  peakSmooDerSmth->SetLineWidth(1);
  peakSmooDerSmth->Draw("sames");

  line = new TLine(0, 0, 1000, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();
  
  leg = new TLegend(0.45,0.7,0.95,0.97);
  leg->SetHeader("From peak histogram","C");
  leg->AddEntry(peakDer,"First derivative Peak-H.","f");
  leg->AddEntry(peakSmooDer,"First derivative Smooth-Peak-H.","f");
  leg->AddEntry(peakSmooDerSmth,"Smooth First derivative Smooth-Peak-H.","f");
  leg->SetTextSize(0.04);
  leg->Draw();
  c2->Print("../plots/peakDerHisto863.pdf");


  // ===============================
  // *** *** *** FITTING *** *** *** 

  TCanvas *c3 = canvasStyle("c3");
  TCanvas *c4 = canvasStyle("c4");

  vector < double > fitRange;
  fitRange = getFitRange( *peakSmooDerSmth );

  vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

  int nXbins = peak->GetXaxis()->GetNbins(); // Number of bins for fitting

	for( int b=0; b<nXbins; b++ ) 
  { 
		ycnts.push_back( peak->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peak->GetBinCenter(b+1) );
	}

  double rangXmin = 0.;
  double rangXmax = 0.;
  double binMax = 0.;

  rangXmin = 1.2*fitRange[0]; // For Exp + LogNorm
  rangXmax = 1.4*fitRange[1];
  binMax = fitRange[2];

	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

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

  getResiduals( peak, fitFcn, rangXmin, rangXmax, xResid, yResid, errResid );

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

  c3->cd();
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

  c3->Print("../plots/peakFittedHisto863.pdf");

  c4->cd();
  residGraph->SetTitle("");
  residGraph->GetXaxis()->SetRangeUser(50, 500);
  residGraph->GetYaxis()->SetRangeUser(-8, 4);
  residGraph->SetLineColor(kBlack);
  residGraph->SetLineWidth(1);
  residGraph->GetXaxis()->SetTitle("[FADC]");
  residGraph->GetYaxis()->SetTitle("Residuals");
  residGraph->SetMarkerStyle(20);
  residGraph->SetMarkerColor(kBlack);
  residGraph->SetMarkerSize(2);
  histoStyle(residGraph);
  residGraph->GetYaxis()->SetTitleOffset(0.8);
  residGraph->Draw("APL");

  TString chindf;
  TString valpeak;
  chindf.Form("%.2f", fitFcn->GetChisquare() / fitFcn->GetNDF() );
  valpeak.Form("%.2f", peakVal);

  leg = new TLegend(0.61,0.19,0.96,0.40);
  leg->AddEntry(residGraph,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","ep"); 
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

  line = new TLine(rangXmin-5, 1, rangXmax+5, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-5, -1, rangXmax+5, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c4->Print("../plots/peakFitResiduals863.pdf");


  // =============================
  // *** *** Fitting Poly2 *** ***

  TCanvas *c5 = canvasStyle("c5");
  TCanvas *c6 = canvasStyle("c6");

  rangXmin = fitRange[3]; // For Poly2
  rangXmax = fitRange[1];

  cerr << rangXmin << " " << rangXmax << endl;

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

  getResiduals( peak, poly2, rangXmin, rangXmax, xResid, yResid, errResid );

  TGraphErrors* residGraphPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );

  TH1F *poly2Resid = new TH1F("poly2Resid", "", 14, -7, 7);
  for ( int val=0; val<yResid.size(); val++ )
    poly2Resid->Fill( yResid[val] );


  c5->cd();
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

  c5->Print("../plots/peakFitPoly2863.pdf");

  c6->cd();
  residGraphPoly2->SetTitle("");
  residGraphPoly2->GetXaxis()->SetRangeUser(50, 500);
  residGraphPoly2->GetYaxis()->SetRangeUser(-8, 4);
  residGraphPoly2->SetLineColor(kBlack);
  residGraphPoly2->SetLineWidth(1);
  residGraphPoly2->GetXaxis()->SetTitle("[FADC]");
  residGraphPoly2->GetYaxis()->SetTitle("Residuals");
  residGraphPoly2->SetMarkerSize(3);
  residGraphPoly2->SetMarkerColor(kBlack);
  residGraphPoly2->SetMarkerStyle(23);
  histoStyle(residGraphPoly2);
  residGraphPoly2->GetYaxis()->SetTitleOffset(0.8);
  residGraphPoly2->Draw("APL");

  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF() );
  valpeak.Form("%.2f", peakVal);

  leg = new TLegend(0.61,0.19,0.96,0.40);
  leg->AddEntry(residGraphPoly2,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{Peak val.: "+valpeak+"}","ep");
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

  line = new TLine(rangXmin-5, 1, rangXmax+5, 1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  line = new TLine(rangXmin-5, -1, rangXmax+5, -1);
  line->SetLineStyle(9);
  line->SetLineWidth(1);
  line->Draw();

  c6->Print("../plots/peakFitResidualsPoly2863.pdf");

  // ==============================
  // *** *** Comparing Chi2 *** ***
  TCanvas *c7 = canvasStyle("c7");
  TCanvas *c8 = canvasStyle("c8");

  // ===========
  // Using the same range for Poly2
  rangXmin = fitRange[3];
  rangXmax = fitRange[1];

  xResid.clear();
  yResid.clear();
  errResid.clear();

  getResiduals( peak, fitFcn, rangXmin, rangXmax, xResid, yResid, errResid );

  TH1F *logNormResid2 = new TH1F("logNormResid2", "", 14, -7, 7);
  for ( int val=0; val<yResid.size(); val++ )
    logNormResid2->Fill( yResid[val] );

  double sum = 0.;
  for ( int kk=0; kk<xResid.size(); kk++ )
    sum += yResid[kk]*yResid[kk];

  TGraphErrors *residGraph2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );
  peakVal = peakMaxPk(fitFcn);

  c7->cd();
  residGraph2->SetTitle("");
  residGraph2->GetYaxis()->SetRangeUser(-13, 5);
  residGraph2->SetLineColor(kRed);
  residGraph2->SetMarkerColor(kRed);
  residGraph2->SetLineWidth(2);
  residGraph2->GetXaxis()->SetTitle("[FADC]");
  residGraph2->GetYaxis()->SetTitle("Residuals");
  residGraph2->SetMarkerStyle(20);
  residGraph2->SetMarkerSize(2);
  histoStyle(residGraph2);
  residGraph2->GetYaxis()->SetTitleOffset(0.8);
  residGraph2->Draw("AP same");

  residGraphPoly2->SetLineColor(kGreen+3);
  residGraphPoly2->SetLineWidth(2);
  residGraphPoly2->SetMarkerStyle(23);
  residGraphPoly2->SetMarkerColor(kGreen+3);
  residGraphPoly2->SetMarkerSize(3);
  residGraphPoly2->Draw("p same");
  
  leg = new TLegend(0.12,0.16,0.75,0.48);
  TString chi2;
  chi2.Form("%.2f", sum );
  TString ndf;
  ndf.Form("%lu",  xResid.size()-5);
  valpeak.Form("%.2f", peakVal);
  leg->AddEntry(residGraph2,"#splitline{Exp+Lognormal: #chi^{2}/ndf = "+chi2+"/"+ndf+" = 3.52}{Peak val.: "+valpeak+"}","ep"); 
  chindf.Form("%.2f", poly2->GetChisquare() / poly2->GetNDF()  );
  peakVal = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0)); 
  valpeak.Form("%.2f", peakVal);
  chi2.Form("%.2f", poly2->GetChisquare() );
  ndf.Form("%d", poly2->GetNDF());
  leg->AddEntry(residGraphPoly2,"#splitline{Poly2: #chi^{2}/ndf = "+chi2+"/"+ndf+" = 3.09}{Peak val.: "+valpeak+"}","ep"); 
  leg->SetTextSize(0.05);
  leg->Draw();
  
  line = new TLine(rangXmin-5, 0, rangXmax+5, 0);
  line->SetLineStyle(4);
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

  c7->Print("../plots/peakResidualsCompar863.pdf");

  xbins.clear();
  ycnts.clear();
  yerrs.clear();

  for( int b=0; b<nXbins; b++ )
  {
		ycnts.push_back( peak->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peak->GetBinCenter(b+1) );
	}

  TGraphErrors* chFit2 = new TGraphErrors( xbins.size(), &xbins.front(),
      &ycnts.front(), 0, &yerrs.front() );

  c8->cd();
  gStyle->SetOptStat(0);
  chFit2->SetTitle("");
  chFit2->GetXaxis()->SetRangeUser(30, 350);
  chFit2->GetYaxis()->SetRangeUser(1000, 3700);
  chFit2->SetLineColor(kBlue);
  chFit2->SetLineWidth(1);
  chFit2->GetXaxis()->SetTitle("[FADC]");
  chFit2->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(chFit2);
  chFit2->Draw();

  poly2->SetLineColor(kGreen+3);
  poly2->Draw("sames");

  fitFcn->SetLineColor(kRed);
  fitFcn->Draw("sames");

  leg = new TLegend(0.71,0.75,0.96,0.96);
  leg->AddEntry(residGraph2,"Exp+Lognormal","l");
  leg->AddEntry(residGraphPoly2,"Poly2","l"); 
  leg->SetTextSize(0.05);
  leg->Draw();

  c8->Print("../plots/peakFitCompar863.pdf");

  // ======================================
  // *** *** Residuals distribution *** ***

  TCanvas *c9 = canvasStyle("c9");
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);

  c9->cd();
  logNormResid2->SetStats(1);
  logNormResid2->SetTitle("");
  logNormResid2->SetLineColor(kRed);
  logNormResid2->SetFillColor(kRed);
  logNormResid2->SetFillStyle(3354);
  logNormResid2->SetLineWidth(1);
  logNormResid2->GetXaxis()->SetTitle("Residual");
  logNormResid2->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(logNormResid2);
  logNormResid2->Draw();

  poly2Resid->SetTitle("");
  poly2Resid->SetLineColor(kGreen+3);
  poly2Resid->SetFillColor(kGreen+3);
  poly2Resid->SetFillStyle(3345);
  poly2Resid->SetLineWidth(1);
  poly2Resid->Draw("sames");

  ptstats = new TPaveStats(0.13, 0.7, 0.36, 0.97,"brNDC");
  ptstats->SetTextColor(kRed);
  logNormResid2->SetName("Residuals Exp.+ LogNorm");
  logNormResid2->GetListOfFunctions()->Add(ptstats);

  ptstats = new TPaveStats(0.72, 0.7, 0.96, 0.97,"brNDC");
  ptstats->SetTextColor(kGreen+3);
  poly2Resid->SetName("Residual Poly2");
  poly2Resid->GetListOfFunctions()->Add(ptstats);

  c9->Print("../plots/peakResidualsDist863.pdf");
}
