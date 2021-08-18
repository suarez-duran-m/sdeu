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
  int binavocero = 0;
  for ( int kk=rawbinMax; kk>0; kk-- ) // 20 FADC after 0 FADC
  {
    if ( tmpBinMin < fabs(h.GetBinContent( kk ) ) )
    {
      tmpBinMin = fabs( h.GetBinContent( kk ) ); 
      tmpMin = h.GetBinCenter(kk-1); // Change of concavity
    }
    if ( h.GetBinContent( kk ) < 0 && nroot==0 && binavocero>3 )
    {
      nroot = 1;
      rawbinMin = h.GetBinCenter(kk); // Bin for local minimum
    }
    if ( h.GetBinContent( kk ) > 0 )
      binavocero++;
    else if ( h.GetBinContent( kk ) > 0 && nroot==1 )
      break;
  }

  minmax.push_back( tmpMin );
  minmax.push_back( tmpMax );
  minmax.push_back( tmpBinMax );
  minmax.push_back( rawbinMin );

  return minmax;
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

void readPeakForFitting()
{
  TString cdasmeName = "uubAoPPMT3St863Mthdec.root"; 
  TFile *cdasmeF = TFile::Open(cdasmeName);
  TTree *cdasmeInfo = (TTree*)cdasmeF->Get("PeakData");

  TH1F *peak = new TH1F();
  TGraphErrors *graphFitted = new TGraphErrors();
  int evt;

  cdasmeInfo->SetBranchAddress("peakForFit", &peak);
  cdasmeInfo->SetBranchAddress("graph", &graphFitted);
  cdasmeInfo->SetBranchAddress("eventId", &evt);

  cdasmeInfo->GetEntry(0);

  int nXbins = peak->GetXaxis()->GetNbins();
  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins; b++ )
    xfadc[b] = peak->GetBinCenter(b);

  TH1F *peakDer = histDerivative(*peak, xfadc); // Peak raw Derivative
  TH1F *peakSmooth = getSmooth(*peak, xfadc); // Smooth Peak
  TH1F *peakSmooDer = histDerivative(*peakSmooth, xfadc); // Peak Smooth Derivative
  TH1F *peakSmooDerSmth = getSmooth(*peakSmooDer, xfadc); // Smooth Peak-Smooth-Derivative

  TLine *line;
  TLegend *leg;
  TPaveStats *ptstats;

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  peak->SetStats(0);
  peak->SetTitle("");
  peak->SetLineColor(kBlue);
  peak->GetXaxis()->SetRangeUser(0, 1000);
  peak->SetLineWidth(1);
  peak->GetXaxis()->SetTitle("[FADC]");
  peak->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(peak);
  peak->Draw();

  peakSmooth->SetLineColor(kOrange+9);
  peakSmooth->SetLineWidth(1);
  peakSmooth->Draw("same");

  TString strEvt;
  strEvt.Form("%d", evt); 

  leg = new TLegend(0.62,0.65,0.95,0.96);
  leg->SetHeader("#splitline{Peak histogrma Station 863}{(Event "+strEvt+")}");
  leg->AddEntry(peak,"Peak histogram","f");
  leg->AddEntry(peakSmooth,"Smooth peak histogram","f");
  leg->Draw();  
  //c1->Print("../plots/peakHisto863.pdf");

  TCanvas *c2 = canvasStyle("c2");
  c2->cd();
  peakDer->SetStats(0);
  peakDer->SetLineColor(kGray);
  peakDer->SetLineWidth(1);
  peakDer->GetXaxis()->SetRangeUser(-10, 1000);
  peakDer->GetXaxis()->SetTitle("[FADC]");
  peakDer->GetYaxis()->SetRangeUser(-15, 10);
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
  line->SetLineWidth(1);
  line->Draw();


  line = new TLine(148, -15, 148, 10);
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
  //c2->Print("../plots/peakDerHisto863.pdf");

  // ===============
  // *** Fitting ***
  double rangXmin = 0.;
  double rangXmax = 0.;
  vector < double > fitRange;

  fitRange = getFitRange( *peakSmooDerSmth );
  rangXmin = fitRange[3];
  rangXmax = fitRange[1];

  vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( int b=0; b<nXbins; b++ ) 
  { 
		ycnts.push_back( peak->GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( peak->GetBinCenter(b+1) );
	}

	TGraphErrors* pkToFit = new TGraphErrors(xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );
  TF1 *poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);

  pkToFit->Fit("poly2","QR");
  double vemPk = -1.*poly2->GetParameter(1) / (2.*poly2->GetParameter(0));
  
  TCanvas *c3 = canvasStyle("c3");
  c3->cd();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  ptstats = new TPaveStats(0.63, 0.67, 0.96, 0.97,"brNDC");
  pkToFit->GetListOfFunctions()->Add(ptstats);
  pkToFit->SetTitle("");
  pkToFit->GetXaxis()->SetRangeUser(10, 500);
  pkToFit->GetYaxis()->SetRangeUser(0, 1e3);
  pkToFit->SetLineColor(kBlue);
  pkToFit->SetLineWidth(1);
  pkToFit->GetXaxis()->SetTitle("[FADC]");
  pkToFit->GetYaxis()->SetTitle("Counts [au]");
  histoStyle(pkToFit);
  pkToFit->Draw();

  line = new TLine(vemPk, 0, vemPk, 1e3);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  //c3->Print("../plots/peakFitPoly2863.pdf");

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;

  getResiduals( peak, poly2, rangXmin, rangXmax, xResid, yResid, errResid );

  TGraphErrors* residPoly2 = new TGraphErrors( xResid.size(), &xResid.front(),
      &yResid.front(), 0, &errResid.front() );


  TCanvas *c4 = canvasStyle("c4");
  c4->cd();
  residPoly2->SetTitle("");
  residPoly2->GetXaxis()->SetRangeUser(50, 500);
  residPoly2->GetYaxis()->SetRangeUser(-8, 4);
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

  leg = new TLegend(0.58,0.19,0.96,0.40);
  leg->AddEntry(residPoly2,"#splitline{ #frac{#chi^{2}}{ndf} = "+chindf+"}{VEM-Peak: "+strVemPk+"}","ep");
  leg->SetTextSize(0.065);
  leg->SetBorderSize(0);
  leg->Draw();

  line = new TLine(rangXmin-5, 0, rangXmax+5, 0);
  line->SetLineStyle(4);
  line->SetLineWidth(2);
  line->Draw();

  line = new TLine(vemPk, -8, vemPk, 4);
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

  //c4->Print("../plots/peakFitResidualsPoly2863.pdf");
}
