#include <iostream>
#include <time.h>

#include <TF1.h>

#include "fitcharge.h"


// Initializing class
fitcharge::fitcharge(TH1F *inHist, int corr, int binsLeftRigh) : rangXmin(0), 
  rangXmax(0), ndf(0.), prob(0.), par0(0.), par1(0.), par2(0.), chi2(0.),
  qpkPos(0.), qpkPosDeri(0.) {
    
  hstCrr = new TH1F();
  histSmooth = new TH1F();
  histSmooDer = new TH1F();
  histSmooDerSmth = new TH1F();

  SetCrr(inHist, corr);
  GetFit(binsLeftRigh); 
}

void fitcharge::SetCrr(TH1F *inHist, int corr) {
	unsigned int nb = 601;
	double xb [nb];

	for ( unsigned int b=1; b<nb; b++ )
		xb[b] = inHist->GetBinCenter(b) - 4 - inHist->GetBinCenter(0); //corr;
  xb[0] = 0;
	xb[nb] = xb[nb-1] + inHist->GetBinWidth(1);
	hstCrr = new TH1F("hstCrr", "", nb, xb);

	for ( unsigned b=0; b<nb; b++ ) {
    if ( b < 401 )
      hstCrr->SetBinContent( b, inHist->GetBinContent(b)/inHist->GetBinWidth(b) );
    else 
      hstCrr->SetBinContent( b, (4.*inHist->GetBinContent(b))/inHist->GetBinWidth(b) );
  }
}


void fitcharge::GetFit(int binsLeftRigh) {
  TString tmpname;
  tmpname.Form("%d",rand());

	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	unsigned nXbins = hstCrr->GetXaxis()->GetNbins(); // Number of bins for fitting
  vector < double > rangeValues;

  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins; b++ )
    xfadc[b] = hstCrr->GetBinCenter(b);

  histSmooth = SetSmooth(hstCrr, xfadc);
  histSmooDer = SetHistDerivative(xfadc);
  // histSmooDerSmth Smoothing twice
  histSmooDerSmth = SetSmooth(histSmooDer, xfadc);
  histSmooDerSmth = SetSmooth(histSmooDerSmth, xfadc);

  // Setting the range for fitting using the histo derivative
  rangeValues = SetRangeFit(*histSmooDerSmth, binsLeftRigh);
  rangXmin = rangeValues[0];
  rangXmax = rangeValues[1];

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( unsigned int b=0; b<nXbins; b++ ) {
		ycnts.push_back( hstCrr->GetBinContent(b) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( hstCrr->GetBinCenter(b) );
  }
	
	TGraphErrors* histToFit = new TGraphErrors();
  histToFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2;

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;
  
  poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  histToFit->Fit("poly2","QR");
  
  chi2 = poly2->GetChisquare();
  ndf = poly2->GetNDF();
  prob = poly2->GetProb();
	qpkPos = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0));
  qpkPosDeri = rangeValues[2];	
  par0 = poly2->GetParameter(0);
  par1 = poly2->GetParameter(1);
  par2 = poly2->GetParameter(2);
  fittedGraph = (TGraphErrors*)histToFit->Clone();  

	delete histToFit;
	delete poly2;
}

TH1F *fitcharge::SetHistDerivative(double xb[]) {// Central differences
  int nbins = 600;
  int h = 0;
  //TString tmpname;
  //tmpname.Form("%d",rand());

  TH1F *derihist = new TH1F();
  derihist = new TH1F("histSmooDer", "", nbins, xb);
  double der = 0.;
  for ( int kk=1; kk<nbins-1; kk++ )
  {
    h = ( hstCrr->GetBinCenter( kk+1 ) - hstCrr->GetBinCenter( kk ) );
    der = ( hstCrr->GetBinContent( kk+1 ) -  hstCrr->GetBinContent( kk-1 ) )/ (2.*h);
    derihist->SetBinContent( kk, der );

  }
  return derihist;
}


TH1F *fitcharge::SetSmooth(TH1F *hist, double xb[]) {
	unsigned int nb = 600;
  double yi = 0.;

  TString tmpname = Form("%d",rand());

	TH1F *hstSmooth = new TH1F();
  hstSmooth = new TH1F(tmpname, "", nb, xb);

	for ( unsigned b=0; b<nb; b++ ) {
    if ( b>6 && b<nb-7 ) {
      yi = hist->GetBinContent(b+1 - 7) 
        + hist->GetBinContent(b+1 - 6) 
        + hist->GetBinContent(b+1 - 5) 
        + hist->GetBinContent(b+1 - 4) 
        + hist->GetBinContent(b+1 - 3) 
        + hist->GetBinContent(b+1 - 2)
        + hist->GetBinContent(b+1 - 1) 
        + hist->GetBinContent(b+1) 
        + hist->GetBinContent(b+1 + 1) 
        + hist->GetBinContent(b+1 + 2)
        + hist->GetBinContent(b+1 + 3)
        + hist->GetBinContent(b+1 + 4)
        + hist->GetBinContent(b+1 + 5)
        + hist->GetBinContent(b+1 + 6)
        + hist->GetBinContent(b+1 + 7);
      yi = yi/15.;
    }
    else
      yi = hist->GetBinContent(b+1);
    
    hstSmooth->SetBinContent(b+1, yi);
  }
  return hstSmooth;
}


vector<double> fitcharge::SetRangeFit( TH1F &hDer, int rightleftBins ) {
  vector < double > minmax;
  int minRng = 0; // Min for fitting
	int maxRng = 0; // Max for fitting
  double binMax = 0.;
  double rawbinMax = 0.;

  for ( int kk=275; kk>100; kk-- ) { // from ~2000 FADC backward
    if ( hDer.GetBinContent(kk) > 0 ) {
      if ( hDer.GetBinCenter(kk) < fabs(hDer.GetBinCenter(kk-1)) ) {
        binMax = hDer.GetBinCenter(kk); // tmp FADC for VEM
        rawbinMax = kk; // Bin for tmp VEM
      }
      else {
        binMax = hDer.GetBinCenter(kk-1);
        rawbinMax = kk; // Bin for tmp VEM
      }
      break;
    }
  } 
  minRng = hDer.GetBinCenter(rawbinMax-rightleftBins);
  maxRng = hDer.GetBinCenter(rawbinMax+rightleftBins);
  minmax.push_back( minRng );
  minmax.push_back( maxRng );
  minmax.push_back( binMax ); 

  return minmax;
}

TGraphErrors* fitcharge::GetFitGraph() {
	return fittedGraph;
}


TH1F *fitcharge::GetChCrr() {
	return hstCrr;
}

void fitcharge::DeleteObjts() {
  fittedGraph->Delete();
  hstCrr->Delete();
  histSmooth->Delete();
  histSmooDer->Delete();
  histSmooDerSmth->Delete();
}
