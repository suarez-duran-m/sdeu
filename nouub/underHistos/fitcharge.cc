#include <iostream>
#include <time.h>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TVirtualFFT.h"

#include "fitcharge.h"

// =======================
// *** Local Functions ***
// =======================

TH1F *histDerivativeCh(TH1F &hist, double xb[]) { // Central differences
  int nbins = 600;
  int h = 0;
  TString tmpname;
  tmpname.Form("%d",rand());

  TH1F *derihist = new TH1F(tmpname, "", nbins, xb);
  double der = 0.;
  for ( int kk=1; kk<nbins-1; kk++ ) {
    h = ( hist.GetBinCenter( kk+1 ) - hist.GetBinCenter( kk ) );
    der = ( hist.GetBinContent( kk+1 ) -  hist.GetBinContent( kk-1 ) )/ (2.*h);
    derihist->SetBinContent( kk, der );
  }
  return derihist;
}

TH1F *getSmoothCh(TH1F &hist, double xb[]) {
  unsigned int nb = 600;
  double yi = 0.;

  TString tmpname;
  tmpname.Form("%d",rand());

	TH1F *hstSmooth = new TH1F(tmpname, "", nb, xb);

	for ( unsigned b=0; b<nb; b++ ) {
    if ( b>6 && b<nb-7 ) {
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

vector < double > getFitRangeCh( TH1F &h, int rightleftBins ) {
  vector < double > minmax;
  int minRng = 0; // Min for fitting
	int maxRng = 0; // Max for fitting

  double binMax = 0.;
  double rawbinMax = 0.;

  for ( int kk=300; kk>80; kk-- ) {// from 300 FADC backward
    if ( h.GetBinContent(kk) > 0 ) {
      if ( h.GetBinCenter(kk) < fabs(h.GetBinCenter(kk-1)) )
      {
        binMax = h.GetBinCenter(kk); // tmp FADC for VEM
        rawbinMax = kk; // Bin for tmp VEM
      }
      else 
      {
        binMax = h.GetBinCenter(kk-1);
        rawbinMax = kk; // Bin for tmp VEM
      }
      break;
    }
  }
  
  minRng = h.GetBinCenter(rawbinMax-rightleftBins);
  maxRng = h.GetBinCenter(rawbinMax+rightleftBins);

  minmax.push_back( minRng );
  minmax.push_back( maxRng );
  minmax.push_back( binMax ); 
  //minmax.push_back( rawbinMax-rightleftBins );
  //minmax.push_back( rawbinMax );

  return minmax;
}


fitcharge::fitcharge() {
	vemPosCh = 0.;
  vemPosDeri = 0.;
	getGraph = false;
	fitChOk = false;

  chisCharge = 0.;
  par0 = 0.;
  par1 = 0.;
  par2 = 0.;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	critGoodFit = 0.;
}

void fitcharge::setChCrr(TH1F &hist, const int corr, TString name) {
	unsigned int nb = 601;
	double xb [nb];

	for ( unsigned int b=1; b<nb; b++ )
		xb[b] = b; 
  xb[0] = 0;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrCh = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
    hstCrrCh->SetBinContent( b, hist.GetBinContent(b) );
}


void fitcharge::getFitCh(TH1F &hist, int nblr) {
  TString tmpname;
  tmpname.Form("%d",rand());

	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	critGoodFit = 0.; // Criterium for "good" fitting
  vector < double > rangeValues;

  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins; b++ )
    xfadc[b] = hist.GetBinCenter(b);

  TH1F *chargeSmooth = getSmoothCh(hist, xfadc);
  TH1F *chargeSmooDer = histDerivativeCh(*chargeSmooth, xfadc);
  TH1F *chargeSmooDerSmth = getSmoothCh(*chargeSmooDer, xfadc);
  chargeSmooDerSmth = getSmoothCh(*chargeSmooDerSmth, xfadc);

  rangeValues = getFitRangeCh(*chargeSmooDerSmth, nblr);
  rangXmin = rangeValues[0];
  rangXmax = rangeValues[1];

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( unsigned int b=0; b<nXbins; b++ ) 
  {
		ycnts.push_back( hist.GetBinContent(b) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( hist.GetBinCenter(b) );
  }
	
	TGraphErrors* chToFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *poly2;

  vector < double > xResid;
  vector < double > yResid;
  vector < double > errResid;
  
  poly2 = new TF1("poly2","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
  chToFit->Fit("poly2","QR");
  
  chisCharge = poly2->GetChisquare();
  ndfCharge = poly2->GetNDF();
  probCharge = poly2->GetProb();
	vemPosCh = -poly2->GetParameter(1) / (2.*poly2->GetParameter(0));
  vemPosDeri = rangeValues[2];
	critGoodFit = 5.5;
  par0 = poly2->GetParameter(0);
  par1 = poly2->GetParameter(1);
  par2 = poly2->GetParameter(2);
  fitGraphCh = (TGraphErrors*)chToFit->Clone();

  if ( (chisCharge/ndfCharge) > critGoodFit )
    vemPosCh = 0.;

	delete chToFit;
	delete poly2;
}


TGraphErrors* fitcharge::getFitGraphCh() {
	return fitGraphCh;
}


TH1F *fitcharge::getChCrr() {
	return hstCrrCh;
}


