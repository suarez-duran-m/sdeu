#include <iostream>
#include <time.h>
#include <vector>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "fitpeak.h"

// =======================
// *** Local Functions ***
// =======================

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


vector < double > minRangMax( TH1F &h )
{
  vector < double > minmax;
  double tmpMin = 0.;
  double tmpMax = 0.;

  double tmpBinMin = 0.;
  double tmpBinMax = 0.;
  double rawbinMin = 0;
  double rawbinMax = 0;

  for ( int kk=77; kk>27; kk-- ) // from 300 FADC backward
    if ( h.GetBinContent(kk) < 0 )
    {
      if ( tmpBinMax < fabs(h.GetBinContent( kk ) ) )
      {
        tmpBinMax = fabs(h.GetBinContent(kk));
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
  for ( int kk=rawbinMax; kk>0; kk-- )
  {
    if ( tmpBinMin < h.GetBinContent( kk ) )
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


fitpeak::fitpeak() 
{
  srand (time(NULL));
	vemPosPk = 0.;
	getGraph = false;
	fitPkOk = false;

  chisPeak = 0.;
  ndfPeak = 0.;
  par0 = 0.;
  par1 = 0.;
  par2 = 0.;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	checkMax = false;
	critGoodFit = 0.;
}

void fitpeak::getCrr(TH1F &hist, const int corr, TString name) 
{
	unsigned int nb = 151;
	double xb [nb];

	for ( unsigned int b=1; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b) - corr; // Aplying correction for calib-baseline
  xb[0] = 0;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrPk = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b < 101 )
      hstCrrPk->SetBinContent(b, hist.GetBinContent(b)/hist.GetBinWidth(b));
    else
      hstCrrPk->SetBinContent(b, (4.*hist.GetBinContent(b))/hist.GetBinWidth(b));
  }
}


void fitpeak::getFitPk(TH1F &hist) 
{
	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	critGoodFit = 0.; // Criterium for "good" fitting

  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins; b++ )
    xfadc[b] = hist.GetBinCenter(b);

  vector < double > rangeValues;
  TH1F *peakSmooth = getSmooth(hist, xfadc);
  TH1F *peakSmooDer = histDerivative(*peakSmooth, xfadc);
  TH1F *peakSmooDerSmth = getSmooth(*peakSmooDer, xfadc);

  rangeValues = minRangMax( *peakSmooDerSmth );
  rangXmin = rangeValues[3]; //1.1*tmp[0];
  rangXmax = rangeValues[1]; //1.3*tmp[1];
  //double binMax = rangeValues[2];
  //cerr << "MSDpkfit: " << rangXmin << " " << rangXmax << endl;

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( unsigned int b=0; b<nXbins; b++ )
  {
		ycnts.push_back( hist.GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( hist.GetBinCenter(b+1) );
	}

  TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *fitFcn = new TF1("fitFcn","[0]*x*x+[1]*x+[2]",rangXmin,rangXmax);
	chFit->Fit("fitFcn","QR");

  chisPeak = fitFcn->GetChisquare();
  ndfPeak = fitFcn->GetNDF();
  probPeak = fitFcn->GetProb();
  vemPosPk = -fitFcn->GetParameter(1) / (2.*fitFcn->GetParameter(0));
  par0 = rangXmin;//fitFcn->GetParameter(0);
  par1 = rangXmax;//fitFcn->GetParameter(1);
  par2 = fitFcn->GetParameter(2);
  critGoodFit = 3.5; // OffLine Default value.

  fitGraphPk = (TGraphErrors*)chFit->Clone();
  if ( (chisPeak/ndfPeak) < critGoodFit )
    fitPkOk = true;

  if ( (chisPeak/ndfPeak) > critGoodFit || ndfPeak==0 )
    vemPosPk = 0.;
	delete chFit;
	delete fitFcn;
}


TGraphErrors* fitpeak::getFitGraphPk() {
	return fitGraphPk;
}


TH1F *fitpeak::getPkCorr() {
	return hstCrrPk;
}
