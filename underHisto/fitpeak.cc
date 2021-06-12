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

#include "fitpeak.h"

// =======================
// *** Local Functions ***
// =======================
double fitFunctionPk(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
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

double peakMaxPk(TF1* f) {
	TCanvas* c = new TCanvas("c","c",800,600);
  TGraph* der = (TGraph*)f->DrawDerivative("goff");
  double * x = der->GetX();
  double * yd1 = der->GetY();
  int n = der->GetN();
  double d0x = 0;
  double d0y = 1000;
  for (int j = 1; j < n; ++j){
    const double yd2 = f->Derivative2(x[j]);
    if (yd2 <0 && d0y > fabs(yd1[j])){
      d0y = fabs(yd1[j]);
      d0x = x[j];
    }
  }
  delete der;
  delete c;
  return d0x;
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


fitpeak::fitpeak() {
  srand (time(NULL));
	vemPosPk = 0.;
	getGraph = false;
	fitPkOk = false;

  chisPeak = 0.;
  ndfPeak = 0.;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	checkMax = false;
	critGoodFit = 0.;
}

void fitpeak::getCrr(TH1F &hist, const int corr, TString name) 
{
	unsigned int nb = 150;
	double xb [nb];

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - corr; // Aplying correction for calib-baseline
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrPk = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
    hstCrrPk->SetBinContent(b+1, hist.GetBinContent(b+1));
}


void fitpeak::getFitPk(TH1F &hist, const double vempk ) 
{
  TString tmpname;
  tmpname.Form("%d",rand());

	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	critGoodFit = 0.; // Criterium for "good" fitting

  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins+1; b++ )
    xfadc[b] = hist.GetBinCenter(b+1);


  TH1 *test = 0;
  TVirtualFFT::SetTransform(0);
  test = hist.FFT(test, "PH");

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
  
  double xc[151];
  for (int j = 0; j < 102; j++)
    xc[j] = 4.*j - 12;
  for (int j = 0; j < 51; j++)
    xc[100 + j] = 100*4. + 4.*4.*j;
  
  TH1F *peakFFT = new TH1F(tmpname,"", 150, xc);
  for (int j = 0; j < 102; j++)
    peakFFT->SetBinContent(j + 1, hbC->GetBinContent(j)/150.);
  for (int j = 0; j < 50; j++)
    peakFFT->SetBinContent(j + 1 + 100, hbC->GetBinContent(j+100)/150.);
  
  TH1 *peakFFTDer = histDerivative(*peakFFT, xfadc);


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

	TString parName; // For Fitted plot title

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

  TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(12.7, vempk, -3., 5., -266.2); //Set  init. fit par. 

	chFit->Fit("fitFcn","QR");
  chisPeak = fitFcn->GetChisquare();
  ndfPeak = fitFcn->GetNDF();
  vemPosPk = peakMaxPk(fitFcn);
  critGoodFit = 5.;

  fitGraphPk = (TGraphErrors*)chFit->Clone();
  if ( (chisPeak/ndfPeak) < critGoodFit )
    fitPkOk = true;
  else
    vemPosPk = 0.;

  delete test;
  delete hbC;
	delete chFit;
	delete fitFcn;
}


TGraphErrors* fitpeak::getFitGraphPk() {
	return fitGraphPk;
}


TH1F *fitpeak::getPkCorr() {
	return hstCrrPk;
}
