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
double fitFunctionCh(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
}


TH1F *histDerivativeCh(TH1F &hist, double xb[]) // Central differences
{
  int nbins = 600;
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


double peakMaxCh(TF1* f) {
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


vector < double > minRangMaxFFTCh( double x[], TH1F h )
{
  vector < double > minmax;
  double tmpMax = 0.;
  double tmpMin = 0.;
  double tmpBinMax = 0.;
  double tmpBinMin = 0.;
  TString tmpname;
  tmpname.Form("%d",rand());

  TH1 *test = 0;
  TVirtualFFT::SetTransform(0);
  test = h.FFT(test, "PH");

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
  
  TH1F *peakFFT = new TH1F(tmpname,"", 600, xc);
  for (int j = 0; j < 402; j++)
    peakFFT->SetBinContent(j + 1, hbC->GetBinContent(j)/600.);
  for (int j = 0; j < 200; j++)
    peakFFT->SetBinContent(j + 1 + 400, hbC->GetBinContent(j+400)/600.);
  
  TH1 *peakFFTDer = histDerivativeCh(*peakFFT, x);

  for ( int kk=312; kk>125; kk-- ) // from 2960 FADC backward
    if ( peakFFTDer->GetBinContent(kk) < 0 )
    {
      if ( tmpBinMax < fabs(peakFFTDer->GetBinContent( kk ) ) )
      {
        tmpBinMax = fabs(peakFFTDer->GetBinContent(kk));
        tmpMax = peakFFTDer->GetBinCenter(kk);
      }
    }
    else
    {
      tmpBinMax = peakFFTDer->GetBinCenter(kk);
      break;
    }

  int tmpneg = 0;
  for ( int kk=25; kk<130; kk++ ) // 200 FADC after 0 FADC
  {
    if ( peakFFTDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( peakFFTDer->GetBinContent( kk ) < 0 )
      if ( tmpBinMin < fabs( peakFFTDer->GetBinContent( kk ) ) )
      {
        tmpMin = peakFFTDer->GetBinCenter(kk);
        tmpBinMin = fabs( peakFFTDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  }
  minmax.push_back( tmpMin );
  minmax.push_back( tmpMax );
  minmax.push_back( tmpBinMax );

  delete test;
  delete hbC;
  return minmax;
}


vector < double > minRangMaxSmoothCh( double x[], TH1F h, bool smooth )
{
  vector < double > minmax;

	unsigned int nb = 600;
  double yi = 0.;
  double tmpMax = 0.;
  double tmpMin = 0.;
  double tmpBinMax = 0.;
  double tmpBinMin = 0.;

  TString tmpname;
  tmpname.Form("%d",rand());

	TH1F *hsmooth = new TH1F(tmpname, tmpname, nb, x);

	for ( unsigned b=0; b<nb; b++ )
  {
    if ( b > 1 && b<nb-2 )
    {
      yi = h.GetBinContent(b+1 - 2) 
        + 2*h.GetBinContent(b+1 - 1) 
        + 3*h.GetBinContent(b+1) 
        + 2*h.GetBinContent(b+1 + 1) 
        + h.GetBinContent(b+1 + 2);
      yi = yi/9.;
    }
    else
      yi = h.GetBinContent(b+1);

    hsmooth->SetBinContent(b+1, yi);
  }

  TH1 *hsmoothDer = histDerivativeCh(*hsmooth, x);

  if ( smooth )
    hsmoothDer->Smooth(700);

  for ( int kk=250; kk>125; kk-- ) // from 2200 FADC backward
  {
    if ( hsmoothDer->GetBinContent(kk) < 0 )
    {
      if ( tmpBinMax < fabs(hsmoothDer->GetBinContent( kk ) ) )
      {
        tmpBinMax = fabs(hsmoothDer->GetBinContent(kk));
        tmpMax = hsmoothDer->GetBinCenter(kk);
      }
    }
    else
    {
      tmpBinMax = hsmoothDer->GetBinCenter(kk);
      break;
    }
  }

  int tmpneg = 0;
  for ( int kk=25; kk<tmpBinMax; kk++ ) // 200 FADC after 0 FADC
  {
    if ( hsmoothDer->GetBinContent( kk ) > 0 && tmpneg == 1 )
      break;
    if ( hsmoothDer->GetBinContent( kk ) < 0 )
      if ( tmpBinMin < fabs( hsmoothDer->GetBinContent( kk ) ) )
      {
        tmpMin = hsmoothDer->GetBinCenter(kk);
        tmpBinMin = fabs( hsmoothDer->GetBinContent( kk ) );
        tmpneg = 1;
      }
  }
  minmax.push_back( tmpMin );
  minmax.push_back( tmpMax );
  minmax.push_back( tmpBinMax );

  return minmax;
}


fitcharge::fitcharge() 
{
	vemPosCh = 0.;
	getGraph = false;
	fitChOk = false;

  chisCharge = 0.;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	critGoodFit = 0.;
}


void fitcharge::getChCrr(TH1F &hist, const int corr, TString name) 
{

	unsigned int nb = 600;
	double xb [nb];

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - 20.*corr;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrCh = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
    hstCrrCh->SetBinContent( b+1, hist.GetBinContent(b+1) ); // Not smooth
}


void fitcharge::getFitCh(TH1F &hist) 
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

  vector < double > tmp;
  tmp = minRangMaxFFTCh( xfadc, hist );

  rangXmin = tmp[0];
  rangXmax = tmp[1];
  double binMax = tmp[2];

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

	TF1 *fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
	fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76); //Set  init. fit par.

	chFit->Fit("fitFcn", "QR");
  chisCharge = fitFcn->GetChisquare();
  ndfCharge = fitFcn->GetNDF();
  probCharge = fitFcn->GetProb();
	vemPosCh = peakMaxCh(fitFcn);
	critGoodFit = 2;

  if ( (chisCharge/ndfCharge) < critGoodFit )
    fitChOk = true;
  else if ( (chisCharge/ndfCharge) > critGoodFit )
  {
    rangXmax += 4;
    rangXmin += 4;
    binMax += 4;
    fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76); //Set  init. fit par.
  	chFit->Fit("fitFcn", "QR");
    chisCharge = fitFcn->GetChisquare();
    ndfCharge = fitFcn->GetNDF();
	  vemPosCh = peakMaxCh(fitFcn);
  }

  if ( (chisCharge/ndfCharge) > critGoodFit )
  {
    tmp.clear();
    tmp = minRangMaxSmoothCh( xfadc, hist, 0);
    
    rangXmin = tmp[0];
    rangXmax = tmp[1];
    double binMax = tmp[2];

    TF1 *fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76);

	  chFit->Fit("fitFcn","QR");
    chisCharge = fitFcn->GetChisquare();
    ndfCharge = fitFcn->GetNDF();
    vemPosCh = peakMaxCh(fitFcn);
  }
  
  if ( (chisCharge/ndfCharge) > critGoodFit )
  {
    tmp.clear();
    tmp = minRangMaxSmoothCh( xfadc, hist, 1 );
    
    rangXmin = 2.5*tmp[0];
    rangXmax = 1.5*tmp[1];
    double binMax = tmp[2];

    TF1 *fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(13.38, binMax, 4.34, 4.81, -1659.76);

	  chFit->Fit("fitFcn","QR");
    chisCharge = fitFcn->GetChisquare();
    ndfCharge = fitFcn->GetNDF();
    vemPosCh = peakMaxCh(fitFcn);
  }

  if ( (chisCharge/ndfCharge) > critGoodFit )
    vemPosCh = 0.;

	delete chFit;
	delete fitFcn;
}


TGraphErrors* fitcharge::getFitGraphCh() {
	return fitGraphCh;
}


TH1F *fitcharge::getChCrr() {
	return hstCrrCh;
}


