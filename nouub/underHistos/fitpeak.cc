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

vector < double > minRangMaxSmooth( double x[], TH1F h )
{
  vector < double > minmax;

	unsigned int nb = 150;
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
    if ( b > 1 && b<148 )
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

  TH1 *hsmoothDer = histDerivative(*hsmooth, x);
  hsmoothDer->Smooth(700);

  for ( int kk=75; kk>27; kk-- ) // from 75 FADC backward
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
  for ( int kk=7; kk<tmpBinMax; kk++ ) // 7 FADC after 0 FADC
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


fitpeak::fitpeak() 
{
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


void fitpeak::getFitPk(TH1F &hist) 
{
	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	critGoodFit = 0.; // Criterium for "good" fitting

  double xfadc[nXbins+1];
  for( unsigned int b=0; b<nXbins+1; b++ )
    xfadc[b] = hist.GetBinCenter(b+1)-.5;

  vector < double > tmp;
  tmp = minRangMaxSmooth( xfadc, hist );

  rangXmin = 1.5*tmp[0]; //0.9*tmp[0];
  rangXmax = 1.5*tmp[1]; //1.2*tmp[1];
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

  //cout << "MSD0: " << rangXmin << " " << rangXmax << " " << binMax << endl;
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

  TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
  fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par.
	chFit->Fit("fitFcn","QR");

  chisPeak = fitFcn->GetChisquare();
  ndfPeak = fitFcn->GetNDF();
  probPeak = fitFcn->GetProb();
  vemPosPk = peakMaxPk(fitFcn);
  critGoodFit = 4.;

  fitGraphPk = (TGraphErrors*)chFit->Clone();

  if ( (chisPeak/ndfPeak) < critGoodFit )
    fitPkOk = true;
  else if ( chisPeak/ndfPeak > critGoodFit )
  {    
    rangXmax = 1.3*tmp[1];

    //cout << "MSD1: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.2*tmp[0];
    rangXmax = 1.2*tmp[1];
    
    //cout << "MSD2: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.1*tmp[0]-1;
    rangXmax = 1.5*tmp[1];
    
    //cout << "MSD3: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

//    TCanvas c1;
//    c1.cd();
//    chFit->Draw();
//    c1.Print("kk.pdf");

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.2*tmp[0];
    rangXmax = 1.3*tmp[1];

    ycnts.clear();
    yerrs.clear();
    xbins.clear();
    hist.Smooth(10);
  	for( unsigned int b=0; b<nXbins; b++ )
    {
		  ycnts.push_back( hist.GetBinContent( b+1 ) );
  		yerrs.push_back( sqrt( ycnts[b] ) );
	  	xbins.push_back( hist.GetBinCenter(b+1) );
	  }
  	chFit = new TGraphErrors( xbins.size(), &xbins.front(),
	  		&ycnts.front(), 0, &yerrs.front() );

    //cout << "MSD4: " << rangXmin << " " << rangXmax << " " << binMax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
    //cout << "MSD4: " << chisPeak << " " << ndfPeak << " " << chisPeak/ndfPeak << endl;
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.1*tmp[0];
    rangXmax = 1.3*tmp[1];
    
    //cout << "MSD5: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.05*tmp[0];
    rangXmax = 0.87*tmp[1];
    
    //cout << "MSD6: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.3*tmp[0];
    rangXmax = 1.5*tmp[1];
    
    cout << "MSD7: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
  if ( chisPeak/ndfPeak > critGoodFit )
  {
    rangXmin = 1.3*tmp[0];
    rangXmax = 1.4*tmp[1];
    
    cout << "MSD8: " << rangXmin << " " << rangXmax << endl;
    TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);
    fitFcn->SetParameters(11., binMax, 0.3, 6., binMax/8.); //Set  init. fit par. 

	  chFit->Fit("fitFcn","QR");
    chisPeak = fitFcn->GetChisquare();
    ndfPeak = fitFcn->GetNDF();
    vemPosPk = peakMaxPk(fitFcn);
  }
 
  if ( (chisPeak/ndfPeak) > critGoodFit )
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
