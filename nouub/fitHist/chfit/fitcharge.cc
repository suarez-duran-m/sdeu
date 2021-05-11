#include <iostream>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

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


fitcharge::fitcharge() {
	vemPosCh = 0.;
	getGraph = false;
	fitChOk = false;

  chisCharge = 0.;
	emPkb = 0;
	emPkc = 0;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	checkMax = false;
	critGoodFit = 0.;
}


void fitcharge::getCrrSmooth(TH1F &hist, const int corr, TString name) {

	unsigned int nb = 600;
	double xb [nb];
  double yi = 0.;

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - 20.*corr;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrSmoothCh = new TH1F(name, name, nb, xb);
  hstCrrCh2 = new TH1F(name+"2", name+"2", nb, xb);

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

    hstCrrSmoothCh->SetBinContent( b+1, yi );
    //hstCrrSmoothCh->SetBinContent( b+1, hist.GetBinContent(b+1) );
		//hstCrrCh2->SetBinContent(b+1, hist.GetBinContent(b+1));
  }
}


int fitcharge::getValidHisto( TH1F &hist )
{
  int yi = 0;
  int tmpcnt = 0;
  int diff = 0;
  int npks = 0;
  int slope = -1;


  for ( unsigned int b=hist.GetNbinsX()-200; b>0; b-- )
  {
    yi = hist.GetBinContent(b+1);

    diff += yi - hist.GetBinContent(b);
    tmpcnt++;
    if ( tmpcnt==30)
    {
      if ( diff > slope )
      {
        npks++;
        slope *= -1;
      }
      diff = 0;
      tmpcnt = 0; 
    }
  }
  npks = 5; 
  return npks;
}


void fitcharge::getFitCh(TH1F &hist, const double frac, const int fstbinFit, const double vemch ) {
	emPkb = 0; // Bin for EM peak
	emPkc = 0; // Counts for EM peak
	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	checkMax = true;
	critGoodFit = 0.; // Criterium for "good" fitting

	for ( unsigned b=0; b<50; b++ ) // Set for EM peak, it works for Pk and Ch
		if ( hist.GetBinContent(b) > emPkc ) {
			emPkc = hist.GetBinContent(b);
			emPkb = hist.GetBinLowEdge(b);
		}

	rangXmin = emPkb + fstbinFit*hist.GetXaxis()->GetBinWidth(1);
	TString parName; // For Fitted plot title

	vector < double > xbins; // X bins for fit-function
	vector < double > ycnts; // Y counts for fit-function
	vector < double > yerrs; // Y errors for fit-function

	for( unsigned int b=0; b<nXbins; b++ ) {
		ycnts.push_back( hist.GetBinContent( b+1 ) );
		yerrs.push_back( sqrt( ycnts[b] ) );
		xbins.push_back( hist.GetBinCenter(b+1) );
		if ( hist.GetBinContent( nXbins - b+1 ) > frac*emPkc && checkMax ) { //Searching for Max-fitting
			rangXmax = hist.GetBinCenter(nXbins - b+1);
			checkMax = false;
		}
	}
	
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

	TF1 *fitFcn = new TF1("fitFcn", fitFunctionCh, rangXmin, rangXmax, 5);

	fitFcn->SetParameters(11.61, vemch, 0.3, 7.8, 45.5); //Set  init. fit par.

	chFit->Fit("fitFcn", "QR");
	vemPosCh = peakMaxCh(fitFcn);
	critGoodFit = 15;

  chisCharge = fitFcn->GetChisquare()/fitFcn->GetNDF();
  cntvemCh = -1;
	
	if ( chisCharge < critGoodFit && vemPosCh > 0) {
		fitChOk = true;
    cntvemCh = fitFcn->Eval(vemPosCh);
		if ( getGraph ) {
			chFit->SetLineColor(kBlue);
			parName = "UUB Charge-VEM Pos.: ";
			parName += to_string( vemPosCh ).substr(0,7) + " +/- " + to_string( fitFcn->GetParError(1) ).substr(0,4);
			parName += "\nXi^2/NDF = " + to_string( fitFcn->GetChisquare()/fitFcn->GetNDF() ).substr(0,4);
			chFit->SetTitle( parName );
			fitGraphCh = (TGraphErrors*)chFit->Clone();
		}
	}
	else {
		fitChOk = false;
		vemPosCh = 0.;
	}

	//fitChOk = true;

	delete chFit;
	delete fitFcn;
}


TGraphErrors* fitcharge::getFitGraphCh() {
	return fitGraphCh;
}


TH1F *fitcharge::getChCorrSmooth() {
	return hstCrrSmoothCh;
}


TH1F *fitcharge::getChCorr2() {
	return hstCrrCh2;
}


double fitcharge::getRms(vector < vector < double > > values, vector < double > average, int st)
{
  double rms = 0.;
  for( unsigned int val = 0; val<values[st].size(); val++ )
    rms += (average[st] - values[st][val])*(average[st] - values[st][val]);

  return sqrt(rms/values[st].size());
}


double fitcharge::getRms(vector < vector < double > > values, vector < int > average, int st)
{
  double rms = 0.;
  for( unsigned int val = 0; val<values[st].size(); val++ )
    rms += (average[st] - values[st][val])*(average[st] - values[st][val]);

  return sqrt(rms/values[st].size());
}
