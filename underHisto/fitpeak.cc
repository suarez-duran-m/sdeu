#include <iostream>

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
double fitFunctionPk(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
	return lognormal+expo;
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


fitpeak::fitpeak() {
	vemPosPk = 0.;
	getGraph = false;
	fitPkOk = false;

  chisPeak = 0.;
	emPkb = 0;
	emPkc = 0;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	checkMax = false;
	critGoodFit = 0.;
}

void fitpeak::getCrrSmooth(TH1F &hist, const int corr, TString name) {

	unsigned int nb = 150;
	double xb [nb];
  double yi = 0.;

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - corr; // Aplying correction for calib-baseline
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrSmoothPk = new TH1F(name, name, nb, xb);
	//hstCrrPk2 = new TH1F(name+"2", name+"2", nb, xb);

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
    
    //hstCrrSmoothPk->SetBinContent(b+1, yi);
    hstCrrSmoothPk->SetBinContent(b+1, hist.GetBinContent(b+1));
    //hstCrrPk2->SetBinContent(b+1, hist.GetBinContent(b+1));
  }
}


int fitpeak::getValidHisto( TH1F &hist )
{
  int yi = 0;
  int tmpcnt = 0;
  int diff = 0;
  int npks = 0;
  int slope = -1;


  for ( unsigned int b=hist.GetNbinsX()-50; b>0; b-- )
  {
    yi = hist.GetBinContent(b+1);

    diff += yi - hist.GetBinContent(b);
    tmpcnt++;
    if ( tmpcnt==5 )
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
  npks = 10;
  return npks;
}


void fitpeak::getFitPk(TH1F &hist, const double frac, const int fstbinFit, const double vempk ) {
	emPkb = 0; // Bin for EM peak
	emPkc = 0; // Counts for EM peak
	rangXmin = 0; // Min for fitting
	rangXmax = 0; // Max for fitting
	nXbins = hist.GetXaxis()->GetNbins(); // Number of bins for fitting
	checkMax = true;
	critGoodFit = 0.; // Criterium for "good" fitting

	for ( unsigned b=5; b<50; b++ ) // Set for EM peak, it works for Pk and Ch
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

	TF1 *fitFcn = new TF1("fitFcn", fitFunctionPk, rangXmin, rangXmax, 5);

	fitFcn->SetParameters(12.7, vempk, -3., 5., -266.2); //Set  init. fit par.
 
	chFit->Fit("fitFcn","QR");
	vemPosPk = peakMaxPk(fitFcn);

	critGoodFit = 10.;
  chisPeak = fitFcn->GetChisquare()/fitFcn->GetNDF();
  fitGraphPk = (TGraphErrors*)chFit->Clone();
  if ( chisPeak < critGoodFit && fitFcn->GetParameter(1) > 100. ) // mean of Gauss < 100.
    fitPkOk = true;
  else
    vemPosPk = 0.;
 
  /*
  if ( chisPeak > 4.6 && chisPeak < 5. )
  {
    cerr << "MSD: " << fitFcn->GetChisquare() << " "
      << fitFcn->GetNDF() << " "
      << fitFcn->GetParameter(1) << " "
      << vempk << endl;
  }
  */

  /*
	if ( chisPeak < critGoodFit && vemPosPk > 0) {
		fitPkOk = true;
		if ( getGraph ) {
			chFit->SetLineColor(kBlue);
			parName = "UB Peak-VEM Pos.: ";
			parName += to_string( vemPosPk ).substr(0,7) + " +/- " + to_string( fitFcn->GetParError(1) ).substr(0,4);
			parName += "\nXi^2/NDF = " + to_string( fitFcn->GetChisquare()/fitFcn->GetNDF() ).substr(0,4);
			chFit->SetTitle( parName );
			fitGraphPk =  (TGraphErrors*)chFit->Clone();
		}
	}
	else {
		fitPkOk = false;
		vemPosPk = 0.;
	}
  */

	delete chFit;
	delete fitFcn;
}


TGraphErrors* fitpeak::getFitGraphPk() {
	return fitGraphPk;
}


TH1F *fitpeak::getPkCorrSmooth() {
	return hstCrrSmoothPk;
}

TH1F *fitpeak::getPkCorr2() {
	return hstCrrPk2;
}

