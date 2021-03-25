#include <iostream>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "readHistos.h"

// =======================
// *** Local Functions ***
// =======================
double fitFunction(double *x0, double *par) {
	const double x = x0[0]; //1./x0[0];
  //const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]), 2 ) );
  //const double expo = exp( par[3] - par[4]*x );
  const double lognormal = exp(par[0])/x*exp(-0.5*pow(((log(x)-log(par[1]))/par[2]),2));
  const double expo = exp(par[3]-1./par[4]*x);
	return lognormal+expo;
}


double peakMax(TF1* f) {
	TCanvas* c = new TCanvas("c","c",800,600);
  TGraph* der = (TGraph*)f->DrawDerivative("goff");
  double * x = der->GetX();
  double * yd1 = der->GetY();
  int n = der->GetN();
  double d0x = 0;
  double d0y = 1000;
  for (int j = 1; j < n; ++j) {
    const double yd2 = f->Derivative2(x[j]);
    if (yd2 <0 && d0y > fabs(yd1[j])) {
      d0y = fabs(yd1[j]);
      d0x = x[j];
    }
  }
  delete der;
  delete c;
  return d0x;
}


readHistos::readHistos() {
	vemPos = 0.;
	getGraph = false;
	fitChOk = false;
	fitPkOk = false;

}


void readHistos::getFullFit(TH1F &hist, const bool ifch, const double frac, const int fstbinFit, const double bs) {
	unsigned int emPkb = 0; // Bin for EM peak
	unsigned int emPkc = 0; // Counts for EM peak
	unsigned int rangXmin = 0;
	double tmpBin = 0.;

	for ( unsigned b=10; b<400; b++ )
		if ( hist.GetBinContent(b) > emPkc ) {
			emPkc = hist.GetBinContent(b);
			emPkb = b;
		}

	tmpBin = ( hist.GetBinCenter(emPkb + fstbinFit+1) - hist.GetBinCenter(emPkb + fstbinFit) )/2.0
		+ hist.GetBinCenter(emPkb + fstbinFit);
	//if ( ifch )
		//tmpBin -= bs;
	rangXmin = hist.GetBinCenter(emPkb + fstbinFit) - hist.GetBinCenter(0); //tmpBin;

	unsigned int rangXmax = 0;
	unsigned int nXbins = hist.GetXaxis()->GetNbins();
	bool checkMax = true;
	TString parName;

	vector < double > xbins;
	vector < double > ycnts;
	vector < double > yerrs;
/*
	ycnts.push_back( hist.GetBinContent(0) );
	yerrs.push_back( sqrt( ycnts[0] ) );
	tmpBin = hist.GetBinCenter(0) - (hist.GetBinCenter(1) - hist.GetBinCenter(0)) / 2.0;
	if ( ifch )
		tmpBin -= bs;
	xbins.push_back( tmpBin );
	cerr << tmpBin << " " << hist.GetBinCenter(0) << " " << hist.GetBinCenter(1) << " " << bs << endl;
*/
	for( unsigned int b=1; b<nXbins+1; b++ ) {
		ycnts.push_back( hist.GetBinContent( b ) );
		yerrs.push_back( sqrt( ycnts[b-1] ) );
		//tmpBin = ( hist.GetBinCenter(b+1) - hist.GetBinCenter(b) )/2.0 
			//	+ hist.GetBinCenter(b);
		//if ( ifch )
			//tmpBin -= bs;
		tmpBin = hist.GetBinCenter(b) - hist.GetBinCenter(0);
		xbins.push_back( tmpBin );
		if ( hist.GetBinContent( nXbins-b ) > frac*emPkc && checkMax ) {
			//tmpBin = ( hist.GetBinCenter(nXbins-b-1) - hist.GetBinCenter(nXbins-b) )/2.0
				//+ hist.GetBinCenter(nXbins-b);
			//if ( ifch )
				//tmpBin -= bs;
			rangXmax = hist.GetBinCenter(nXbins-b) - hist.GetBinCenter(0); //tmpBin;
			checkMax = false;
		}
	}
/*
	tmpBin = hist.GetBinCenter(nXbins - 1) + (hist.GetBinCenter(nXbins - 1) - hist.GetBinCenter(nXbins - 2)) / 2.0;
	if ( ifch )
		tmpBin -= bs;
	xbins.push_back( tmpBin );
*/
	//cerr << emPkc << " " << frac*emPkc << endl;
	//cerr << rangXmin << " " << rangXmax << endl;
	
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);
	double mean = (rangXmin+rangXmax)/2.;

	if ( ifch )
		fitFcn->SetParameters(11.61, mean, 0.3, 7.8, 45.5);
	else
		fitFcn->SetParameters(11, mean, 0.3, 6, mean/4);

	chFit->Fit("fitFcn","QR");

	if ( fitFcn->GetChisquare()/fitFcn->GetNDF() < 10 ) {
		if ( ifch )
			fitChOk = true;
		else 
			fitPkOk = true;
		vemPos = peakMax(fitFcn);
		if ( getGraph ) {
			chFit->SetLineColor(kBlue);
			if ( ifch )
				parName = "Charge-VEM Pos.: ";
			else
				parName = "Peak-VEM Pos.: ";
			parName += to_string( vemPos ).substr(0,7) + " +/- " + to_string( fitFcn->GetParError(1) ).substr(0,4);
			parName += "\nXi^2/NDF = " + to_string( fitFcn->GetChisquare()/fitFcn->GetNDF() ).substr(0,4);
			chFit->SetTitle( parName );
			fitGraph =  (TGraphErrors*)chFit->Clone();
		}
	}
	else {
		fitChOk = false;
		fitPkOk = false;
		vemPos = 0.;
	}
/*
	if ( ifch ) {
		fitGraph =  (TGraphErrors*)chFit->Clone();
		fitChOk = true;
		//for (int b=0; b<5; b++ )
			//cerr << fitFcn->GetParameter(b) << endl;
	}
*/
	delete chFit;
	delete fitFcn;
}


TGraphErrors* readHistos::getFitGraph() {
	return fitGraph;
}
