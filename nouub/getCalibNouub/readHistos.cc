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


void readHistos::getFullFit(TH1F &hist, const bool ifch, const double frac, const int fstbinFit) {
	unsigned int emPkb = 0; // Bin for EM peak
	unsigned int emPkc = 0; // Counts for EM peak
	unsigned int rangXmin = 0;

	for ( unsigned b=10; b<400; b++ )
		if ( hist.GetBinContent(b) > emPkc ) {
			emPkc = hist.GetBinContent(b);
			emPkb = b;
		}

	rangXmin =  emPkb + fstbinFit; //hist.GetBinLowEdge( emPkb )+ fstbinFit;

	unsigned int rangXmax = 0;
	unsigned int nXbins = hist.GetXaxis()->GetNbins();
	bool checkMax = true;
	TString parName;

	vector < double > xbins;
	vector < double > ycnts;
	vector < double > yerrs;

	for( unsigned int b=1; b<nXbins+1; b++ ) {
		ycnts.push_back( hist.GetBinContent( b ) );
		yerrs.push_back( sqrt( ycnts[b-1] ) );
		xbins.push_back( b );
		if ( hist.GetBinContent( nXbins-b ) > frac*emPkc && checkMax ) {
			rangXmax = nXbins - b;
			checkMax = false;
		}
	}
	
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);
	double mean = (rangXmin+rangXmax)/2.;

	if ( ifch )
		//fitFcn->SetParameters(13.81, (rangXmin+rangXmax)/2., 23.02, -33.7, -46272.5);
		fitFcn->SetParameters(11.61, mean, 0.3, 7.8, 45.5);
	else
		fitFcn->SetParameters(11, mean, 0.3, 6, mean/4);
		//fitFcn->SetParameters(12.22, (rangXmin+rangXmax)/2., 5.2, -0.13, -622.935);

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

/*	if ( ifch ) {
		fitGraph =  (TGraphErrors*)chFit->Clone();
		fitChOk = true;
		for (int b=0; b<5; b++ )
			cerr << fitFcn->GetParameter(b) << endl;
	}
	*/	
	delete chFit;
	delete fitFcn;
}


TGraphErrors* readHistos::getFitGraph() {
	return fitGraph;
}
