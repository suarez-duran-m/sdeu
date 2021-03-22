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
	const double x = 1./x0[0];
  const double lognormal = exp(par[0])*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
  const double expo = exp( par[3] - par[4]*x );
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
	
	if ( ifch )
		rangXmin = 8*(emPkb + fstbinFit); // Times 8 because fo the binwidth.
	else {
		rangXmin = 4*(emPkb + fstbinFit); //hist.GetBinCenter(0) + 4*(emPkb + fstbinFit); // bin 0 for offset and times 4 because fo the binwidth.
	}

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
		if ( ifch )
			xbins.push_back( 8*b ); //hist.GetBinCenter(b) - 4 ); // -4 to fit the first bin to 0.
		else
			xbins.push_back( 4*b ); //hist.GetBinCenter(b) - 2 ); // -2 to fit the first bin to 0.
		if ( hist.GetBinContent( nXbins-b ) > frac*emPkc && checkMax ) {
			if ( ifch )
				rangXmax = 8*(nXbins - b); //hist.GetBinLowEdge( nXbins-b );
			else
				rangXmax = 4*(nXbins - b);
			checkMax = false;
		}
	}
	
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);

	if ( ifch )
		fitFcn->SetParameters(13.38, (rangXmax-rangXmin)/2., 4.34, 4.81, -1059.76); 
	else
		fitFcn->SetParameters(12.5, (rangXmax-rangXmin)/2., 3., 5.8, -180.);

	chFit->Fit("fitFcn","QR");
	vemPos = peakMax(fitFcn);

	if ( fitFcn->GetChisquare()/fitFcn->GetNDF() < 5 && vemPos > 0) {
		if ( ifch )
			fitChOk = true;
		else
			fitPkOk = true;
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
	if ( !ifch ) {
		fitGraph =  (TGraphErrors*)chFit->Clone();
		fitPkOk = true;
		for( int b=0; b<5; b++ )
			cerr << fitFcn->GetParameter(b) << endl;
	}
	*/
		
	delete chFit;
	delete fitFcn;
}


TGraphErrors* readHistos::getFitGraph() {
	return fitGraph;
}
