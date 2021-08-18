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
	vemPosPk = 0.;
	vemPosCh = 0.;
	getGraph = false;
	fitChOk = false;
	fitPkOk = false;

  chisPeak = 0.;
  chisCharge = 0.;
	emPkb = 0;
	emPkc = 0;
	rangXmin = 0;
	rangXmax = 0;
	nXbins = 0;
	checkMax = false;
	critGoodFit = 0.;
}

void readHistos::getPkCrrOst(TH1F &hist, const int corr, TString name) {

	unsigned int nb = 150;
	double xb [nb];

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - corr;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrPk = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
		hstCrrPk->SetBinContent(b+1, hist.GetBinContent(b+1));	
}


void readHistos::getChCrrOst(TH1F &hist, const int corr, TString name) {

	unsigned int nb = 600;
	double xb [nb];

	for ( unsigned int b=0; b<nb; b++ )
		xb[b] = hist.GetBinCenter(b+1) - 20.*corr;
	xb[nb] = xb[nb-1] + hist.GetBinWidth(1);

	hstCrrCh = new TH1F(name, name, nb, xb);

	for ( unsigned b=0; b<nb; b++ )
		hstCrrCh->SetBinContent(b+1, hist.GetBinContent(b+1));	
}


void readHistos::getFitPk(TH1F &hist, const double frac, const int fstbinFit, const double vempk ) {
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

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);

	fitFcn->SetParameters(11., vempk, 0.3, 6., vempk/8.); //Set  init. fit par.

	chFit->Fit("fitFcn","QR");
	vemPosPk = peakMax(fitFcn);

	critGoodFit = 7;
  chisPeak = fitFcn->GetChisquare()/fitFcn->GetNDF();

	if ( chisPeak < critGoodFit && vemPosPk > 0) {
		fitPkOk = true;
		if ( getGraph ) {
			chFit->SetLineColor(kBlue);
			parName = "UUB Peak-VEM Pos.: ";
			parName += to_string( vemPosPk ).substr(0,7) + " +/- " + to_string( fitFcn->GetParError(1) ).substr(0,4);
			parName += "\nXi^2/NDF = " + to_string( fitFcn->GetChisquare()/fitFcn->GetNDF() ).substr(0,4);
			chFit->SetTitle( parName );
			fitGraphPk =  (TGraphErrors*)chFit->Clone();
		}
	}
	else {
		fitChOk = false;
		fitPkOk = false;
		vemPosPk = 0.;
	}

	delete chFit;
	delete fitFcn;
}


void readHistos::getFitCh(TH1F &hist, const double frac, const int fstbinFit, const double vemch ) {
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

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);

	fitFcn->SetParameters(11.61, vemch, 0.3, 7.8, 45.5); //Set  init. fit par.

	chFit->Fit("fitFcn", "QR");
	vemPosCh = peakMax(fitFcn);
	critGoodFit = 15;

  chisCharge = fitFcn->GetChisquare()/fitFcn->GetNDF();
	
	if ( chisCharge < critGoodFit && vemPosCh > 0) {
		fitChOk = true;
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
		fitPkOk = false;
		vemPosCh = 0.;
	}

	//fitChOk = true;

	delete chFit;
	delete fitFcn;
}


TGraphErrors* readHistos::getFitGraphPk() {
	return fitGraphPk;
}


TGraphErrors* readHistos::getFitGraphCh() {
	return fitGraphCh;
}


TH1F *readHistos::getPkCorr() {
	return hstCrrPk;
}


TH1F *readHistos::getChCorr() {
	return hstCrrCh;
}


void readHistos::resetHstCrrs() {
	hstCrrPk->Reset();
	hstCrrPk = new TH1F();

	hstCrrCh->Reset();
	hstCrrCh = new TH1F();
}
