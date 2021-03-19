#include <iostream>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>

#include "readHistos.h"

double fitFunction(double *x, double *par);
double peakMax(TF1* f);


readHistos::readHistos() {
	binVm = 0;
	area = 0.;
	peak = 0.;
}


unsigned int readHistos::getBinVem(TH1F &hist) {
	int lclMn = hist.GetMaximum();
	int lclPk = 0;
	int bc = 0;
	binVm = 0;
	int count = 0;

	int refbin = hist.GetMaximumBin() + 2;
	for(int i=hist.GetMaximumBin()+30; i<hist.GetXaxis()->GetNbins(); i++) {
		bc = hist.GetBinContent(i);
		if ( lclMn > bc && fabs(refbin-i) < 10) {
			lclMn = bc;
			lclPk = lclMn;
			refbin = i;
		}
		if ( lclPk < bc ) {
			lclPk = bc;
			binVm = i;
		}
	}
	return binVm;
}


double readHistos::getFitVem(TH1F &hist, const bool ifch) {
	double mean = 0.;

	if ( ifch ) {
		hist.Fit("gaus","","", hist.GetXaxis()->GetBinCenter(binVm)-(8*50), 
				hist.GetXaxis()->GetBinCenter(binVm) + (8*50) // 8 for binwidth; 50 bins fordward
				);
		gaussFit = (TF1*)hist.GetListOfFunctions()->FindObject("gaus");
		if ( gaussFit->GetParError(1) < 40 && gaussFit->GetParameter(1) > 0 ) {
			mean = gaussFit->GetParameter(1);
			area = mean;
			gaussFit->SetLineColor(kRed);
			//cerr << "GaussMean: " << gaussFit->GetParError(1) << endl;
		}
	}
	else if ( !ifch ) {
		hist.Fit("gaus","Q","",
				hist.GetXaxis()->GetBinCenter(binVm) - (4*12),
				hist.GetXaxis()->GetBinCenter(binVm) + (4*12)
				);
		gaussFit = (TF1*)hist.GetListOfFunctions()->FindObject("gaus");
		if ( gaussFit->GetParError(1) < 40 && gaussFit->GetParameter(1) > 0 ) {
			mean = gaussFit->GetParameter(1);
			peak = mean;
		}
	}
	else
		mean = 0.;
	return mean;
}


TGraphErrors *readHistos::getFullFit(TH1F &hist, const bool ifch, const double frac, const int fstbinFit) {
	double max = 0.;
	unsigned int emPkb = hist.GetMaximumBin(); // Bin for EM peak
	unsigned int emPkc = hist.GetMaximum(); // Counts for EM peak
	unsigned int rangXmin = 8*(emPkb + fstbinFit);
	unsigned int rangXmax = 0;
	unsigned int nXbins = hist.GetXaxis()->GetNbins();
	bool checkMax = true;

	vector < double > xbins;
	vector < double > ycnts;
	vector < double > yerrs;

	for( unsigned int b=1; b<nXbins+1; b++ ) {
		ycnts.push_back( hist.GetBinContent( b ) );
		yerrs.push_back( sqrt( ycnts[b-1] ) );
		xbins.push_back( hist.GetBinCenter(b) - 4 ); // -4 to fit the first bin on 0.
		if ( hist.GetBinContent( nXbins-b ) > frac*emPkc && checkMax ) {
			rangXmax = hist.GetBinLowEdge( nXbins-b );
			checkMax = false;
		}
	}
	
	TGraphErrors* chFit = new TGraphErrors( xbins.size(), &xbins.front(),
			&ycnts.front(), 0, &yerrs.front() );

	TF1 *fitFcn = new TF1("fitFcn", fitFunction, rangXmin, rangXmax, 5);

	fitFcn->SetParameters(801970, 1790.76, 3.41137, 3.75441, -1772.71); // Par. coming from 1st Try for Entry 55->49 Station 1223
	chFit->Fit("fitFcn","QR");
	chFit->SetLineColor(kBlue);
	//max = peakMax(fitFcn);
	hist.SetLineColor(kBlue);
	TString parName = to_string( fitFcn->GetParError(1) );
	parName += " Xi^2/NDF" + to_string( fitFcn->GetChisquare()/fitFcn->GetNDF() );
	chFit->SetTitle( parName );
	//cerr << fitFcn->GetParError(1) << endl;

	//for ( int i=0; i<5; i++ )
		//cerr << i << " " << fitFcn->GetParameter(i) << endl;

	//delete chFit;
	delete fitFcn;

	return chFit;
}


double fitFunction(double *x0, double *par) {
	const double x = 1./x0[0];
  const double lognormal = par[0]*x*exp( -0.5*pow( ((log(x)+log(par[1]))*par[2]),2 ) );
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
