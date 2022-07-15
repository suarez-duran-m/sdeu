#ifndef FITCHARGE_H
#define FITCHARGE_H

#include <iostream>
#include "TH1.h"
#include "TGraphErrors.h"

using namespace std;

class fitcharge {
  public:
		fitcharge(TH1F *inHist, int corr, int binsLeftRigh);
    ~fitcharge();

		int rangXmin;
		int rangXmax;
    double ndf;
    double prob;
    double par0;
    double par1;
    double par2;
    double chi2;
    double qpkPos;
    double qpkPosDeri;

		TGraphErrors *GetFitGraph();
		TH1F *GetChCrr();
    void DeleteObjts();

	private:
		TGraphErrors *fittedGraph;
		TH1F *hstCrr;
    TH1F *histSmooth;
    TH1F *histSmooDer;
    TH1F *histSmooDerSmth;
    
		void SetCrr(TH1F *inHist, int corr);
		void GetFit(int binsLeftRigh);

    TH1F *SetSmooth(TH1F *hist, double xb[]);
    TH1F *SetHistDerivative(double xb[]);
    vector<double> SetRangeFit( TH1F &hDer, int rightleftBins );
};

#endif
