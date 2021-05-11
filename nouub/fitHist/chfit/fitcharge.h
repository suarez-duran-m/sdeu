#ifndef FITCHARGE_H
#define FITCHARGE_H

#include <iostream>
//#include <TH1.h>
class TH1;
class TF1;
class TGraphErrors;
class TGH1F;
class TList;

using namespace std;

class fitcharge {
  public:
		fitcharge();
    double vemPosCh; 
    double chisCharge;
    int cntvemCh;
		bool getGraph;
		bool fitChOk;

		void getCrrSmooth(TH1F &hist, const int corr, TString name);
    int getValidHisto( TH1F &hist );

		void getFitCh(TH1F &hist, const double frac, const int fstbinFit, const double vemch);
    double getRms(vector < vector < double > > values, vector < double > average, int st);
    double getRms(vector < vector < double > > values, vector < int > average, int st);

		TGraphErrors *getFitGraphCh();
		TH1F *getChCorrSmooth();
    TH1F *getChCorr2();

	private:
		unsigned int emPkb;
		unsigned int emPkc;
		unsigned int rangXmin;
		unsigned int rangXmax;
		unsigned int nXbins;
		bool checkMax;
		double critGoodFit;

		TGraphErrors *fitGraphCh;
		TH1F *hstCrrSmoothCh;
		TH1F *hstCrrCh2;
};

#endif
