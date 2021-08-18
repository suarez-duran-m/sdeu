#ifndef FITPEAK_H
#define FITPEAK_H

#include <iostream>
//#include <TH1.h>
class TH1;
class TF1;
class TGraphErrors;
class TGH1F;
class TList;

using namespace std;

class fitpeak {
  public:
		fitpeak();
    double vemPosPk;
    double chisPeak;
    int cntvemPk;
		bool getGraph;
		bool fitPkOk;

		void getCrrSmooth(TH1F &hist, const int corr, TString name);
    int getValidHisto( TH1F &hist );

		void getFitPk(TH1F &hist, const double frac, const int fstbinFit, const double vempk);

		TGraphErrors *getFitGraphPk();
		TH1F *getPkCorrSmooth();
		TH1F *getPkCorr2();

    double getRms(vector < vector < double > > values, vector < double > average, int st);
    double getRms(vector < vector < double > > values, vector < int > average, int st);

	private:
		unsigned int emPkb;
		unsigned int emPkc;
		unsigned int rangXmin;
		unsigned int rangXmax;
		unsigned int nXbins;
		bool checkMax;
		double critGoodFit;

		TGraphErrors *fitGraphPk;
		TH1F *hstCrrSmoothPk;
		TH1F *hstCrrPk2;
};

#endif
