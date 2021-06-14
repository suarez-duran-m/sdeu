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
    double ndfPeak;
		bool getGraph;
		bool fitPkOk;

		void getCrr(TH1F &hist, const int corr, TString name);
		void getFitPk(TH1F &hist);

		TGraphErrors *getFitGraphPk();
		TH1F *getPkCorr();

	private:
		unsigned int rangXmin;
		unsigned int rangXmax;
		unsigned int nXbins;
		bool checkMax;
		double critGoodFit;

		TGraphErrors *fitGraphPk;
		TH1F *hstCrrPk;
};

#endif
