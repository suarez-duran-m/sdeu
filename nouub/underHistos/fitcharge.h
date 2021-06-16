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
    double ndfCharge;
    double probCharge;
		bool getGraph;
		bool fitChOk;

		void getChCrr(TH1F &hist, const int corr, TString name);
		void getFitCh(TH1F &hist);

		TGraphErrors *getFitGraphCh();
		TH1F *getChCrr();

	private:
		unsigned int rangXmin;
		unsigned int rangXmax;
		unsigned int nXbins;
		double critGoodFit;

		TGraphErrors *fitGraphCh;
		TH1F *hstCrrCh;
};

#endif
