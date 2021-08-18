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
    double vemPosDeri;
    double chisCharge;
    double ndfCharge;
    double probCharge;
		int rangXmin;
		int rangXmax;
    double par0;
    double par1;
    double par2;
		bool getGraph;
		bool fitChOk;

		void setChCrr(TH1F &hist, const int corr, TString name);
		void getFitCh(TH1F &hist);

		TGraphErrors *getFitGraphCh();
		TH1F *getChCrr();

	private:
		unsigned int nXbins;
		double critGoodFit;

		TGraphErrors *fitGraphCh;
		TH1F *hstCrrCh;
};

#endif
