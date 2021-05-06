#ifndef READHISTOS_H
#define READHISTOS_H

#include <iostream>
//#include <TH1.h>
class TH1;
class TF1;
class TGraphErrors;
class TGH1F;
class TList;

using namespace std;

class readHistos {
  public:
		readHistos();
    double vemPosPk;
    double vemPosCh;
    double chisPeak;
    double chisCharge;
		bool getGraph;
		bool fitChOk;
		bool fitPkOk;

		void getPkCrrOst (TH1F &hist, const int offset, TString name);
		void getChCrrOst (TH1F &hist, const int offset, TString name);

		void getFitPk(TH1F &hist, const double frac, const int fstbinFit, const double vempk);
		void getFitCh(TH1F &hist, const double frac, const int fstbinFit, const double vemch);
		void resetHstCrrs();

		TGraphErrors *getFitGraphPk();
		TGraphErrors *getFitGraphCh();
		TH1F *getPkCorr();
		TH1F *getChCorr();

	private:
		unsigned int emPkb;
		unsigned int emPkc;
		unsigned int rangXmin;
		unsigned int rangXmax;
		unsigned int nXbins;
		bool checkMax;
		double critGoodFit;

		TGraphErrors *fitGraphPk;
		TGraphErrors *fitGraphCh;
		TH1F *hstCrrPk;
		TH1F *hstCrrCh;
};

#endif
