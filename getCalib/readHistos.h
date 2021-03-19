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
    unsigned int binVm;
		TF1 *gaussFit;
    unsigned int getBinVem(TH1F &hist);

		double getFitVem(TH1F &hist, const bool ifch);
		TGraphErrors *getFullFit(TH1F &hist, const bool ifch, const double frac, const int fstbinFit);

		double area;
		double peak;

	//private:
		//double fitFunction( double *x, double *par);
};

#endif
