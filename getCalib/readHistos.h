#ifndef READHISTOS_H
#define READHISTOS_H

#include <iostream>
//#include <TH1.h>
class TH1;
class TF1;
class TGH1F;
class TList;

using namespace std;

class readHistos {
  public:
		readHistos();
    unsigned int binVm;
		TF1 *gMean;
    unsigned int getBinVem(TH1F *hist);

		double getFitVem(TH1F *hist, bool ifch);

		double area;
		double peak;
};

#endif
