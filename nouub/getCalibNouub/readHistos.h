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
    double vemPos;
		bool getGraph;
		bool fitChOk;
		bool fitPkOk;

		void getFullFit(TH1F &hist, const bool ifch, const double frac, const int fstbinFit, const double bs);
		TGraphErrors *getFitGraph();

	private:
		TGraphErrors *fitGraph;
};

#endif
