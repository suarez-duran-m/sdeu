#include <iostream>

#include <TH1.h>
#include <TF1.h>
#include <TH1F.h>
#include <TList.h>

#include "readHistos.h"

readHistos::readHistos() {
	binVm = 0;
	gMean = NULL;
	area = 0.;
	peak = 0.;
}

unsigned int readHistos::getBinVem(TH1F *hist) {
	int lclMn = hist->GetMaximum();
	int lclPk = 0;
	int bc = 0;
	binVm = 0;

	int refbin = hist->GetMaximumBin() + 2;
	for(int i=hist->GetMaximumBin()+10; i<hist->GetXaxis()->GetNbins(); i++) {
		bc = hist->GetBinContent(i);
		if ( lclMn > bc && fabs(refbin-i) < 10) {
			lclMn = bc;
			lclPk = lclMn;
			refbin = i;
		}
		if ( lclPk < bc ) {
			lclPk = bc;
			binVm = i;
		}
	}
	return binVm;
}


double readHistos::getFitVem(TH1F *hist, bool ifch) {
	gMean = NULL;
	double mean = 0.;

	if ( ifch ) {
		hist->Fit("gaus","","", hist->GetXaxis()->GetBinCenter(binVm)-(8*50),
				hist->GetXaxis()->GetBinCenter(binVm) + (8*50)
				);
		gMean = (TF1*)hist->GetListOfFunctions()->FindObject("gaus");
		if ( gMean->GetParError(1) < 40 && gMean->GetParameter(1) > 0 ) {
			mean = gMean->GetParameter(1);
			area = mean;
		}
	}
	else if ( !ifch ) {
		hist->Fit("gaus","","",
				hist->GetXaxis()->GetBinCenter(binVm) - (4*12),
				hist->GetXaxis()->GetBinCenter(binVm) + (4*12)
				);
		gMean = (TF1*)hist->GetListOfFunctions()->FindObject("gaus");
		if ( gMean->GetParError(1) < 40 && gMean->GetParameter(1) > 0 ) {
			mean = gMean->GetParameter(1);
			peak = mean;
		}
	}
	else
		mean = 0.;
	return mean;
}
