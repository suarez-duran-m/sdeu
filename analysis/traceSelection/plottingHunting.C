void plottingHunting(int pmtId){  
  TString filename;
  TString printname;
  if ( pmtId > 0 && pmtId < 4 ){
    filename.Form("zooTraces100binsPMT%d.root", pmtId);
    printname.Form("PMT%d", pmtId);
  }
  else if ( pmtId == 4 ){
    filename = "zooTraces100binsSPMT.root";
    printname = "SPMT";
  }
  else if ( pmtId == 5 ){
    filename = "zooTraces100binsPMTSSD.root";
    printname = "PMTSSD";
  }
  else{
    cout << "==================================================" << endl;
    cout << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
    cout << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " 
      << "4 For SPMT; " << "5 For PMTSSD" << endl;
    cout << "==================================================" << endl;
    exit(0);
  }
  cout << "You have selected " << filename << endl;
    
  auto hFile = TFile::Open(filename);
  
  auto T = (TTree*)hFile->Get("T");

  auto okhRmsDiff = (TH2D*)hFile->Get("okDiffRmsh");
	auto dischRmsfDiff = (TH2D*)hFile->Get("discDiffRmsh");

	auto oklRmsDiff = (TH2D*)hFile->Get("okDiffRmsl");
	auto disclRmsfDiff = (TH2D*)hFile->Get("discDiffRmsl");


	auto pmt1 = (TH1D*)hFile->Get("stRms1");
	auto pmt2 = (TH1D*)hFile->Get("stRms2");
	auto pmt3 = (TH1D*)hFile->Get("stRms3");

	//auto wideHist = (TH2D*)hFile->Get("wide");
	//cerr << wideHist->GetBinContent(2, 23) << endl;

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};

  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

	auto copy1 = (TH2D*)okhRmsDiff->Clone();
	copy1->Reset();
	Int_t bin = 0;
	Double_t maxVal = okhRmsDiff->GetMaximum();
	
	for(Int_t i = 1; i <= okhRmsDiff->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= okhRmsDiff->GetNbinsY(); ++j){
			bin = okhRmsDiff->GetBin(i, j);
			if(!okhRmsDiff->GetBinContent(bin))
				copy1->SetBinContent(bin, maxVal);
		}
	}
/*
	for (int i=0; i<okhRmsDiff->GetNbinsX(); i++ ) {
		cout << i << " " << okhRmsDiff->GetBinContent(i, 10) << endl;
		cout << i << " " << copy1->GetBinContent(i, 10) << endl;
	}
*/
  TCanvas *c1 = new TCanvas("c1", "2D", 1, 35,3600,2400);
  c1->cd();
  c1->SetLeftMargin(0.11);
  c1->SetRightMargin(0.135);

	copy1->SetFillStyle(3004);
	copy1->SetStats(0);
	copy1->GetZaxis()->SetRangeUser(-2.,2.);
	copy1->GetZaxis()->SetTickLength(0.);
	copy1->GetZaxis()->SetLabelSize(0.02);
	for ( int i=0; i<19; i++)
		copy1->GetYaxis()->SetBinLabel(i+1, stIds[i]);
	copy1->GetYaxis()->SetTitle("Station ID");
	copy1->GetYaxis()->SetTickLength(0.);
	copy1->GetYaxis()->SetLabelSize(.05);
	copy1->GetYaxis()->SetTitleFont(42);
	copy1->GetYaxis()->SetTitleOffset(1.3);
	copy1->GetXaxis()->SetTitle("Days since December 1st 2020");
	copy1->GetXaxis()->SetNdivisions(1020*4, "kTRUE");
	copy1->GetXaxis()->SetTickLength(0.);
 	copy1->GetXaxis()->SetLabelSize(.02);
	copy1->SetFillColor(1);
	copy1->Draw("BOX SAME");
	okhRmsDiff->GetZaxis()->SetTitle("Diff. RMS / FADC");
  okhRmsDiff->GetZaxis()->SetTickLength(0.);
  okhRmsDiff->SetStats(0);
  okhRmsDiff->GetYaxis()->SetTickLength(0.);
  plain->SetPalette(87);
  c1->SetGridy();
  c1->SetGridx();
  okhRmsDiff->Draw("COLZ1 SAME");
  c1->Print("../../plots/trOk"+printname+"Hg.pdf");


	auto copy2 = (TH2D*)oklRmsDiff->Clone();
	copy2->Reset();
	maxVal = oklRmsDiff->GetMaximum();
	
	for(Int_t i = 1; i <= oklRmsDiff->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= oklRmsDiff->GetNbinsY(); ++j){
			bin = oklRmsDiff->GetBin(i, j);
			if(!oklRmsDiff->GetBinContent(bin))
				copy2->SetBinContent(bin, maxVal);
		}
	}


  TCanvas *c2 = new TCanvas("c2", "2D", 1, 35,3600,2400);
  c2->cd();
	c2->Clear();
	c2->ResetDrawn();
  c2->SetLeftMargin(0.11);
  c2->SetRightMargin(0.135);

	copy2->SetFillStyle(3004);
	copy2->SetStats(0);
	copy2->GetZaxis()->SetRangeUser(-1.,1.);
	copy2->GetZaxis()->SetTickLength(0.);
	copy2->GetZaxis()->SetLabelSize(0.02);
	for ( int i=0; i<19; i++)
		copy2->GetYaxis()->SetBinLabel(i+1, stIds[i]);
	copy2->GetYaxis()->SetTitle("Station ID");
	copy2->GetYaxis()->SetTickLength(0.);
	copy2->GetYaxis()->SetLabelSize(.05);
	copy2->GetYaxis()->SetTitleFont(42);
	copy2->GetYaxis()->SetTitleOffset(1.3);
	copy2->GetXaxis()->SetTitle("Days since December 1st 2020");
	copy2->GetXaxis()->SetNdivisions(1020*4, "kTRUE");
	copy2->GetXaxis()->SetTickLength(0.);
 	copy2->GetXaxis()->SetLabelSize(.02);
	copy2->SetFillColor(1);
	copy2->Draw("BOX SAME");
	oklRmsDiff->GetZaxis()->SetTitle("Diff. RMS / FADC");
  oklRmsDiff->GetZaxis()->SetTickLength(0.);
  oklRmsDiff->SetStats(0);
  oklRmsDiff->GetYaxis()->SetTickLength(0.);
  plain->SetPalette(87);
  c2->SetGridy();
  c2->SetGridx();
  oklRmsDiff->Draw("COLZ1 SAME");
  c2->Print("../../plots/trOk"+printname+"Lg.pdf");

	TString pmt1m;
	TString pmt2m;
	TString pmt3m;
	TString pmt1r;
	TString pmt2r;
	TString pmt3r;

	pmt1m.Form("%.2g", pmt1->GetMean());
	pmt2m.Form("%.2g", pmt2->GetMean());
	pmt3m.Form("%.2g", pmt3->GetMean());
	pmt1r.Form("%.2g", pmt1->GetRMS());
	pmt2r.Form("%.2g", pmt2->GetRMS());
	pmt3r.Form("%.2g", pmt3->GetRMS());


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  TCanvas c3("c3", "2D", 0,0,3600,2400);
  c3.SetBottomMargin(0.18);
  c3.cd();
  pmt1->GetYaxis()->SetTitle("Counts/au");
  pmt1->GetXaxis()->SetTitle("Difference of Mean / FADC");
  pmt1->SetLineColor(880);
	pmt1->SetLineWidth(2);
	pmt1->GetXaxis()->SetRangeUser(-2,8);
  pmt1->Draw();
  pmt2->SetLineColor(2);
	pmt2->SetLineWidth(2);
  pmt2->Draw("SAME");
  pmt2->SetLineColor(kBlue-4);
	pmt3->SetLineColor(kRed);
	pmt3->SetLineWidth(2);
  pmt3->Draw("SAME");
  auto legend2 = new TLegend(.8,0.7,0.5,0.9);
  legend2->AddEntry(pmt1,"PMT1 HG");
  legend2->AddEntry((TObject*)0, "Mean: "+pmt1m+"; RMS: "+pmt1r+"", "");
  legend2->AddEntry(pmt2,"PMT2 HG");
  legend2->AddEntry((TObject*)0, "Mean: "+pmt2m+"; RMS: "+pmt2r+"", "");
  legend2->AddEntry(pmt3,"PMT3 HG");
  legend2->AddEntry((TObject*)0, "Mean: "+pmt3m+"; RMS: "+pmt3r+"", "");
  legend2->Draw();
  c3.Print("../../plots/rmsDistriPmt1851.pdf");

}
