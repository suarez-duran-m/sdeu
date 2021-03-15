void plotting(int pmtId){  
  TString filename;
  TString printname;
  if ( pmtId > 0 && pmtId < 4 ){
    filename.Form("calibHistPMT%d.root", pmtId);
    printname.Form("PMT%d", pmtId);
  }
  else if ( pmtId == 4 ){
    filename = "calibHistSPMT.root";
    printname = "SPMT";
  }
  else if ( pmtId == 5 ){
    filename = "calibHistPMTSSD.root";
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
  
  auto charge = (TH2F*)hFile->Get("hCh");
  auto peak = (TH2F*)hFile->Get("hPk");
  auto ap = (TH2F*)hFile->Get("hap");

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};

  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

	auto copy1 = (TH2D*)charge->Clone();
	copy1->Reset();
	Int_t bin = 0;
	Double_t maxVal = charge->GetMaximum();
	
	for(Int_t i = 1; i <= charge->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= charge->GetNbinsY(); ++j){
			bin = charge->GetBin(i, j);
			if(!charge->GetBinContent(bin))
				copy1->SetBinContent(bin, maxVal);
		}
	}

  TCanvas *c1 = new TCanvas("c1", "2D", 1, 35,3600,2400);
  c1->cd();
  c1->SetLeftMargin(0.11);
  c1->SetRightMargin(0.135);

	copy1->SetFillStyle(3004);
	copy1->SetStats(0);
	if ( pmtId == 3 )
		copy1->GetZaxis()->SetRangeUser(10000.,30000.);
	else if ( pmtId == 4 )
		copy1->GetZaxis()->SetRangeUser(16000.,34000.);
	else if ( pmtId == 5 )
		copy1->GetZaxis()->SetRangeUser(1000.,3200.);
	else
		copy1->GetZaxis()->SetRangeUser(1000.,2000.);
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
	charge->GetZaxis()->SetTitle("Position of VEM-Charge / FADC");
  charge->GetZaxis()->SetTickLength(0.);
  charge->GetZaxis()->SetTitleOffset(1.4);
  charge->SetStats(0);
  charge->GetYaxis()->SetTickLength(0.);
  plain->SetPalette(87);
  c1->SetGridy();
  c1->SetGridx();
  charge->Draw("COLZ1 SAME");
  c1->Print("../plots/charge"+printname+"Hg.pdf");


	auto copy2 = (TH2D*)peak->Clone();
	copy2->Reset();
	bin = 0;
	maxVal = peak->GetMaximum();
	
	for(Int_t i = 1; i <= peak->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= peak->GetNbinsY(); ++j){
			bin = peak->GetBin(i, j);
			if(!peak->GetBinContent(bin))
				copy2->SetBinContent(bin, maxVal);
		}
	}

  TCanvas *c2 = new TCanvas("c2", "2D", 1, 35,3600,2400);
  c2->cd();
  c2->SetLeftMargin(0.11);
  c2->SetRightMargin(0.135);

	copy2->SetFillStyle(3004);
	copy2->SetStats(0);
	if ( pmtId == 4 )
		copy2->GetZaxis()->SetRangeUser(400.,500.);
	else if ( pmtId == 5 )
		copy2->GetZaxis()->SetRangeUser(300.,1000.);
	else
		copy2->GetZaxis()->SetRangeUser(200.,600.);
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
	peak->GetZaxis()->SetTitle("Position of VEM-Peak / FADC");
  peak->GetZaxis()->SetTickLength(0.);
	peak->GetZaxis()->SetTitleOffset(1.3);
  peak->SetStats(0);
  peak->GetYaxis()->SetTickLength(0.);
  plain->SetPalette(87);
  c2->SetGridy();
  c2->SetGridx();
  peak->Draw("COLZ1 SAME");
  c2->Print("../plots/peak"+printname+"Hg.pdf");


	auto copy3 = (TH2D*)ap->Clone();
	copy3->Reset();
	bin = 0;
	maxVal = ap->GetMaximum();
	
	for(Int_t i = 1; i <= ap->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= ap->GetNbinsY(); ++j){
			bin = ap->GetBin(i, j);
			if(!ap->GetBinContent(bin))
				copy3->SetBinContent(bin, maxVal);
		}
	}

  TCanvas *c3 = new TCanvas("c3", "2D", 1, 35,3600,2400);
  c3->cd();
  c3->SetLeftMargin(0.11);
  c3->SetRightMargin(0.135);

	copy3->SetFillStyle(3004);
	copy3->SetStats(0);
	if ( pmtId==3 )
		copy3->GetZaxis()->SetRangeUser(30.,80.);
	else if ( pmtId==4 )
		copy3->GetZaxis()->SetRangeUser(10.,50.);
	else if ( pmtId == 5 )
		copy3->GetZaxis()->SetRangeUser(1.5,7.5);
	else 
		copy3->GetZaxis()->SetRangeUser(2.,5.);
	copy3->GetZaxis()->SetTickLength(0.);
	copy3->GetZaxis()->SetLabelSize(0.02);
	for ( int i=0; i<19; i++)
		copy3->GetYaxis()->SetBinLabel(i+1, stIds[i]);
	copy3->GetYaxis()->SetTitle("Station ID");
	copy3->GetYaxis()->SetTickLength(0.);
	copy3->GetYaxis()->SetLabelSize(.05);
	copy3->GetYaxis()->SetTitleFont(42);
	copy3->GetYaxis()->SetTitleOffset(1.3);
	copy3->GetXaxis()->SetTitle("Days since December 1st 2020");
	copy3->GetXaxis()->SetNdivisions(1020*4, "kTRUE");
	copy3->GetXaxis()->SetTickLength(0.);
 	copy3->GetXaxis()->SetLabelSize(.02);
	copy3->SetFillColor(1);
	copy3->Draw("BOX SAME");
	ap->GetZaxis()->SetTitle("Position of VEM-Peak / FADC");
  ap->GetZaxis()->SetTickLength(0.);
	ap->GetZaxis()->SetTitleOffset(1.3);
  ap->SetStats(0);
  ap->GetYaxis()->SetTickLength(0.);
  plain->SetPalette(87);
  c3->SetGridy();
  c3->SetGridx();
  ap->Draw("COLZ1 SAME");
  c3->Print("../plots/ap"+printname+"Hg.pdf");

}
