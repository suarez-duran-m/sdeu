void plottingBl100bins(int pmtId, int hg){  
  TString filename;
  TString printname;
  TString hgname;
  if ( pmtId > 0 && pmtId < 4 ){
    filename.Form("bl100binsPMT%d.root", pmtId);
    printname.Form("PMT%d", pmtId);
  }
  else if ( pmtId == 4 ){
    filename = "bl100binsSPMT.root";
    printname = "SPMT";
  }
  else if ( pmtId == 5 ){
    filename = "bl100binsPMTSSD.root";
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
  if ( hg !=0 && hg !=1 ){
    cout << hg << endl;
    cerr << "Please select HG (1) or LG (0) option" 
      << endl;
    exit(0);
  }
  cout << "You have selected " << filename << " "
    << "with HG option: " << hg << endl;
 
  auto hFile = TFile::Open(filename);

  TH2F *pmtmf;
  TH2F *pmtml;
  TH2F *pmtrf;
  TH2F *pmtrl;
  TH2F *pmtdm;
  TH2F *pmtdr;
  
  if ( hg==1 ){
    pmtmf = (TH2F*)hFile->Get("pmthmeanf");
    pmtml = (TH2F*)hFile->Get("pmthmeanl");

    pmtrf = (TH2F*)hFile->Get("pmthrmsf");
    pmtrl = (TH2F*)hFile->Get("pmthrmsl");

    pmtdm = (TH2F*)hFile->Get("pmthdiffmean");
    pmtdr = (TH2F*)hFile->Get("pmthdiffrms");

    hgname = "Hg";
  }
  else{
    pmtmf = (TH2F*)hFile->Get("pmtlmeanf");
    pmtml = (TH2F*)hFile->Get("pmtlmeanl");

    pmtrf = (TH2F*)hFile->Get("pmtlrmsf");
    pmtrl = (TH2F*)hFile->Get("pmtlrmsl");

    pmtdm = (TH2F*)hFile->Get("pmtldiffmean");
    pmtdr = (TH2F*)hFile->Get("pmtldiffrms");
    
    hgname = "Lg";
  }

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};

  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

	auto copy = (TH2F*)pmtdm->Clone();
	copy->Reset();
	copy->SetFillStyle(3004);
	Int_t bin = 0;
	Double_t maxVal = pmtdm->GetMaximum();
	
	for(Int_t i = 1; i <= pmtdm->GetNbinsX(); ++i){
		for(Int_t j = 1; j <= pmtdm->GetNbinsY(); ++j){
			bin = pmtdm->GetBin(i, j);
			if(!pmtdm->GetBinContent(bin))
				copy->SetBinContent(bin, maxVal);
		}
	}

  TCanvas c1("c5", "2D", 0,0,3600,2400);
  c1.cd();
  c1.SetLeftMargin(0.11);
	c1.SetRightMargin(0.135);

	copy->SetFillStyle(3004);
	copy->SetStats(0);
	copy->GetZaxis()->SetRangeUser(-2.,4.5);
	copy->GetZaxis()->SetTickLength(0.);
	copy->GetZaxis()->SetLabelSize(0.02);
	copy->SetTickLength(0.);
	copy->GetZaxis()->SetLabelSize(0.02);
	for ( int i=0; i<19; i++)
		copy->GetYaxis()->SetBinLabel(i+1, stIds[i]);
	copy->GetYaxis()->SetNdivisions(1020, "kTRUE");
	copy->GetYaxis()->SetTickLength(0.);
	copy->GetYaxis()->SetTitle("Station ID");
	copy->GetYaxis()->SetTitleFont(42);
	copy->GetYaxis()->SetTitleOffset(1.3);
	copy->GetYaxis()->SetLabelSize(.05);
	copy->SetFillColor(1);
 	copy->GetXaxis()->SetTitle("Days since December 1st 2020");
 	copy->GetXaxis()->SetNdivisions(1020*4, "kTRUE");
 	copy->GetXaxis()->SetRangeUser(0,90);
 	copy->GetXaxis()->SetTickLength(0.);
 	copy->GetXaxis()->SetLabelSize(.02);
	copy->Draw("BOX SAME");
  pmtdm->SetStats(0);
	pmtdm->GetZaxis()->SetTitle("Diff. Mean / FADC");
	pmtdm->GetZaxis()->SetTickLength(0.);

  plain->SetPalette(56);
  c1.SetGridy();
  c1.SetGridx();
  pmtdm->Draw("COLZ1 SAME");
  c1.Print("../../plots/bl"+printname+"Diffmean100"+hgname+".pdf");
}
