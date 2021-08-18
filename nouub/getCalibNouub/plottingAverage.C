void plottingAverage(int pmtId, int stId) {
	TString filename;
	TString printname;
	
	if ( pmtId > 0 && pmtId < 4 ) {
		filename.Form("ubCalibHistPMT%dSt%d.root", pmtId, stId);
		printname.Form("Pmt%dSt%d", pmtId, stId);
	}
	else if ( pmtId == 4 ){
		filename = "calibHistSPMT.root";
		printname = "SPMT";
	}
	else if ( pmtId == 5 ){
		filename = "calibHistPMTSSD.root";
		printname = "PMTSSD";
	}
	else {
		cout << "==================================================" << endl;
		cout << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
		cout << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " 
			<< "4 For SPMT; " << "5 For PMTSSD" << endl;
		cout << "==================================================" << endl;
		exit(0);
	}
	
	cout << "You have selected " << filename << endl;
	auto hFile = TFile::Open(filename);

	auto pmt = (TH1F*)hFile->Get("stckH");
	auto Ch = (TTree*)hFile->Get("forEntries");
	int ne = 0;
	Ch->SetBranchAddress("nentry", &ne);
  Ch->GetEntry(0);	
	pmt->Scale(1./ne);

	TCanvas *c1 = new TCanvas("c1", "Average for Charge", 1, 35,3600,2400);
	c1->cd();

	//pmt->GetXaxis()->SetRangeUser(1000, 2200);
	pmt->GetXaxis()->SetTitle("FADCq");
	pmt->GetYaxis()->SetTitle("Counts / au");
	pmt->SetMarkerSize(3);
	pmt->SetMarkerStyle(7);
	pmt->SetMarkerColor(4);
	pmt->Draw();
	c1->Print("../../plots/ubAveChHisto"+printname+".pdf");
}
