void plottingAverage(int pmtId) {
	TString filename;
	TString printname;
	
	if ( pmtId > 0 && pmtId < 4 ) {
		filename.Form("uubCalibHistPMT%d.root", pmtId);
		printname.Form("Pmt%dSt%d", pmtId);
	}
	else if ( pmtId == 4 ){
		filename = "uubCalibHistSPMT.root";
		printname = "SPMT";
	}
	else if ( pmtId == 5 ){
		filename = "uubCalibHistPMTSSD.root";
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

	auto distAp = (TH1F*)hFile->Get("apDist");

	TCanvas *c1 = new TCanvas("c1", "A/P RMS", 1, 35,3600,2400);
	c1->cd();

	//stckCh->GetXaxis()->SetRangeUser(0, 5000);
	//stckCh->GetXaxis()->SetTitle("FADCq");
	stckCh->GetYaxis()->SetTitle("Counts / au");
	stckCh->SetMarkerSize(3);
	stckCh->SetMarkerStyle(7);
	stckCh->SetMarkerColor(4);
	stckCh->Draw("P");
	c1->Print("../plots/uubAveChHisto"+printname+".pdf");
}
