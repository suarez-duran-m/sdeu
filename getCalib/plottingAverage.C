void plottingAverage(int pmtId, int stId) {
	TString filename;
	TString printname;
	
	if ( pmtId > 0 && pmtId < 4 ) {
		filename.Form("uubCalibHistPMT%dSt%d.root", pmtId, stId);
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

	auto stckCh = (TH1F*)hFile->Get("stckHc");
	auto stckPk = (TH1F*)hFile->Get("stckHp");
	auto Ch = (TTree*)hFile->Get("forEntries");
	int ne = 0;
	Ch->SetBranchAddress("nentry", &ne);
	Ch->GetEntry(0);
	stckCh->Scale(1./ne);
	stckPk->Scale(1./ne);


	TCanvas *c1 = new TCanvas("c1", "Average for Charge", 1, 35,3600,2400);
	c1->cd();

	stckCh->GetXaxis()->SetRangeUser(0, 5000);
	stckCh->GetXaxis()->SetTitle("FADCq");
	stckCh->GetYaxis()->SetTitle("Counts / au");
	stckCh->SetMarkerSize(3);
	stckCh->SetMarkerStyle(7);
	stckCh->SetMarkerColor(4);
	stckCh->Draw("P");
	c1->Print("../plots/uubAveChHisto"+printname+".pdf");


	TCanvas *c2 = new TCanvas("c2", "Average for Peak", 1, 35,3600,2400);
	c2->cd();

	stckPk->GetXaxis()->SetRangeUser(0, 2000);
	stckPk->GetXaxis()->SetTitle("FADCq");
	stckPk->GetYaxis()->SetTitle("Counts / au");
	stckPk->SetMarkerSize(3);
	stckPk->SetMarkerStyle(7);
	stckPk->SetMarkerColor(4);
	stckPk->Draw("P");
	c2->Print("../plots/uubAvePkHisto"+printname+".pdf");

}
