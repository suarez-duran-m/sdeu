void plotting(int pmt) { 
	
	TString pmtId; 
	
	pmtId.Form("PMT%d", pmt); 

	//TFile f("uubCalibHistAll"+pmtId+".root");
	auto f = TFile::Open("uubCalibHistAll"+pmtId+".root");
	cout << "Reading file: " << "uubCalibHistAll"+pmtId+".root" <<endl;
	
	TTree *histos = (TTree*)f->Get("Histograms");

	TH1F *charge = new TH1F();
	TH1F *peak = new TH1F();
	TH1F *chHisto = new TH1F();
	TH1F *pkHisto = new TH1F();
	TH1F *baseline = new TH1F();
	int offsetch;
	int offsetpk;
	int eventId;

	histos->SetBranchAddress("setCh", &charge);

	histos->SetBranchAddress("setPk", &peak);
	histos->SetBranchAddress("setChHisto", &chHisto);
	histos->SetBranchAddress("setPkHisto", &pkHisto);
	histos->SetBranchAddress("base", &baseline);
	histos->SetBranchAddress("entryEvt", &eventId);
	histos->SetBranchAddress("offSetCh", &offsetch);
	histos->SetBranchAddress("offSetPk", &offsetpk);

	TH1F *diffOffsetP = (TH1F*)f->Get("diffOffsetStationP");
	TH1F *diffOffsetC = (TH1F*)f->Get("diffOffsetStationC");

	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};


	// ============================================
	// *** Checking Offset differences for Peak ***
	
	c1->cd();

	diffOffsetP->SetYTitle("Difference / au");
	diffOffsetP->GetYaxis()->SetTitleSize(0.05);
	diffOffsetP->GetYaxis()->SetTitleOffset(1.);
	diffOffsetP->GetYaxis()->SetLabelSize(0.04);
	diffOffsetP->GetYaxis()->SetRangeUser(-1,1);
	diffOffsetP->SetXTitle("Station Id");
	diffOffsetP->GetXaxis()->SetTitleSize(0.045);
	diffOffsetP->GetXaxis()->SetLabelSize(0.04);
	for ( int i=0; i<19; i++)
		diffOffsetP->GetXaxis()->SetBinLabel(i+1, stIds[i]);
	diffOffsetP->SetLineColor(kBlue);
	diffOffsetP->SetLineWidth(2.0);
	diffOffsetP->Draw();
	c1->Print("../../plots/uubDiffOffsetPeak"+pmtId+".pdf");

	// ============================================
	// *** Checking Offset differences for Peak ***

	c1->cd();

	diffOffsetC->SetYTitle("Difference / au");
	diffOffsetC->GetYaxis()->SetTitleSize(0.05);
	diffOffsetC->GetYaxis()->SetTitleOffset(1.);
	diffOffsetC->GetYaxis()->SetLabelSize(0.04);
	diffOffsetC->GetYaxis()->SetRangeUser(-1,1);
	diffOffsetC->SetXTitle("Station Id");
	diffOffsetC->GetXaxis()->SetTitleSize(0.045);
	diffOffsetC->GetXaxis()->SetLabelSize(0.04);
	for ( int i=0; i<19; i++)
		diffOffsetC->GetXaxis()->SetBinLabel(i+1, stIds[i]);
	diffOffsetC->SetLineColor(kBlue);
	diffOffsetC->SetLineWidth(2.0);
	diffOffsetC->Draw();
	c1->Print("../../plots/uubDiffOffsetCharge"+pmtId+".pdf");

}
