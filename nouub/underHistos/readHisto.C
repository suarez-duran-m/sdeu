void readHisto(int pmt, int st) { 
	
	TString pmtId; 
	TString stat; 
	
	pmtId.Form("PMT%d", pmt); 
	stat.Form("St%d", st); 

	TFile f("ubCalibHist"+pmtId+stat+".root");
	cout << "Reading file: " << "ubCalibHist"+pmtId+stat+".root" <<endl;
	
	TTree *histos = (TTree*)f.Get("Histograms");

	TH1F *charge = new TH1F();
	TH1F *peak = new TH1F();
	TH1F *chHisto = new TH1F();
	TH1F *pkHisto = new TH1F();
	TH1F *baseline = new TH1F();
	TH1F *peakBl = new TH1F();
	TH1F *peakOff = new TH1F();
	int offsetch;
	int offsetpk;
	int fstPkraw;
	int fstChraw;
	int eventId;

	histos->SetBranchAddress("setCh", &charge);

	histos->SetBranchAddress("setPk", &peak);
	histos->SetBranchAddress("setChHisto", &chHisto);
	histos->SetBranchAddress("setPkHisto", &pkHisto);
	histos->SetBranchAddress("base", &baseline);
	histos->SetBranchAddress("entryEvt", &eventId);
	histos->SetBranchAddress("offSetCh", &offsetch);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("fstRawBinPk", &fstPkraw);
	histos->SetBranchAddress("fstRawBinCh", &fstChraw);
	histos->SetBranchAddress("pkCorrBl", &peakBl);
	histos->SetBranchAddress("pkCorrOff", &peakOff);

	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);
/*
	cout << histos->GetEntries() << endl;
	histos->GetEntry(0);
	cout << eventId << endl;

	TString nameEvt;
	nameEvt.Form("%d", eventId);

	// ==========================================
	// *************** For Peak *****************
	
	c1->cd();
	pkHisto->SetTitle("UB Peak-Raw-Histogram "+pmtId+" "+stat+" event: "+nameEvt);
	pkHisto->SetYTitle("Counts / au");
	pkHisto->GetYaxis()->SetTitleSize(0.05);
	pkHisto->GetYaxis()->SetTitleOffset(1.1);
	pkHisto->GetYaxis()->SetLabelSize(0.04);
	pkHisto->SetXTitle("1 / ua");
	pkHisto->GetXaxis()->SetTitleSize(0.045);
	pkHisto->GetXaxis()->SetLabelSize(0.04); // >SetTitleSize(1.0);
	pkHisto->GetXaxis()->SetRangeUser(0,150);
	pkHisto->SetStats(kFALSE);
	pkHisto->SetLineColor(kBlue);
	pkHisto->Draw();
	c1->Print("../../plots/ubRawPeak"+pmtId+stat+".pdf");

	c1->cd();
	pkHisto->SetTitle("UB Peak-Raw-Histogram "+pmtId+" "+stat+" event: "+nameEvt);
	pkHisto->SetYTitle("Counts / au");
	pkHisto->GetYaxis()->SetTitleSize(0.05);
	pkHisto->GetYaxis()->SetTitleOffset(1.1);
	pkHisto->GetYaxis()->SetLabelSize(0.04);
	pkHisto->GetYaxis()->SetRangeUser(0, 3800);
	pkHisto->SetXTitle("1 / ua");
	pkHisto->GetXaxis()->SetTitleSize(0.045);
	pkHisto->GetXaxis()->SetLabelSize(0.04); // >SetTitleSize(1.0);
	pkHisto->GetXaxis()->SetRangeUser(0,20);
	pkHisto->SetStats(kFALSE);
	pkHisto->SetLineColor(kBlue);
	pkHisto->Draw();
	c1->Print("../../plots/ubRawPeakZoom"+pmtId+stat+".pdf");


	TString nameOffset;
	nameOffset.Form("%d", offsetpk);

	c1->cd();
	peak->SetTitle("UB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Offset: "+nameOffset);
	peak->SetYTitle("Counts / au");
	peak->GetYaxis()->SetTitleSize(0.05);
	peak->GetYaxis()->SetTitleOffset(1.1);
	peak->GetYaxis()->SetLabelSize(0.04);
	peak->SetXTitle("1 / FADCp");
	peak->GetXaxis()->SetTitleSize(0.045);
	peak->GetXaxis()->SetLabelSize(0.04);
	//peak->GetXaxis()->SetRangeUser(0, 150);
	peak->SetStats(kFALSE);
	peak->SetLineColor(kBlue);
	peak->Draw();
	c1->Print("../../plots/ubPeak"+pmtId+stat+".pdf");
	cout << "First bin for Peak histogram: " << peak->GetBinLowEdge(1) << endl;

	c1->cd();
	baseline->SetTitle("Baseline "+pmtId+" "+stat+" event: "+nameEvt);
	baseline->SetYTitle("Counts / au");
	baseline->GetYaxis()->SetTitleSize(0.05);
	baseline->GetYaxis()->SetTitleOffset(1.1);
	baseline->GetYaxis()->SetLabelSize(0.04);
	baseline->SetXTitle("1 / FADC");
	baseline->GetXaxis()->SetTitleSize(0.045);
	baseline->GetXaxis()->SetLabelSize(0.04);
	baseline->SetStats(kTRUE);
	baseline->SetLineColor(kBlue);
	baseline->SetLineWidth(2.);
	baseline->Draw();
	c1->Print("../../plots/ubBase"+pmtId+stat+".pdf");


	c1->cd();
	peakBl->SetTitle("UB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Corr. BL: 58");
	peakBl->SetYTitle("Counts / au");
	peakBl->GetYaxis()->SetTitleSize(0.05);
	peakBl->GetYaxis()->SetTitleOffset(1.1);
	peakBl->GetYaxis()->SetLabelSize(0.04);
	peakBl->SetXTitle("1 / FADCp");
	peakBl->GetXaxis()->SetTitleSize(0.045);
	peakBl->GetXaxis()->SetLabelSize(0.04);
	peakBl->SetStats(kFALSE);
	peakBl->SetLineColor(kBlue);
	peakBl->Draw();
	c1->Print("../../plots/ubPeakCoorBl"+pmtId+stat+".pdf");


	c1->cd();
	peakOff->SetTitle("UB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Corr. Offset: 57");
	peakOff->SetYTitle("Counts / au");
	peakOff->GetYaxis()->SetTitleSize(0.05);
	peakOff->GetYaxis()->SetTitleOffset(1.1);
	peakOff->GetYaxis()->SetLabelSize(0.04);
	peakOff->SetXTitle("1 / FADCp");
	peakOff->GetXaxis()->SetTitleSize(0.045);
	peakOff->GetXaxis()->SetLabelSize(0.04);
	peakOff->SetStats(kFALSE);
	peakOff->SetLineColor(kBlue);
	peakOff->Draw();
	c1->Print("../../plots/ubPeakCoorOff"+pmtId+stat+".pdf");

	c1->cd();
	peakOff->SetTitle("UB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Zoom Corr. Offset: 57");
	peakOff->SetYTitle("Counts / au");
	peakOff->GetYaxis()->SetTitleSize(0.05);
	peakOff->GetYaxis()->SetTitleOffset(1.1);
	peakOff->GetYaxis()->SetLabelSize(0.04);
	peakOff->SetXTitle("1 / FADCp");
	peakOff->GetXaxis()->SetTitleSize(0.045);
	peakOff->GetXaxis()->SetLabelSize(0.04);
	peakOff->GetXaxis()->SetRangeUser(0,60);
	peakOff->SetStats(kFALSE);
	peakOff->SetLineColor(kBlue);
	peakOff->Draw();
	c1->Print("../../plots/ubPeakCoorOffZoom"+pmtId+stat+".pdf");


	// ==========================================
	// ************** For Charge ****************

	c1->cd();
	chHisto->SetTitle("UB Charge-Raw-Histogram "+pmtId+" "+stat+" event: "+nameEvt);
	chHisto->SetYTitle("Counts / au");
	chHisto->GetYaxis()->SetTitleSize(0.05);
	chHisto->GetYaxis()->SetTitleOffset(1.1);
	chHisto->GetYaxis()->SetLabelSize(0.04);
	chHisto->SetXTitle("1 / au");
	chHisto->GetXaxis()->SetTitleSize(0.045);
	chHisto->GetXaxis()->SetLabelSize(0.04); // >SetTitleSize(1.0);
	chHisto->SetStats(kFALSE);
	chHisto->SetLineColor(kBlue);
	chHisto->Draw();
	c1->Print("../../plots/ubRawCharge"+pmtId+stat+".pdf");

	nameOffset.Form("%d", offsetch);

	c1->cd();
	charge->SetTitle("UB Charge-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Offset: "+nameOffset);
	charge->SetYTitle("Counts / au");
	charge->GetYaxis()->SetTitleSize(0.05);
	charge->GetYaxis()->SetTitleOffset(1.1);
	charge->GetYaxis()->SetLabelSize(0.04);
	charge->SetXTitle("1 / FADCq");
	charge->GetXaxis()->SetTitleSize(0.045);
	charge->GetXaxis()->SetLabelSize(0.04);
	//charge->GetXaxis()->SetRangeUser(0, 150);
	charge->SetStats(kFALSE);
	charge->SetLineColor(kBlue);
	charge->Draw();
	c1->Print("../../plots/ubCharge"+pmtId+stat+".pdf");
	cout << "First bin for Charge histogram: " << charge->GetBinLowEdge(1) << endl;

*/
	// ==================================
	// *** OffSet as funciton of time ***

	Int_t nx = histos->GetEntries();
	Double_t x[nx];
	Double_t y[nx];
	TGraph* gr;

	/*
	// ================
	// *** For Peak ***

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk;	
	}

	TGraph* gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset for Peak Histograms UB "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetPk"+pmtId+stat+".pdf");

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk - peak->GetBinLowEdge(1);
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset minus first-bin Peak Histograms UB "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetDiffPk"+pmtId+stat+".pdf");

	// ==================
	// *** For Charge ***

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetch;
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset for Charge Histograms UB "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	//if ( pmtId == "PMT1" || pmtId == "PMT2")
		//gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetCh"+pmtId+stat+".pdf");

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk - peak->GetBinLowEdge(1);
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset minus first-bin Charge Histograms UB "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetDiffCh"+pmtId+stat+".pdf");
*/

	// =============================================
	// *** Counts in first bin of raw histograms ***

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = eventId;
		y[i] = fstPkraw;
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Counts in first bin for UB Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetXaxis()->SetTitleOffset(1.2);
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.5);
	gr->GetYaxis()->SetRangeUser(-5e3,6e4);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubCntFirstBinPk"+pmtId+stat+".pdf");

/*
	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = fstChraw;
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Counts in first bin for UB Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubCntFirstBinCh"+pmtId+stat+".pdf");
*/


/*
	TCanvas *c2 = new TCanvas("c2","C2",10,10,600,400); 
	c2->cd();

	for(Int_t evt=0; evt<histos->GetEntries(); evt++) { 
		cout << evt << endl; 
		//c2->Clear();
		//c2->cd();
		histos->GetEntry(evt);	
		cout << "Offset: " << offset << endl;
		//baseline->SetStats();
		//baseline->Draw();
		charge->Draw();
		gPad->WaitPrimitive();
	}
*/


/*
	histos->GetEntry(2);
	Int_t nx = 1.e5;
	Double_t x[nx];
	Double_t y[nx];
	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		x[tmp] = 0.;
		y[tmp] = 0.;
	}

	for (Int_t i=0; i<histos->GetEntries(); i++) {
		histos->GetEntry(i);
		for (Int_t b=1; b<charge->GetXaxis()->GetNbins(); b++) {
			x[ (int)charge->GetBinCenter(b) ] = charge->GetBinCenter(b);
			y[ (int)charge->GetBinCenter(b) ] += charge->GetBinContent(b);
		}
	}

	for ( Int_t tmp=0; tmp<nx; tmp++ )
		y[tmp] /= histos->GetEntries();

	TGraph* gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Average Charge Histogram");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("FADCq");
	gr->GetYaxis()->SetTitle("Counts / au");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/ubAveChHisto"+pmtId+stat+".pdf");


	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		x[tmp] = 0.;
		y[tmp] = 0.;
	}

	for (Int_t i=0; i<histos->GetEntries(); i++) {
		histos->GetEntry(i);
		for (Int_t b=0; b<peak->GetXaxis()->GetNbins(); b++) {
			x[ (int)peak->GetBinCenter(b) ] = peak->GetBinCenter(b);
			y[ (int)peak->GetBinCenter(b) ] += peak->GetBinContent(b);
		}
	}

	for ( Int_t tmp=0; tmp<nx; tmp++ )
		y[tmp] /= histos->GetEntries();

	TGraph* gr1 = new TGraph(nx,x,y);
	c1->cd();
	gr1->SetTitle("Average Peak Histogram");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr1->GetXaxis()->SetTitle("FADCp");
	gr1->GetYaxis()->SetTitle("Counts / au");
	gr1->GetYaxis()->SetTimeOffset(1.3);
	gr1->SetFillStyle(1000);
	gr1->SetMarkerStyle(20);
	gr1->SetMarkerSize(0.4);
	gr1->SetMarkerColor(4);
	gr1->Draw("AP");
	c1->Print("../plots/ubAvePkHisto"+pmtId+stat+".pdf");
*/

}
