void readHisto(int pmt, int st) { 
	
	TString pmtId; 
	TString stat; 
	
	pmtId.Form("PMT%d", pmt); 
	stat.Form("St%d", st); 

	TFile f("uubCalibHist"+pmtId+stat+".root");
	cout << "Reading file: " << "uubCalibHist"+pmtId+stat+".root" <<endl;
	
	TTree *histos = (TTree*)f.Get("Histograms");

	TH1F *charge = new TH1F();
	TH1F *peak = new TH1F();
	TH1F *chHisto = new TH1F();
	TH1F *pkHisto = new TH1F();
	TH1F *baseline = new TH1F();
	TH1F *peakBl = new TH1F();
	TH1F *peakOff = new TH1F();
  TH1F *rawPeak = new TH1F();
	int offsetch;
	int offsetpk;
	int eventId;

  histos->SetBranchAddress("rawPk", &rawPeak);
	histos->SetBranchAddress("recePk", &peak);
	histos->SetBranchAddress("offSetCh", &offsetch);
	histos->SetBranchAddress("receCh", &charge);
 /* 
	histos->SetBranchAddress("setChHisto", &chHisto);
	histos->SetBranchAddress("setPkHisto", &pkHisto);
	histos->SetBranchAddress("base", &baseline);
	histos->SetBranchAddress("entryEvt", &eventId);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("pkCorrBl", &peakBl);
	histos->SetBranchAddress("pkCorrOff", &peakOff);
  */

	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);

	cout << histos->GetEntries() << endl;
	histos->GetEntry(0);
	cout << eventId << endl;

	TString nameEvt;
	nameEvt.Form("%d", eventId);

	TString nameOffset;
	nameOffset.Form("%d", offsetpk);

	// ==========================================
	// *************** For Peak *****************
	
	c1->cd();
  histos->GetEntry(0);
  TLegend *legend;
  /*
	rawPeak->SetYTitle("Counts [au]");
	rawPeak->GetYaxis()->SetTitleSize(0.05);
	rawPeak->GetYaxis()->SetTitleOffset(1.);
	rawPeak->GetYaxis()->SetLabelSize(0.04);
	rawPeak->SetXTitle("Bin [au]");
	rawPeak->GetXaxis()->SetTitleSize(0.045);
	rawPeak->GetXaxis()->SetLabelSize(0.04); // >SetTitleSize(1.0);
	rawPeak->SetStats(kFALSE);
	rawPeak->SetLineColor(kBlue);
	rawPeak->Draw();
  legend = new TLegend(0.45,0.77,0.89,0.87);
  //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  legend->AddEntry(rawPeak,"#splitline{UUB Peak-Raw-Histogram PMT1}{St863 event: 61219267}","f");
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.035);
  //legend->AddEntry("f1","Function abs(#frac{sin(x)}{x})","l");
  //legend->AddEntry("gr","Graph with error bars","lep");
  legend->Draw();
	c1->Print("../../plots/uubRawPeak"+pmtId+stat+".pdf");
  */
/*
	c1->cd();
	peak->SetTitle("");
      //"UUB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Offset: "+nameOffset);
	peak->SetYTitle("Counts [au]");
	peak->GetYaxis()->SetTitleSize(0.05);
	peak->GetYaxis()->SetTitleOffset(1.);
	peak->GetYaxis()->SetLabelSize(0.04);
	peak->SetXTitle("Bin [FADC]");
	peak->GetXaxis()->SetTitleSize(0.045);
	peak->GetXaxis()->SetLabelSize(0.04);
	peak->GetXaxis()->SetRangeUser(0, 150);
	peak->SetStats(kFALSE);
	peak->SetLineColor(kBlue);
	peak->Draw();
  legend = new TLegend(0.5,0.77,0.88,0.87);
  legend->AddEntry(rawPeak,"#splitline{UUB Peak-Histogram PMT1}{St863 event: 61219267}","f");
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.035);
  legend->Draw();
	c1->Print("../../plots/uubPeak"+pmtId+stat+".pdf");
	//cout << "First bin for Peak histogram: " << peak->GetBinLowEdge(1) << endl;
*/


  /*
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
	c1->Print("../plots/uubBase"+pmtId+stat+".pdf");
*/
/*
	c1->cd();
	peakBl->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Corr. BL: 284");
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
	//c1->Print("../plots/uubPeakCoorBl"+pmtId+stat+".pdf");
  */

/*
	c1->cd();
	peakOff->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Corr. Offset: 273");
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
	//c1->Print("../plots/uubPeakCoorOff"+pmtId+stat+".pdf");
  */

  /*
	c1->cd();
	peakOff->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Zoom Corr. Offset: 273");
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
	c1->Print("../plots/uubPeakCoorOffZoom"+pmtId+stat+".pdf");




	// ==========================================
	// ************** For Charge ****************

	c1->cd();
	chHisto->SetTitle("UUB Charge-Raw-Histogram "+pmtId+" "+stat+" event: "+nameEvt);
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
	c1->Print("../plots/uubRawCharge"+pmtId+stat+".pdf");
*/

  histos->GetEntry(0);
	nameOffset.Form("%d", offsetch);

	c1->cd();
	charge->SetTitle(""); //UUB Charge-Histogram "+pmtId+" "+stat+" event: "+nameEvt+" Offset: "+nameOffset);
	charge->SetYTitle("Counts [au]");
	charge->GetYaxis()->SetTitleSize(0.05);
	charge->GetYaxis()->SetTitleOffset(1.);
	charge->GetYaxis()->SetLabelSize(0.04);
	charge->SetXTitle("Bin [FADCq]");
	charge->GetXaxis()->SetTitleSize(0.045);
	charge->GetXaxis()->SetLabelSize(0.04);
	charge->SetStats(kFALSE);
	charge->SetLineColor(kBlue);
	charge->Draw();
  legend = new TLegend(0.5,0.77,0.88,0.87);
  legend->AddEntry(charge,"#splitline{UUB Charge-Histogram PMT1}{St863 event: 61219267}","f");
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.035);
  legend->Draw();
	c1->Print("../../plots/uubCharge"+pmtId+stat+".pdf");
	//cout << "First bin for Charge histogram: " << charge->GetBinLowEdge(1) << endl;


	// ==================================
	// *** OffSet as funciton of time ***
/*
	Int_t nx = histos->GetEntries();
	Double_t x[nx];
	Double_t y[nx];

	// ================
	// *** For Peak ***

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk;	
	}

	TGraph* gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset for Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetPk"+pmtId+stat+".pdf");

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk - peak->GetBinLowEdge(1);
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset minus first-bin Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetDiffPk"+pmtId+stat+".pdf");

	// ==================
	// *** For Charge ***

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetch;
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset for Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	if ( pmtId == "PMT1" || pmtId == "PMT2")
		gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetCh"+pmtId+stat+".pdf");

	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		x[i] = i;
		y[i] = offsetpk - peak->GetBinLowEdge(1);
	}

	gr = new TGraph(nx,x,y);
	c1->cd();
	gr->SetTitle("Offset minus first-bin Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1, 1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetDiffCh"+pmtId+stat+".pdf");

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
	c1->Print("../plots/uubAveChHisto"+pmtId+stat+".pdf");


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
	c1->Print("../plots/uubAvePkHisto"+pmtId+stat+".pdf");
*/

}
