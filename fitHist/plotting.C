void plotting(int pmt, int st) {

	TString pmtId;
	TString stat;

	pmtId.Form("PMT%d", pmt);
	stat.Form("St%d", st);

	TFile *f = TFile::Open("uubCalibHist"+pmtId+stat+".root");
	//TFile f ("uubCalibHist"+pmtId+stat+".root");
	cout << "Reading file: " << "uubCalibHist"+pmtId+stat+".root" <<endl;

 	TTree *histos = (TTree*)f->Get("Histograms"); 
	TH1F *hcharge = new TH1F();
	TH1F *hpeak = new TH1F();
	TH1F *peak = new TH1F(); 
	TH1F *charge = new TH1F();
	TH1F *peakCalibBl = new TH1F();
	TH1F *chargeCalibBl = new TH1F();
 	TGraphErrors *fittedPk = new TGraphErrors(); 
 	TGraphErrors *fittedCh = new TGraphErrors(); 

	double chpk;
	double blHbase; 
	double blCalib;
	int offsetpk;
	int offsetch;
	int evtTime;
 
	histos->SetBranchAddress("receCh", &hcharge);	
	histos->SetBranchAddress("recePk", &hpeak);	
	histos->SetBranchAddress("corrPk",&peak); 	
	histos->SetBranchAddress("corrCh",&charge);
	histos->SetBranchAddress("fllHPk", &fittedPk); 
	histos->SetBranchAddress("fllHCh", &fittedCh);
	histos->SetBranchAddress("baselineHbase", &blHbase);
	histos->SetBranchAddress("baselineCalib", &blCalib);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("offSetCh", &offsetch);
	histos->SetBranchAddress("ap", &chpk);
	histos->SetBranchAddress("entryEvt", &evtTime);

	TCanvas *c0 = new TCanvas("c0","C0",500,10,1200,800);
	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);

	int nx = histos->GetEntries();
	TGraph* gr;
	double xbl[nx];
	double ybl[nx];
	double minlimHb = 1e4;
	double maxlimHb = 0.;
	double minlimCl = 1e4;
	double maxlimCl = 0.;

	TString p0nm;
	TString p1nm;
	TPaveStats *ptstats;
	TText *ptstats_LaTex;

	
	histos->GetEntry(0);
	c0->cd();
	fittedCh->Draw();
	cerr << offsetch << endl;
	//fittedPk->Draw();
//	c1->cd();
	//hcharge->Draw();
	//fittedCh->Draw();


	// ==================
	// *** For Offset ***

	/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = offsetpk;
	}


	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Offset for Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetPk"+pmtId+stat+".pdf");
*/

/*	
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = offsetch;
	}

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Offset for Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	if ( pmt==3 )
		gr->GetYaxis()->SetTitleOffset(1.5);
	else
		gr->GetYaxis()->SetTitleOffset(1.3);
	if ( pmt!=3 )
		gr->GetYaxis()->SetRangeUser(-1,1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubOffsetCh"+pmtId+stat+".pdf");
*/



	// ====================
	// *** For Baseline ***

/*
	TH1F *histBlHbase = new TH1F("histBlHbase", "UUB Mean Baseline HBase "+pmtId+" "+stat, 100, 200, 300);
	TH1F *histBlCalib = new TH1F("histBlCalib", "UUB Baseline Calib "+pmtId+" "+stat, 100, 200, 300);

	TF1 *linef = new TF1("linef", "[1]*x + [0]");

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = blHbase;
		histBlHbase->Fill( blHbase );
		if ( blHbase < minlimHb )
			minlimHb = blHbase;
		if ( blHbase > maxlimHb )
			maxlimHb = blHbase;
	}
	minlimHb -= 1.;
	maxlimHb += 2.;

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Mean Baseline HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("Mean Baseline HBase / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	p0nm.Form("%.2f", linef->GetParameters()[0]);
	p1nm.Form("%.2e", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.74,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);	
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c1->Print("../plots/uubBlHbase"+pmtId+stat+".pdf");


	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );		
		xbl[tmp] = evtTime;
		ybl[tmp] = blCalib;
		histBlCalib->Fill( blCalib );
		if ( blCalib < minlimCl )
			minlimCl = blCalib;
		if ( blCalib > maxlimCl )
			maxlimCl = blCalib;
	}

	minlimCl -= 1.;
	maxlimCl += 2.;

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Baseline Calib "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("Calib.Base / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	p0nm.Form("%.2f", linef->GetParameters()[0]);
	p1nm.Form("%.2e", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.74,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);	
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c1->Print("../plots/uubBlCalib"+pmtId+stat+".pdf");

	// ===========================
	// *** For Baseline histos ***

	c0->cd();
	histBlHbase->SetLineWidth(2.0);
	histBlHbase->GetXaxis()->SetRangeUser((int)minlimHb, (int)maxlimHb);	
	histBlHbase->SetStats(kTRUE);
	histBlHbase->GetYaxis()->SetTitle("Counts / au");
	histBlHbase->GetXaxis()->SetTitle("1 / FADC");
	histBlHbase->Draw();
	c0->Print("../plots/uubBlHistHbase"+pmtId+stat+".pdf");

	c0->cd();
	histBlCalib->SetLineWidth(2.0);
	histBlCalib->GetXaxis()->SetRangeUser((int)minlimCl, (int)maxlimCl);
	histBlCalib->SetStats(kTRUE);
	histBlCalib->GetYaxis()->SetTitle("Counts / au");
	histBlCalib->GetXaxis()->SetTitle("1 / FADC");
	histBlCalib->Draw();
	c0->Print("../plots/uubBlHistCalib"+pmtId+stat+".pdf");
*/


	// =================================
	// *** For Corrected Peak Histos ***
/*
	histos->GetEntry(0); // Set for first event
	int cnt = 0;

	for (int i=0; i<peak->GetXaxis()->GetNbins(); i++)
		cnt += peak->GetBinContent(i);

	p0nm.Form("%.2d",cnt);

	cerr << blHbase << " " << blCalib << endl;
	cerr << peak->GetBinCenter(1) << endl;
	*/

/*
	c1->cd();
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/uubPeakCoorHBase"+pmtId+stat+".pdf");
*/
/*
	c1->cd();	
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetXaxis()->SetRangeUser(-60, 60);
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/uubPeakCoorHBaseZoom"+pmtId+stat+".pdf");
*/

	// =====================================
	// *** First Bin for Peak Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = peak->GetBinCenter(1);
	}
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubPeakFirstBinHBase"+pmtId+stat+".pdf");
*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubPeakFirstBinCalibBase"+pmtId+stat+".pdf");

*/
















	// =============================
	// *** For Charge Histograms ***
/*
	histos->GetEntry(0);
	c1->cd();
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/uubChargeCoorHBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 283");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/uubChargeCoorHBaseZoom"+pmtId+stat+".pdf");
*/

	/*
	c1->cd();
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/uubChargeCoorCalibBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UUB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../plots/uubChargeCoorCalibBaseZoom"+pmtId+stat+".pdf");
*/


	// =====================================
	// *** First Bin for Peak Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = peak->GetBinCenter(1);
	}
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Peak histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubPeakFirstBinHBase"+pmtId+stat+".pdf");
*/


	// =======================================
	// *** First Bin for Charge Histograms ***
/*
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = charge->GetBinCenter(1);
	}
	*/

/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Charge histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubChargeFirstBinHBase"+pmtId+stat+".pdf");
	*/
/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB Charge histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubChargeFirstBinCalibBase"+pmtId+stat+".pdf");
	*/


	// =============================
	// *** Area/Peak Calculation ***

	nx = histos->GetEntries();
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		xbl[tmp] = evtTime;
		ybl[tmp] = chpk;
	}
	

	/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" HBase-Offset");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / 8.33*ns");
	gr->GetYaxis()->SetRangeUser(0, 12);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubApOffset"+pmtId+stat+".pdf");
	*/

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" Calib.Base-Offset");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / 8.33*ns");
	gr->GetYaxis()->SetRangeUser(0, 12);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../plots/uubApOffsetCalib"+pmtId+stat+".pdf");


/*
	c1->cd();
	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/uubPeakCoorCalibBase"+pmtId+stat+".pdf");

	peak->SetTitle("UUB Peak-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 233");
	peak->SetStats(0);
	peak->SetLineColor(kBlue);
	peak->GetXaxis()->SetTitle("1. / FADC");
	peak->GetXaxis()->SetRangeUser(40, 60);
	peak->GetYaxis()->SetTitle("Counts / au");
	peak->GetYaxis()->SetTitleOffset(1.3);
	peak->Draw();
	c1->Print("../plots/uubPeakCoorCalibBaseZoom"+pmtId+stat+".pdf");
*/


/*	
	TLatex latex(9e2,3.4e3,"Entries "+p0nm);
	latex.SetTextSize(0.04);
	latex.DrawClone("Same");
*/

	
/*
	int nx = histos->GetEntries();
	double xbl[nx];
	double ybl[nx];

	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) {
		histos->GetEntry(tmp);
		xbl[tmp] = tmp;
		ybl[tmp] = chpk;
	}

	TGraph* gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB First bin for Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d:%H");
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
*/

	/*
	// ================================
	// *** For Baseline Differences ***

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = tmp;
		ybl[tmp] = blHbase - blCalib;
	}

	TF1 *linef = new TF1("linef", "[1]*x + [0]");

	gr = new TGraph(nx,xbl,ybl);
	c0->cd();
	gr->SetTitle("UUB Baseline Difference HBase - Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("(HBase - Calib.Base) / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Fit(linef, "Q");

	TString p0nm;
	TString p1nm;
	p0nm.Form("%f", linef->GetParameters()[0]);
	p1nm.Form("%f", linef->GetParameters()[1]);
	TPaveStats *ptstats = new TPaveStats(0.62,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  TText *ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());
	gr->Draw("AP");
	c0->Print("../plots/uubDiffBlHbaseCalib"+pmtId+stat+".pdf");


	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xbl[tmp] = tmp;
		ybl[tmp] = offsetpk - blCalib;
	}

	gr = new TGraph(nx,xbl,ybl);
	c0->cd();
	gr->SetTitle("UUB Offset-Peak - Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("(OffsetPeak - Calib.Base) / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);

	gr->Fit(linef, "Q");
	p0nm.Form("%f", linef->GetParameters()[0]);
	p1nm.Form("%f", linef->GetParameters()[1]);
	ptstats = new TPaveStats(0.62,0.775,0.98,0.935,"brNDC");
  ptstats->SetName("stats");
  ptstats->SetBorderSize(1);
  ptstats->SetFillColor(0);
  ptstats->SetTextAlign(12);
  ptstats->SetTextFont(42);
  ptstats_LaTex = ptstats->AddText("p1*x + p0");
  ptstats_LaTex = ptstats->AddText("p0 = "+p0nm);
  ptstats_LaTex = ptstats->AddText("p1 = "+p1nm);
  ptstats->SetOptStat(0);
  ptstats->Draw();
  gr->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(gr->GetListOfFunctions());

	gr->Draw("AP");
	c0->Print("../plots/uubDiffBlOffsetPkCalib"+pmtId+stat+".pdf");
*/

	// =================================
	// *** Plotting Fitted histograms ***


	//histos->GetEntry(0);
/*
	c1->cd();
	peakCalibBl->SetTitle("UUB Peak histogram "+pmtId+" "+stat+" Corr. for Calib.Base") ;
	peakCalibBl->SetStats(kFALSE);
	peakCalibBl->GetXaxis()->SetTitle("1 / FADCp");
	peakCalibBl->GetYaxis()->SetTitle("Counts / au");
	peakCalibBl->GetYaxis()->SetTitleOffset(1.3);
	peakCalibBl->SetLineColor(kBlue);
	peakCalibBl->Draw();
	c1->Print("../plots/uubPeakCorrCalib"+pmtId+stat+".pdf");
	*/
/*
	c1->cd();
	peakCalibBl->SetTitle("UUB Peak histogram "+pmtId+" "+stat+" Corr. for Calib.Base") ;
	peakCalibBl->SetStats(kFALSE);
	peakCalibBl->GetXaxis()->SetRange(-50,40);
	peakCalibBl->GetXaxis()->SetTitle("1 / FADCp");
	peakCalibBl->GetYaxis()->SetTitle("Counts / au");
	peakCalibBl->GetYaxis()->SetTitleOffset(1.3);
	peakCalibBl->SetLineColor(kBlue);
	peakCalibBl->Draw();
	c1->Print("../plots/uubPeakCorrCalibZoom"+pmtId+stat+".pdf");
*/

	// ==============================================
	// *** Checking for first bin Peak Histograms ***
	
	/*
	int nx = histos->GetEntries();
	double xbl[nx];
	double ybl[nx];

	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) {
		histos->GetEntry(tmp);
		xbl[tmp] = tmp;
		ybl[tmp] = peakCalibBl->GetBinLowEdge(1);
	}

	TGraph* gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UUB First bin for Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since December 1st, 2020 until March 31st, 2021");
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../plots/uubFirstBinPeakCrrCalib"+pmtId+stat+".pdf");

*/
	//chargeCalibBl->Draw();
	//hcharge->Draw();

	//c1->Print("../plots/uubFitPk"+pmtId+stat+".pdf");

/*
	TCanvas *c2 = new TCanvas("c2","C2",500,10,1200,800);

	histos->GetEntry(0);
	c2->cd();
	//charge->SetLineColor(kBlack);
	//charge->Draw("SAME");
	fittedCh->GetYaxis()->SetTitle("Counts / au");
	fittedCh->GetYaxis()->SetTitleOffset(1.3);
	fittedCh->GetXaxis()->SetTitle("1 / FADCq");
	fittedCh->Draw("AL");
	c2->Print("../plots/uubFitCh"+pmtId+stat+".pdf");

	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xb[tmp] = tmp;
		yb[tmp] = chpk * (8.33/25);
	}

	TGraph* gr = new TGraph(nx,xb,yb);
	c1->cd();
	gr->SetTitle("UUB A/P "+pmtId+" "+stat+" Here A/P*(8.33/25)");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Events");
	gr->GetYaxis()->SetTitle("Counts / au");
	gr->GetYaxis()->SetRangeUser(0, 5);
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../plots/uubAp"+pmtId+stat+".pdf");
	*/
}
