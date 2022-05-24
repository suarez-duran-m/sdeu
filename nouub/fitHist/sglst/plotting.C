void plotting(int pmt, int st) {

	TString pmtId;
	TString stat;

	pmtId.Form("PMT%d", pmt);
	stat.Form("St%d", st);

	TFile *f = TFile::Open("ubCalibHist"+pmtId+stat+".root");
	cout << "Reading file: " << "ubCalibHist"+pmtId+stat+".root" <<endl;

 	TTree *histos = (TTree*)f->Get("Histograms"); 
	TH1F *peak = new TH1F(); 
	TH1F *charge = new TH1F();
	TH1F *rawpeak = new TH1F(); 
	TH1F *rawcharge = new TH1F();
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
 	
	histos->SetBranchAddress("corrPk",&peak);	
	histos->SetBranchAddress("corrCh",&charge);
	histos->SetBranchAddress("recePk",&rawpeak);
	histos->SetBranchAddress("receCh",&rawcharge);
	histos->SetBranchAddress("fllHPk", &fittedPk); 
	histos->SetBranchAddress("fllHCh", &fittedCh);
	histos->SetBranchAddress("corrPkCalibBl",&peakCalibBl);
	histos->SetBranchAddress("corrChCalibBl",&chargeCalibBl);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("offSetCh", &offsetch);
	histos->SetBranchAddress("ap", &chpk);
	histos->SetBranchAddress("baselineHbase", &blHbase);
	histos->SetBranchAddress("baselineCalib", &blCalib);
	histos->SetBranchAddress("entryEvt", &evtTime);

	//TCanvas *c0 = new TCanvas("c0","C0",500,10,1200,800);
	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);

	TGraph *gr;
	int nx = histos->GetEntries();
	double xbl[nx];
	double ybl[nx];

	TString p0nm;
	TString p1nm;
	TPaveStats *ptstats;
	TText *ptstats_LaTex;

	histos->GetEntry(0);
  cout << chpk << endl;
  c1->cd();
  fittedPk->Draw("AP");
  /*
	//c0->cd();
	//rawpeak->Draw();
	c1->cd();
  gStyle->SetOptTitle(0); 
  gStyle->SetOptStat(0);
  charge->GetYaxis()->SetTitle("Counts [au]");
  charge->GetYaxis()->SetLabelFont(42);
  //charge->GetYaxis()->SetTitleOffset(1);
  charge->GetYaxis()->SetTitleSize(1.2);
  charge->GetYaxis()->SetTitleFont(42);
  charge->GetXaxis()->SetTitle("[FADC]");
  charge->Draw("");
  TPaveText *t = new TPaveText(0.3497496,0.922179,0.639399,0.997406, "brNDC"); // left-up
  t->AddText("UB Charge histogram Offset - Baseline "+pmtId);
  t->SetFillColor(0);
  t->SetFillStyle(0);
  t->SetLineColor(0);
  t->SetLineWidth(15);
  t->SetTextFont(42);
  t->SetTextSize(0.05);
  t->Draw();
  c1->Modified();
  c1->cd();
  c1->SetSelected(c1);
*/
  //title->Draw(); 
  //charge->SetTitle("UB Charge histogram Offset - Baseline "+pmtId);

	//rawcharge->Draw();


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
	gr->SetTitle("UB Offset for Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since September 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	if ( pmt==1 )
		gr->GetYaxis()->SetRangeUser(55, 60);
	else if ( pmt==2 )
		gr->GetYaxis()->SetRangeUser(40, 43);
	else
		gr->GetYaxis()->SetRangeUser(50, 54);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetPk"+pmtId+stat+".pdf");
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
	gr->SetTitle("UB Offset for Charge Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since September 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	
	if ( pmt==1 )
		gr->GetYaxis()->SetRangeUser(1150, 1170);
	else if ( pmt==2 )
		gr->GetYaxis()->SetRangeUser(830, 843);
	else
		gr->GetYaxis()->SetRangeUser(1050, 1062);
		
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.8);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubOffsetCh"+pmtId+stat+".pdf");
*/


	// =============================================
	// *** Counts in first bin of raw histograms ***
/*
	for (Int_t i=0; i<nx; i++) {
		histos->GetEntry(i);
		xbl[i] = evtTime;
		ybl[i] = rawpeak->GetBinContent(1);
	}

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("Counts in first bin for UB Peak Histograms "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1,1);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.2);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubCntFirstBinPk"+pmtId+stat+".pdf");
*/


/*
	cout << histos->GetEntries() << endl;
	histos->GetEntry(0);
	c0->cd();
	fittedPk->Draw("AP");
	//rawpeak->Draw();
	c1->cd();
	fittedCh->Draw("AP");
	//rawcharge->Draw();
	
	cout << event << endl; 

	double sum = 0.;
	int cnt = 0;
	for ( int etr=0; etr<histos->GetEntries(); etr++ )
	{
		histos->GetEntry(etr);
		if ( chpk != 0 )
			cnt++;
		sum += chpk;
	}
	cout << sum / cnt << endl;

	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t etr=0; etr<nx; etr++ ) 
	{
		histos->GetEntry( etr );
		xb[etr] = event;
		yb[etr] = chpk;
	}

	gr = new TGraph(nx,xb,yb);
	c1->cd();
	gr->SetTitle("UB First bin Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetXaxis()->SetTitle("Since September 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitleOffset(1.);
	gr->GetYaxis()->SetRangeUser(2.5, 3.5);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
*/

	// ================
	// *** Baseline ***
/*
	double minlimHb = 1e4;
	double maxlimHb = 0.;
	double minlimCl = 1e4;
	double maxlimCl = 0.;

	TString p0nm;
	TString p1nm;
	TPaveStats *ptstats;
	TText *ptstats_LaTex;

	TH1F *histBlHbase = new TH1F("histBlHbase", "UB Mean Baseline HBase "+pmtId+" "+stat, 100, 0, 100);
	TH1F *histBlCalib = new TH1F("histBlCalib", "UB Baseline Calib "+pmtId+" "+stat, 100, 0, 100);

	TF1 *linef = new TF1("linef", "[1]*x + [0]");

	nx = 0;
	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ )
	{
		histos->GetEntry( tmp );
		if ( blHbase > 0 )
			nx++;
	}
	xbl[nx];
	ybl[nx];
	int tmp0 = 0;

	cerr << nx << " " << histos->GetEntries() << endl;

	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) 
	{
		histos->GetEntry( tmp );
		if ( blHbase > 0 )
		{
			xbl[tmp0] = evtTime;
			ybl[tmp0] = blHbase;
			tmp0++;
			histBlHbase->Fill( blHbase );
			if ( blHbase < minlimHb )
				minlimHb = blHbase;
			`if ( blHbase > maxlimHb )
				maxlimHb = blHbase;
		}
	}
	minlimHb -= 1.;
	maxlimHb += 2.;

	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB Mean Baseline HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since September 1st, 2020 (month/day)");
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
	c1->Print("../../plots/ubBlHbase"+pmtId+stat+".pdf");
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
	gr->SetTitle("UB Peak histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1.5, 1.5);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubPeakFirstBinHBase"+pmtId+stat+".pdf");
	*/

/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB Peak histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-1., 2.);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubPeakFirstBinCalibBase"+pmtId+stat+".pdf");
*/

	// =============================
	// *** For Charge Histograms ***

	histos->GetEntry(0);
	//cerr << blHbase << " " << endl;
	/*
	c1->cd();
	charge->SetTitle("UB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 58");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../../plots/ubChargeCoorHBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. HBase: 58");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../../plots/ubChargeCoorHBaseZoom"+pmtId+stat+".pdf");
*/

	/*
	cerr << blCalib << " " << endl;
	c1->cd();
	charge->SetTitle("UB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 57");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../../plots/ubChargeCoorCalibBase"+pmtId+stat+".pdf");

	c1->cd();	
	charge->SetTitle("UB Charge-Histogram "+pmtId+" "+stat+" Ev.: 61219267 Corr. Calib.Base: 57");
	charge->SetStats(0);
	charge->SetLineColor(kBlue);
	charge->GetXaxis()->SetTitle("1. / FADC");
	charge->GetXaxis()->SetRangeUser(charge->GetXaxis()->GetBinLowEdge(1), charge->GetXaxis()->GetBinLowEdge(30));
	charge->GetYaxis()->SetTitle("Counts / au");
	charge->GetYaxis()->SetTitleOffset(1.3);
	charge->Draw();
	c1->Print("../../plots/ubChargeCoorCalibBaseZoom"+pmtId+stat+".pdf");
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
	gr->SetTitle("UB Charge histogram First Bin Center HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-10, 10);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubChargeFirstBinHBase"+pmtId+stat+".pdf");
	*/
	

	/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB Charge histogram First Bin Center Calib.Base "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->GetYaxis()->SetRangeUser(-10, 30);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.6);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubChargeFirstBinCalibBase"+pmtId+stat+".pdf");
*/


	// =============================
	// *** Area/Peak Calculation ***
/*
	nx = histos->GetEntries();
  double tmpap = 0.;
	for ( Int_t tmp=0; tmp<nx; tmp++ ) 
	{
		histos->GetEntry( tmp );
		//xbl[tmp] = evtTime;
		//ybl[tmp] = chpk;
    tmpap += chpk;
	}
  cout << tmpap/nx << endl;
  */
/*	
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB A/P With Offset and Baseline removed "+pmtId+" "+stat+" HBase");
	gr->GetXaxis()->SetTitle("Time since August 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("AoP [25*ns]");
	//gr->GetYaxis()->SetRangeUser(0, 10);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c1->Print("../../plots/ubApOffset"+pmtId+stat+".pdf");
*/

	/*
	gr = new TGraph(nx,xbl,ybl);
	c1->cd();
	gr->SetTitle("UB A/P "+pmtId+" "+stat+" Calib.Base");
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Since December 1st, 2020 (month/day)");
	gr->GetXaxis()->SetTimeDisplay(1);
	gr->GetXaxis()->SetTimeFormat("%m/%d");
	gr->GetYaxis()->SetTitle("1. / 25*ns");
	gr->GetYaxis()->SetRangeUser(0, 10);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c1->Print("../../plots/ubApOffsetCalib"+pmtId+stat+".pdf");
*/




	/*
	// =====================================
	// *** Looking for threshols in 		 ***
	// *** firsts bin of Peak histograms ***

	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t etr=0; etr<nx; etr++ ) 
	{
		histos->GetEntry( etr );
		xb[etr] = etr;
		for ( Int_t b=0; b<20; b++ )
			if ( peakCalibBl->GetBinContent(b+1) != 0 ) 
			{
				yb[etr] = peakCalibBl->GetBinLowEdge(b+1);
				break;
			}
	}

	TGraph* gr = new TGraph(nx,xb,yb);
	c0->cd();
	gr->SetTitle("UB First bin Peak Histogram BL from Calib.Base "+pmtId+" "+stat);
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetYaxis()->SetTitle("1 / FADCp");
	gr->GetXaxis()->SetTitle("Events since September 1st, 2020 until November 30th, 2020");
	gr->GetYaxis()->SetTitleOffset(1.);
	gr->GetYaxis()->SetRangeUser(0, 6);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	c0->Print("../../plots/ubFirstBinPeakCrrCalib"+pmtId+stat+".pdf");
*/
	/*
	c0->cd();
	peakCalibBl->GetXaxis()->SetRangeUser(0,1200);
	peakCalibBl->Draw();
	*/

	/*		
	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);
	c1->cd();
	fittedPk->GetYaxis()->SetTitle("Counts / au");
	fittedPk->GetYaxis()->SetTitleOffset(1.3);
	fittedPk->GetXaxis()->SetTitle("1 / FADCp");
	fittedPk->Draw("AL");
	//c1->Print("../../plots/ubFitPk"+pmtId+stat+".pdf");
	*/

/*
	histos->GetEntry(17);
	cout << "Offset: " << offsetCh << " " << endl  
		<< chargeCalibBl->GetBinContent(0) << " " 
		<< chargeCalibBl->GetBinCenter(0) << " "
		<< chargeCalibBl->GetBinContent(1) << " "
	 	<< chargeCalibBl->GetBinCenter(1) << " "	
		<< endl;
	TCanvas *c2 = new TCanvas("c2","C2",500,10,1200,800);

	c2->cd();
	//chargeCalibBl->Draw();
	fittedCh->Draw("AL");
	//fittedPk->Draw("AL");
	//chargeCalibBl->Draw();
		
	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xb[tmp] = tmp;
		yb[tmp] = chpk;
		if ( chpk > 3.5 && chpk < 3.6 )
			cout << tmp << endl;
	}

	TGraph* gr = new TGraph(nx,xb,yb);
	c0->cd();
	gr->SetTitle("UB A/P "+pmtId+" "+stat);
	//gr->GetXaxis()->SetRangeUser(1e4, 5e4);
	gr->GetXaxis()->SetTitle("Events");
	gr->GetYaxis()->SetTitle("Counts / au");
	gr->GetYaxis()->SetTimeOffset(1.3);
	gr->GetYaxis()->SetRangeUser(0, 5);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(0.4);
	gr->SetMarkerColor(4);
	gr->Draw("AP");
	//c0->Print("../../plots/ubAp"+pmtId+stat+".pdf");
*/
}
