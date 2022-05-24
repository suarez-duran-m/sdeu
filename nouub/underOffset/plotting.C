void plotting(int pmt, int st) {

	TString pmtId;
	TString stat;

	pmtId.Form("PMT%d", pmt);
	stat.Form("St%d", st);

	TFile *f = TFile::Open("ubOffset"+pmtId+stat+".root");
	cout << "Reading file: " << "ubCalibHist"+pmtId+stat+".root" <<endl;

 	TTree *histos = (TTree*)f->Get("Histograms");

	double blHbase; 
	double blCalib; 
	int offsetpk;
	int offsetCh;

	histos->SetBranchAddress("baselineHbase", &blHbase);
	histos->SetBranchAddress("baselineCalib", &blCalib);
	histos->SetBranchAddress("offSetPk", &offsetpk);
	histos->SetBranchAddress("offSetCh", &offsetCh);


	// ====================
	// *** For Baseline ***
	
	TCanvas *c0 = new TCanvas("c0","C0",500,10,1200,800);
	int nx = 0;
	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) // Ignoring outliers values
	{
		histos->GetEntry( tmp );
		if ( blHbase<80 && blHbase>1 )
			nx++;
	}

	double xbl[nx];
	double ybl[nx];
	double minlimHb = 1e4;
	double maxlimHb = 0.;
	double minlimCl = 1e4;
	double maxlimCl = 0.;
	TGraph* gr;

	for (Int_t i=0; i<nx; i++) 
	{






























/*
	TH1F *histBlHbase = new TH1F("histBlHbase", "UB Mean Baseline HBase "+pmtId+" "+stat, 500, 0, 500);
	TH1F *histBlCalib = new TH1F("histBlCalib", "UB Baseline Calib "+pmtId+" "+stat, 500, 0, 500);
	TF1 *linef = new TF1("linef", "[1]*x + [0]"); // Linear Fit function


	// *** HBase ***
	int tmpb = 0;
	for ( Int_t etr=0; etr<histos->GetEntries(); etr++ )
	{
		histos->GetEntry( etr );
		histBlHbase->Fill( blHbase );
		if ( blHbase>80 || blHbase==0 )
			continue;
		xbl[tmpb] = tmpb;
		ybl[tmpb] = blHbase;
		if ( blHbase < minlimHb )
			minlimHb = blHbase;
		if ( blHbase > maxlimHb )
			maxlimHb = blHbase;
		tmpb++;
	}
	minlimHb -= 1.;
	maxlimHb += 2.;

 	gr = new TGraph(nx,xbl,ybl);

	c0->cd();
	gr->SetTitle("UB Mean Baseline HBase "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st, 2020 until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Mean Baseline HBase / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(8);
	gr->SetMarkerSize(0.8);
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
	c0->Print("../../plots/ubBlHbase"+pmtId+stat+".pdf");
	
	
	// *** Calib.Base ***
	nx = 0;
	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) // Ignoring outliers values
	{
		histos->GetEntry( tmp );
		if ( blHbase<80 && blHbase>1 )
			nx++;
	}
	double xblB[nx];
	double yblB[nx];
	tmpb = 0;
	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) {
		histos->GetEntry( tmp );
		histBlCalib->Fill( blCalib );
		if ( blCalib>80 || blCalib==0 )
			continue;
		xblB[tmpb] = tmpb;
		yblB[tmpb] = blCalib;
		if ( blCalib < minlimCl )
			minlimCl = blCalib;
		if ( blCalib > maxlimCl )
			maxlimCl = blCalib;
		tmpb++;
	}

	minlimCl -= 1.;
	maxlimCl += 2.;

	gr = new TGraph(nx, xblB, yblB);
	c0->cd();
	gr->SetTitle("UB Baseline Calib "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events since September 1st, 2020 until November 30th, 2020");
	gr->GetYaxis()->SetTitle("Calib.Base / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(8);
	gr->SetMarkerSize(0.8);
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
	c0->Print("../../plots/ubBlCalib"+pmtId+stat+".pdf");


	// ===========================
	// *** For Baseline histos ***

	c0->cd();
	histBlHbase->SetLineWidth(2.0);
	histBlHbase->GetXaxis()->SetRangeUser((int)minlimHb, (int)maxlimHb);	
	histBlHbase->SetStats(kTRUE);
	histBlHbase->GetYaxis()->SetTitle("Counts / au");
	histBlHbase->GetXaxis()->SetTitle("1 / FADC");
	histBlHbase->Draw();
	c0->Print("../../plots/ubBlHistHbase"+pmtId+stat+".pdf");

	c0->cd();
	histBlCalib->SetLineWidth(2.0);
	histBlCalib->GetXaxis()->SetRangeUser((int)minlimCl, (int)maxlimCl);
	histBlCalib->SetStats(kTRUE);
	histBlCalib->GetYaxis()->SetTitle("Counts / au");
	histBlCalib->GetXaxis()->SetTitle("1 / FADC");
	histBlCalib->Draw();
	c0->Print("../../plots/ubBlHistCalib"+pmtId+stat+".pdf");
	*/

/*
	// ==================
	// *** For Offset ***

	nx = 0;
	for ( Int_t tmp=0; tmp<histos->GetEntries(); tmp++ ) // Ignoring outliers values
	{
		histos->GetEntry( tmp );
		if ( blHbase<80 && blHbase>1 )
			nx++;
	}

	int n = 0;
	int bincnt = 0;
	double cmOffset = 0.;
	double aveEvent = 20.;
	int nxave = int(nx/aveEvent);
	double xosav[nxave];
	double yosav[nxave];
	nx = histos->GetEntries();

	for ( Int_t tmp=0; tmp<nx; tmp++ )
	{
		histos->GetEntry( tmp );
		if ( blCalib>80 || blCalib==0 )
			continue;
		cmOffset += offsetpk;
		n++;
		if ( n==int(aveEvent) )
		{
			xosav[bincnt] = bincnt;
			yosav[bincnt] = cmOffset/aveEvent;
			cmOffset = 0.;
			n = 0;
			bincnt++;
		}
	}

	gr = new TGraph(nxave, xosav, yosav);
	c0->cd();
	gr->SetTitle("UB Offset for Peak histograms Average per 20 events "+pmtId+" "+stat);
	gr->GetXaxis()->SetTitle("Events/20; Since September 1st, 2020 until November 30ht, 2020");
	gr->GetYaxis()->SetTitle("Offset / FADC");
	gr->GetYaxis()->SetTitleOffset(1.3);
	//gr->GetYaxis()->SetRangeUser(272.5,275);
	gr->SetFillStyle(1000);
	gr->SetMarkerStyle(20);
	gr->SetMarkerSize(1.0);
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
	c0->Print("../../plots/ubDiffBlHbaseCalib"+pmtId+stat+".pdf");
*/
}
