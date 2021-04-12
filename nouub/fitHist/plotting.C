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
	TH1F *peakCalibBl = new TH1F();
 	TGraphErrors *fittedPk = new TGraphErrors(); 
 	TGraphErrors *fittedCh = new TGraphErrors(); 

	double chpk;
 	
	histos->SetBranchAddress("corrPk",&peak);	
	histos->SetBranchAddress("corrCh",&charge);
	histos->SetBranchAddress("corrPkCalibBl",&peakCalibBl);
	histos->SetBranchAddress("fllHPk", &fittedPk); 
	histos->SetBranchAddress("fllHCh", &fittedCh);
	histos->SetBranchAddress("ap", &chpk);

	histos->GetEntry(0);

	TCanvas *c0 = new TCanvas("c0","C0",500,10,1200,800);

	c0->cd();
	peakCalibBl->GetXaxis()->SetRangeUser(0,1200);
	peakCalibBl->Draw();


	/*
	TCanvas *c1 = new TCanvas("c1","C1",500,10,1200,800);

	c1->cd();
 	//peak->SetLineColor(kBlack);
  //peak->Draw("SAME"); 
	fittedPk->GetYaxis()->SetTitle("Counts / au");
	fittedPk->GetYaxis()->SetTitleOffset(1.3);
	fittedPk->GetXaxis()->SetTitle("1 / FADCp");
	fittedPk->Draw("AL");
	c1->Print("../../plots/ubFitPk"+pmtId+stat+".pdf");

	TCanvas *c2 = new TCanvas("c2","C2",500,10,1200,800);

	c2->cd();
	//charge->SetLineColor(kBlack);
	//charge->Draw("SAME");
	fittedCh->GetYaxis()->SetTitle("Counts / au");
	fittedCh->GetYaxis()->SetTitleOffset(1.3);
	fittedCh->GetXaxis()->SetTitle("1 / FADCq");
	fittedCh->Draw("AL");
	c2->Print("../../plots/ubFitCh"+pmtId+stat+".pdf");

	int nx = histos->GetEntries();
	double xb[nx];
	double yb[nx];

	for ( Int_t tmp=0; tmp<nx; tmp++ ) {
		histos->GetEntry( tmp );
		xb[tmp] = tmp;
		yb[tmp] = chpk;
	}

	TGraph* gr = new TGraph(nx,xb,yb);
	c1->cd();
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
	c1->Print("../../plots/ubAp"+pmtId+stat+".pdf");
*/
}
