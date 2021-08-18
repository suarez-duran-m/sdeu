void plottingSingle(){  
 
  auto hFile = TFile::Open("zooTraces100binsPMT1.root");
  
  auto T = (TTree*)hFile->Get("T");

	auto sgl = (TH1D*)hFile->Get("sglEvt");

  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

	auto copy1 = (TH1D*)sgl->Clone();
	//copy1->Reset();
	
	for(Int_t i = 101; i <= 1947; ++i)
		copy1->SetBinContent(i, 0);

  TCanvas *c1 = new TCanvas("c1", "2D", 1, 35,3600,2400);
  c1->cd();

	gStyle->SetOptTitle(0);
	sgl->SetStats(0);
	sgl->SetLineWidth(2.0);
	sgl->GetYaxis()->SetTitle("FADC");

	sgl->GetYaxis()->SetTickLength(0.);
	sgl->GetYaxis()->SetLabelSize(.03);
	sgl->GetYaxis()->SetTitleFont(42);
	sgl->GetXaxis()->SetTitle("Time / (8.33 ns)");
	sgl->GetXaxis()->SetRangeUser(0, 2047);
	sgl->GetYaxis()->SetRangeUser(100, 750);
	sgl->GetXaxis()->SetTickLength(0.);
 	sgl->GetXaxis()->SetLabelSize(.03);
	sgl->Draw("SAME");  

  copy1->SetStats(0);
  copy1->GetYaxis()->SetTickLength(0.);
	copy1->SetLineColor(kRed);
	copy1->SetLineWidth(2.0);
  copy1->Draw("SAME");
  c1->Print("../../plots/sglEvt.pdf");

}
