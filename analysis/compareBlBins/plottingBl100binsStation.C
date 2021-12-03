void plottingBl100binsStation(int hg, int StId){ 
  TString filename;
  TString printname;
  TString hgname;

  if ( hg !=0 && hg !=1 ){
    cout << hg << endl;
    cerr << "Please select HG (1) or LG (0) option" 
      << endl;
    exit(0);
  }

	filename.Form("bl100bins%d.root", StId);
	printname.Form("blDiffMeanSt%d", StId);
  cout << "You have selected " << filename << " "
    << "with HG option: " << hg << endl;
 
  auto hFile = TFile::Open(filename);

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};

  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

 /* 
  TH1F *st1740 = new TH1F("st1740","Distribution of The Difference of The Mean for The Station 1740's PMT1 HG", 320,-1, 30);
  for(int i=0; i<pmtdm->GetXaxis()->GetNbins(); i++){
    if ( pmtdm->GetBinContent(i+1, 10)!=0 )
      st1740->Fill(pmtdm->GetBinContent(i+1, 10));
  }
*/

  TH1F *pmt1 = (TH1F*)hFile->Get("stdiffDist");
  TH1F *pmt2 = (TH1F*)hFile->Get("stdiffDist2");
  TH1F *pmt3 = (TH1F*)hFile->Get("stdiffDist3");

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(kFALSE);
  TCanvas c7("c7", "2D", 0,0,3600,2400);

  c7.SetBottomMargin(0.18);
  c7.cd();
  pmt1->GetYaxis()->SetTitle("Counts/au");
  pmt1->GetXaxis()->SetTitle("Difference of Mean / FADC");

  pmt3->SetLineColor(880);
	pmt3->SetLineWidth(2);
  pmt3->Draw();
  pmt2->SetLineColor(kBlue);
	pmt2->SetLineWidth(2);
  pmt2->Draw("SAME");
  pmt1->SetLineColor(kRed);
	pmt1->SetLineWidth(2);
  pmt1->Draw("SAME");

  auto legend = new TLegend(0.1,0.7,0.4,0.9);
	
	TString pmt1m = to_string( pmt1->GetMean() );
	TString pmt2m = to_string( pmt2->GetMean() );
	TString pmt3m = to_string( pmt3->GetMean() );
	TString pmt1r = to_string( pmt1->GetRMS() );
	TString pmt2r = to_string( pmt2->GetRMS() );
	TString pmt3r = to_string( pmt3->GetRMS() );

  legend->AddEntry(pmt1,"PMT1 HG");
  legend->AddEntry((TObject*)0, "Mean: "+pmt1m+"; RMS: 3.32", "");
  legend->AddEntry(pmt2,"PMT2 HG");
  legend->AddEntry((TObject*)0, "Mean: 0.71; RMS: 1.32", "");
  legend->AddEntry(pmt3,"PMT3 HG");
  legend->AddEntry((TObject*)0, "Mean: 0.87; RMS: 1.21", "");
  //legend->AddEntry("gr","Graph with error bars","lep");
  legend->Draw();
  c7.Print("../../plots/"+printname+".pdf");
}
