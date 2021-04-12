void plottingHunting(int pmtId){  
  TString filename;
  TString printname;
  if ( pmtId > 0 && pmtId < 4 ){
    filename.Form("zooTraces100binsPMT%d.root", pmtId);
    printname.Form("PMT%d", pmtId);
  }
  else if ( pmtId == 4 ){
    filename = "zooTraces100binsSPMT.root";
    printname = "SPMT";
  }
  else if ( pmtId == 5 ){
    filename = "zooTraces100binsPMTSSD.root";
    printname = "PMTSSD";
  }
  else{
    cout << "==================================================" << endl;
    cout << "Wrong Id for PMT, please introduce a valid PMT Id:" << endl;
    cout << "1 For PMT1; " << "2 For PMT2; " << "3 For PMT3; " 
      << "4 For SPMT; " << "5 For PMTSSD" << endl;
    cout << "==================================================" << endl;
    exit(0);
  }
  cout << "You have selected " << filename << endl;
    
  auto hFile = TFile::Open(filename);
  
  auto T = (TTree*)hFile->Get("T");

  auto pmtokh = (TH1F*)hFile->Get("pmthok");
  auto pmtokl = (TH1F*)hFile->Get("pmtlok");

  auto pmtunh = (TH2F*)hFile->Get("pmthun");
  auto pmtunl = (TH2F*)hFile->Get("pmtlun");

  auto pmtfrh = (TH2F*)hFile->Get("pmthfr");
  auto pmtfrl = (TH2F*)hFile->Get("pmtlfr");

  const char *stIds[19] = {"863", "1211", "1217", "1219", "1221", "1222", "1223", "1729", "1735", "1740", "1741", "1743", "1745", "1746", "1747", "1791", "1818", "1819", "1851"};

  auto hist = (TH2F*)hFile->Get("badTrac");
  //for(int e=1; e<hist->GetNbinsX(); e++)
    for(int i=1; i<hist->GetNbinsY()+1; i++)
      if ( hist->GetBinContent(5379, i) > 0 )
        cout << i << " " << hist->GetBinContent(5379, i) << endl;


  auto plain  = new TStyle("Plain","Plain Style (no colors/fill areas)");

  TCanvas *c0 = new TCanvas("c0", "2D", 0, 0,3600,2400);
  c0->cd();
  /*
  T->SetBranchAddress("badTrac", &trabad);
  T->GetEntry(13); 
  trabad->SetLineColor(kBlue);
  c0->Update();
  T->GetEntry(2);
  trabad->SetLineColor(kRed);
  trabad->Draw("SAME"); 
  c0->Print("trBad"+printname+"Hg.pdf");

  THStack *hs = new THStack("hs","Stacked 1D histograms");
  T->GetEntry(13);
  hs->Add(trabad);
  T->GetEntry(2);
  hs->Add(trabad);
  hs->Draw("nostack");
  */


  TCanvas *c1 = new TCanvas("c1", "2D", 0, 0,3600,2400);
  c1->cd();
  //c1->SetLeftMargin(0.11);
  //c1->SetRightMargin(0.135);
  //pmth->GetXaxis()->SetTitle("Days since December 1st 2020");
  //pmth->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmth->GetXaxis()->SetRangeUser(0,52);
  //pmth->GetXaxis()->SetTickLength(0.);
  //pmth->GetXaxis()->SetLabelSize(.02);
  //pmth->GetZaxis()->SetTitle("FADC");
  //pmth->GetZaxis()->SetRangeUser(180,300);
  //pmth->GetZaxis()->SetTickLength(0.);
  pmtokh->SetStats(0);
  //pmth->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmth->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtokh->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  //pmth->GetYaxis()->SetTitle("Station");
  //pmth->GetYaxis()->SetTitleFont(42);
  //pmth->GetYaxis()->SetTitleOffset(1.3);
  //pmth->GetYaxis()->SetLabelSize(.05);
  //plain->SetPalette(56);
  c1->SetGridy();
  c1->SetGridx();
  pmtokh->Draw("COLZ");
  //c1->Print("testMeanFirst100Hg.pdf");
  c1->Print("../../plots/trOk"+printname+"Hg.pdf");

  TCanvas c2("c2", "2D", 0,0,3600,2400);
  c2.cd();
  //c2.SetLeftMargin(0.11);
  //c2.SetRightMargin(0.135);
  //pmtl->GetXaxis()->SetTitle("Days since December 1st 2020");
  //pmtl->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmtl->GetXaxis()->SetRangeUser(0,52);
  //pmtl->GetXaxis()->SetTickLength(0.);
  //pmtl->GetXaxis()->SetLabelSize(.02);
  //pmtl->GetZaxis()->SetTitle("FADC");
  //pmtl->GetZaxis()->SetRangeUser(180,300);
  //pmtl->GetZaxis()->SetTickLength(0.);
  pmtokl->SetStats(0);
  //pmtl->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmtl->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtokl->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  //pmtl->GetYaxis()->SetTitle("Station");
  //pmtl->GetYaxis()->SetTitleFont(42);
  //pmtl->GetYaxis()->SetTitleOffset(1.3);
  //pmtl->GetYaxis()->SetLabelSize(.05);
  plain->SetPalette(56);
  c2.SetGridy();
  c2.SetGridx();
  pmtokl->Draw("COLZ");
  //c2.Print("testMeanLast100Hg.pdf");
  c2.Print("../../plots/trOk"+printname+"Lg.pdf");

  TCanvas c3("c3", "2D", 0,0,3600,2400);
  c3.cd();
  //c3.SetLeftMargin(0.11);
  //c3.SetRightMargin(0.135);
  //pmtunh->GetXaxis()->SetTitle("Days since December 1st 2020");
  //pmtunh->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmtunh->GetXaxis()->SetRangeUser(0,52);
  //pmtunh->GetXaxis()->SetTickLength(0.);
  //pmtunh->GetXaxis()->SetLabelSize(.02);
  //pmtunh->GetYaxis()->SetTitle("Station");
  //pmtunh->GetZaxis()->SetTitle("FADC");
  //pmtunh->GetZaxis()->SetRangeUser(0.,12.);
  pmtunh->SetStats(0);
  //pmtunh->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmtunh->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtunh->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  plain->SetPalette(56);
  c3.SetGridy();
  c3.SetGridx();
  pmtunh->Draw("COLZ");
  //c3.Print("testRmsFirst100Hg.pdf");
  c3.Print("../../plots/trUn"+printname+"Hg.pdf");

  TCanvas c4("c4", "2D", 0,0,3600,2400);
  c4.cd();
  //c4.SetLeftMargin(0.11);
  //c4.SetRightMargin(0.135);
  //pmtunl->GetXaxis()->SetTitle("Days since December 1st 2020");
  //pmtunl->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmtunl->GetXaxis()->SetRangeUser(0,52);
  //pmtunl->GetXaxis()->SetTickLength(0.);
  //pmtunl->GetXaxis()->SetLabelSize(.02);
  //pmtunl->GetYaxis()->SetTitle("Station");
  //pmtunl->GetZaxis()->SetTitle("FADC");
  //pmtunl->GetZaxis()->SetRangeUser(0.,12.);
  pmtunl->SetStats(0);
  //pmtunl->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmtunl->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtunl->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  plain->SetPalette(56);
  c4.SetGridy();
  c4.SetGridx();
  pmtunl->Draw("COLZ");
  //c4.Print("testRmsLast100Hg.pdf");
  c4.Print("../../plots/trUn"+printname+"Lg.pdf");

  TCanvas c5("c5", "2D", 0,0,3600,2400);
  c5.cd();
  //c5.SetLeftMargin(0.11);
  //c5.SetRightMargin(0.135);
  //pmtfrh->GetXaxis()->SetTitle("Days since December 1st 2020.");
  //pmtfrh->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmtfrh->GetXaxis()->SetRangeUser(0,52);
  //pmtfrh->GetXaxis()->SetTickLength(0.);
  //pmtfrh->GetXaxis()->SetLabelSize(.02);
  //pmtfrh->GetYaxis()->SetTitle("Station");
  //pmtfrh->GetZaxis()->SetTitle("FADC");
  //pmtfrh->GetZaxis()->SetRangeUser(0.,12.);
  pmtfrh->SetStats(0);
  //pmtfrh->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmtfrh->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtfrh->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  plain->SetPalette(56);
  c5.SetGridy();
  c5.SetGridx();
  pmtfrh->Draw("COLZ");
  //c5.Print("testRmsLast100Hg.pdf");
  c5.Print("../../plots/trFr"+printname+"Hg.pdf");


  TCanvas c6("c6", "2D", 0,0,3600,2400);
  c6.cd();
  //c6.SetLeftMargin(0.11);
  //c6.SetRightMargin(0.135);
  //pmtfrl->GetXaxis()->SetTitle("Days since December 1st 2020.");
  //pmtfrl->GetXaxis()->SetNdivisions(1020*3, "kTRUE");
  //pmtfrl->GetXaxis()->SetRangeUser(0,52);
  //pmtfrl->GetXaxis()->SetTickLength(0.);
  //pmtfrl->GetXaxis()->SetLabelSize(.02);
  //pmtfrl->GetYaxis()->SetTitle("Station");
  //pmtfrl->GetZaxis()->SetTitle("FADC");
  //pmtfrl->GetZaxis()->SetRangeUser(0.,12.);
  pmtfrl->SetStats(0);
  //pmtfrl->GetYaxis()->SetNdivisions(1020, "kTRUE");
  //pmtfrl->GetYaxis()->SetTickLength(0.);
  for ( int i=0; i<19; i++)
    pmtfrl->GetXaxis()->SetBinLabel(i+1, stIds[i]);
  plain->SetPalette(56);
  c6.SetGridy();
  c6.SetGridx();
  pmtfrl->Draw("COLZ");
  //c6.Print("testRmsLast100Hg.pdf");
  c6.Print("../../plots/trFr"+printname+"Lg.pdf");
}
