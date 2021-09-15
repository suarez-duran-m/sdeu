void runningThroughHistos()
{
  TFile *apr = TFile::Open("uubChPkPMT3St1223lrb35aug.root");
  //TFile *apr = TFile::Open("uubChPkPMT2St1740lrb35apr.root"); 
  //TFile *may = TFile::Open("uubChPkPMT2St1740lrb35may.root");

  TTree *infoApr = (TTree*)apr->Get("ChargeData"); 
  //TTree *infoMay = (TTree*)may->Get("ChargeData");

  TGraphErrors *grApr = new TGraphErrors(); 
  //TGraphErrors *grMay = new TGraphErrors();
  int evt = 0;

  infoApr->SetBranchAddress("graph", &grApr); 
  infoApr->SetBranchAddress("eventId", &evt); 
  //infoMay->SetBranchAddress("graph", &grMay);
  
  TLegend *leg;
  TString strEvt;
  TCanvas *c0 = new TCanvas("c0", "c0", 1600, 900);

  for (int getry=0; getry<infoApr->GetEntries(); getry++ )
  {
    c0->cd();
    c0->Clear();
    infoApr->GetEntry(getry);
    strEvt.Form("%d", evt);
    grApr->GetYaxis()->SetRangeUser(0, 160);
    grApr->Draw();

    TF1 *poly;
    poly = grApr->GetFunction("poly2");        
    poly->SetLineColor(kRed);
    poly->SetLineWidth(2);
    poly->Draw("same");    
  
    leg = new TLegend(0.5,0.31,0.92,0.5);
    leg->AddEntry(infoApr, "Event number: "+strEvt,"");
    leg->SetTextSize(0.06);
    leg->SetBorderSize(0);
    leg->Draw();
    gPad->WaitPrimitive();
    /*
    for (int etry=450; etry<458; etry++ )
    { 
      infoMay->GetEntry(etry); 
      grMay->SetLineColor(kBlue); 
      grMay->Draw("same"); 
      gPad->WaitPrimitive();
    }
    */
  }
}
