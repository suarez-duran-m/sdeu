void runningThroughHistos()
{
  TFile *apr = TFile::Open("uubChPkPMT2St1740lrb35apr.root"); 
  TFile *may = TFile::Open("uubChPkPMT2St1740lrb35may.root");

  TTree *infoApr = (TTree*)apr->Get("ChargeData"); 
  TTree *infoMay = (TTree*)may->Get("ChargeData");

  TGraphErrors *grApr = new TGraphErrors(); 
  TGraphErrors *grMay = new TGraphErrors();

  infoApr->SetBranchAddress("graph", &grApr); 
  infoMay->SetBranchAddress("graph", &grMay);

  for (int getry=360; getry<365; getry++ )
  { 
    infoApr->GetEntry(getry); 
    grApr->Draw(); 
    for (int etry=450; etry<458; etry++ )
    { 
      infoMay->GetEntry(etry); 
      grMay->SetLineColor(kBlue); 
      grMay->Draw("same"); 
      gPad->WaitPrimitive();
    }
  }
}
