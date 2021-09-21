TCanvas *canvasStyle(TString name) {
  TCanvas *canvas = new TCanvas(name, name, 1600, 900);
  canvas->SetBorderMode(0);
  canvas->SetBorderSize(2);
  canvas->SetLeftMargin(0.11); 
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.02); 
  canvas->SetBottomMargin(0.15);
  canvas->SetFrameBorderMode(0);
  return canvas;
}

void readAscii()
{
  //diffHighUb.dat  diffHighUub.dat  diffLowUb.dat  diffLowUub.dat
  TH1F *diffLowUb = new TH1F("diffLowUb", "chargeRaw - chargeLow", 3000, 0., 3000.);
  TH1F *diffLowUub = new TH1F("diffLowUub", "chargeRaw - chargeLow", 3000, 0., 3000.);
  TH1F *diffHighUb = new TH1F("diffHighUb", "chargeRaw - chargeHigh", 3000, 0., 3000.);
  TH1F *diffHighUub = new TH1F("diffHighUub", "chargeRaw - chargeHigh", 3000, 0., 3000.);

  ifstream inp;
  double x;
  double y;
  int nlines = 2000;
  inp.open("diffLowUb.dat");
  for (int i=1; i<=nlines; i++) {
    inp >> x >> y;
    diffLowUb->SetBinContent(x,y);
  }
  inp.close();
  inp.open("diffLowUub.dat");
  for (int i=1; i<=nlines; i++) {
    inp >> x >> y;
    diffLowUub->SetBinContent(x,y);
  }
  inp.close();

  inp.open("diffHighUb.dat");
  for (int i=1; i<=nlines; i++) {
    inp >> x >> y;
    diffHighUb->SetBinContent(x,y);
  }
  inp.close();
  inp.open("diffHighUub.dat");
  for (int i=1; i<=nlines; i++) {
    inp >> x >> y;
    diffHighUub->SetBinContent(x, y);
  }
  inp.close();

  TLegend *leg;
  TString strMu;
  TString strRms;

  TCanvas *c1 = new TCanvas("c1", "Low");
  c1->cd();
  diffLowUb->SetStats(kFALSE);
  diffLowUb->Draw();

  diffLowUub->SetLineColor(kGreen+3);
  diffLowUub->Draw("same");


  leg = new TLegend(0.4,0.6,0.62,0.85);
  strMu.Form("%.2f", diffLowUb->GetMean());
  strRms.Form("%.2f", diffLowUb->GetRMS());
  leg->AddEntry(diffLowUb, "For UB, #mu: "+strMu+" RMS: "+strRms);
  strMu.Form("%.2f", diffLowUub->GetMean());
  strRms.Form("%.2f", diffLowUub->GetRMS());
  leg->AddEntry(diffLowUub, "For UUB, #mu: "+strMu+" RMS: "+strRms);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2", "High");
  c2->cd();
  diffHighUb->SetStats(kFALSE);
  diffHighUb->GetYaxis()->SetRangeUser(0, 1.7e3);
  diffHighUb->Draw();

  diffHighUub->SetLineColor(kGreen+3);
  diffHighUub->Draw("same");
  
  leg = new TLegend(0.4,0.6,0.62,0.85);
  strMu.Form("%.2f", diffHighUb->GetMean());
  strRms.Form("%.2f", diffHighUb->GetRMS());
  leg->AddEntry(diffLowUb, "For UB, #mu: "+strMu+" RMS: "+strRms);
  strMu.Form("%.2f", diffHighUub->GetMean());
  strRms.Form("%.2f", diffHighUub->GetRMS());
  leg->AddEntry(diffLowUub, "For UUB, #mu: "+strMu+" RMS: "+strRms);
  leg->SetTextSize(0.04);
  leg->SetBorderSize(0);
  leg->Draw();
}
