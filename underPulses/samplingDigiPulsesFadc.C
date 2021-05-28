#include "Riostream.h"
void samplingDigiPulsesFadc() 
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("samplingDigiPulsesFadc.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%spulse_muon.dat",dir.Data()));

  Float_t x,y;
  TFile *f = new TFile("basic.root","RECREATE");
  TH1D *h = new TH1D("","",400, 0, 400);

  vector < double > xbins;
  vector < double > ycnts;
  vector < double > yerrs;

  for ( int line = 0; line<40; line++ )
  {
    in >> x >> y;
    if (!in.good())
      break;
    xbins.push_back( x );
    ycnts.push_back( y );
    yerrs.push_back( sqrt(ycnts[line]) );
    h->Fill(x, y);
  }
  in.close();

  TGraphErrors* chFit = new TGraphErrors( 
      xbins.size(), &xbins.front(), &ycnts.front(), 0, &yerrs.front() );

  chFit->Fit("landau", "","", 80, 350);
  TF1 *myfunc;
  myfunc = chFit->GetFunction("landau");
  double meanp0 = myfunc->GetParameter(0);
  double sigmp0 = myfunc->GetParError(0);
  double meanp1 = myfunc->GetParameter(1);
  double sigmp1 = myfunc->GetParError(1);
  double meanp2 = myfunc->GetParameter(2);
  double sigmp2 = myfunc->GetParError(2);

  double tmpp0 = 0.;
  double tmpp1 = 0.;
  double tmpp2 = 0.;
  
  meanp0 = myfunc->GetParameter(0);
  sigmp0 = myfunc->GetParError(0);
  meanp1 = 9.55307e+01;
  sigmp1 = myfunc->GetParError(1);
  meanp2 = 1.5537e+01;
  sigmp2 = myfunc->GetParError(2);
    
  int Npulses = 1e4;
  double fadc2mvUb = 2000./1024.;//1.95;//2000./1024.; // mV
  double fadc2mvUub = 0.49;//2000./4096.;//0.49;//2000./4096.; // mV

  double factorChUb = 25./50.; // ns / Ohms
  double factorChUub = 8.33/50.; // ns / Ohms

  TH1F **pulse = new TH1F *[Npulses];
  TString pulId;

  TH1F *histPkUb = new TH1F("histPkUb","UB", 200, 0, 200); // In mV
  TH1F *histChUb = new TH1F("histChUb","UB", 2000, 0, 2000); // In pC

  TH1F *histPkUub = new TH1F("histPkUub","UUB", 200, 0, 200); // In mV
  TH1F *histChUub = new TH1F("histChUub","UUB", 2000, 0, 2000); // In pC

  TH1F *histPk = new TH1F("histPk","1. ns", 200, 0, 200); // In mV
  TH1F *histCh = new TH1F("histCh","1. ns", 7000, 0, 7000); // In pC

  double tmpsgnl = 0.;
  double tmpPk = 0;
  int tmpCh = 0;

  double tmppk1ns = 0.;
  double tmpch1ns = 0.;

  for (int pl = 0; pl<Npulses; pl++) // N pulses
  {
    tmpp0 = gRandom->Gaus(meanp0, sigmp0/2.);
    tmpp1 = gRandom->Gaus(meanp1, sigmp1/2.);
    tmpp2 = gRandom->Gaus(meanp2, sigmp2/2.);
    pulId.Form("%d", pl);
    pulse[pl] = new TH1F(pulId,"", 400, 0 , 400);
  
    tmpPk = 0;
    tmpCh = 0;
    for (int dt = 0; dt<425; dt+=25)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      if( tmpPk < tmpsgnl )
        tmpPk = tmpsgnl; // Storing in mV
      tmpCh += int(tmpsgnl/fadc2mvUb); // Storing in FADC
    }
    tmpPk /= fadc2mvUb; // Converting to FADC
    histPkUb->Fill( int(tmpPk) ); // Storing in FADC
    histChUb->Fill( tmpCh ); // Sotoring in FADC

    tmpPk = 0;
    tmpCh = 0;    
    for (double dt=0; dt<400; dt+=8.33)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      if( tmpPk < tmpsgnl )
        tmpPk = tmpsgnl; // Storing in mV
      tmpCh += int(tmpsgnl/fadc2mvUub); // Sotoring in FADC
    }
    tmpPk /= fadc2mvUub; // Converting to FADC
    histPkUub->Fill( int(tmpPk) ); // Storing in FADC
    histChUub->Fill( tmpCh ); // Storing in FADC

    tmppk1ns = 0;
    tmpch1ns = 0;    
    for (double dt=0; dt<400; dt+=12.5)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      pulse[pl]->SetBinContent( int(dt), tmpsgnl );
      if( tmppk1ns < tmpsgnl )
        tmppk1ns = tmpsgnl; // Storing in mV
      tmpch1ns += tmpsgnl; // Sotoring in FADC
    }
    histPk->Fill( tmppk1ns ); // Storing in FADC
    histCh->Fill( tmpch1ns ); // Storing in FADC
  }

  TPaveStats *ptstats;

  TCanvas *c1 = new TCanvas("c1","c1",1600, 900);
  c1->cd();

  histPk->SetLineColor(kBlack);
  histPk->SetLineWidth(2);
  histPk->SetFillColor(kBlack);
  histPk->SetFillStyle(3001);
  histPk->SetTitle("");
  histPk->GetXaxis()->SetTitle("Peak [FADC]");
  histPk->GetXaxis()->SetRangeUser(20, 170);
  histPk->GetYaxis()->SetTitle("Counts [au]");
  histPk->GetYaxis()->SetRangeUser(0, 2350);
  histPk->Draw("");

  histPkUb->SetLineColor(kBlue);
  histPkUb->SetLineWidth(2);
  histPkUb->SetFillColor(kBlue);
  histPkUb->SetFillStyle(3001);
  histPkUb->Draw("SAMES");

  histPkUub->SetLineColor(kRed);
  histPkUub->SetLineWidth(2);
  histPkUub->SetFillColor(kRed);
  histPkUub->SetFillStyle(3001);
  histPkUub->Draw("SAMES");

  ptstats = new TPaveStats(0.68,0.75,0.89,0.89,"brNDC");
  histPk->SetName("Ideal (1.0 ns)");
  ptstats->SetParent(histPk);
  ptstats->SetTextColor(kBlack);
  histPk->GetListOfFunctions()->Add(ptstats);
  cout << histPk->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.58,0.89,0.72,"brNDC");
  histPkUb->SetName("UB: FADC");
  ptstats->SetParent(histPkUb);
  ptstats->SetTextColor(kBlue);
  histPkUb->GetListOfFunctions()->Add(ptstats);
  cout << histPkUb->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.4,0.89,0.54,"brNDC");
  histPkUub->SetName("UUB: FADC");
  ptstats->SetTextColor(kRed);
  histPkUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histPkUub);
  cout << histPkUub->GetMean() << endl;

	c1->Print("../plots/samplingDigiPkHistosFadc.pdf");


  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 900);
  c2->cd();
  gPad->SetLogy();

  histCh->SetLineColor(kBlack);
  histCh->SetLineWidth(2);
  histCh->SetFillColor(kBlack);
  histCh->SetFillStyle(3001);
  histCh->SetTitle("");
  histCh->GetXaxis()->SetTitle("Charge [FADC]");
  histCh->GetXaxis()->SetRangeUser(50, 6500);
  histCh->GetYaxis()->SetTitle("Counts [au]");
  histCh->GetYaxis()->SetRangeUser(1, 750);
  histCh->Draw("");

  histChUb->SetLineColor(kBlue);
  histChUb->SetFillColor(kBlue);
  histChUb->SetFillStyle(3001);
  histChUb->SetLineWidth(2);
  histChUb->Draw("SAMES");

  histChUub->SetLineColor(kRed);
  histChUub->SetFillColor(kRed);
  histChUub->SetFillStyle(3001);
  histChUub->SetLineWidth(2);
  histChUub->Draw("SAMES");

  ptstats = new TPaveStats(0.38,0.71,0.59,0.89, "brNDC"); //0.67,0.71,0.89,0.89,"brNDC");
  histCh->SetName("1.0 ns sampling");
  ptstats->SetParent(histCh);
  ptstats->SetTextColor(kBlack);
  histCh->GetListOfFunctions()->Add(ptstats);
  cout << histCh->GetMean() << endl;

  ptstats = new TPaveStats(0.38,0.48,0.59,0.66,"brNDC");
  histChUb->SetName("UB: FADC)");
  ptstats->SetParent(histChUb);
  ptstats->SetTextColor(kBlue);
  histChUb->GetListOfFunctions()->Add(ptstats);
  cout << histChUb->GetMean() << endl;

  ptstats = new TPaveStats(0.38,0.25,0.59,0.43,"brNDC");
  histChUub->SetName("UUB: FADC"); //ns/50 #Omega)(0.49 mV)");
  ptstats->SetTextColor(kRed);
  histChUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histChUub);
  cout << histChUub->GetMean() << endl;

	c2->Print("../plots/samplingDigiChHistosFadc.pdf");

  f->Write();
}
