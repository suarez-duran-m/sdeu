#include "Riostream.h"
void samplingDigiPulses() 
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("samplingDigiPulses.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%spulse_muon.dat",dir.Data()));

  Float_t x,y;
  TFile *f = new TFile("basic.root","RECREATE");
  TH1D *h = new TH1D("","",400, 0, 400);

  TCanvas *c0 = new TCanvas("c0","c0",1600, 900);
  TCanvas *c1 = new TCanvas("c1","c1",1600, 900);

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
  double fadc2mvUub = 2000./4096.;//0.49;//2000./4096.; // mV

  double factorChUb = 25./50.; // ns / Ohms
  double factorChUub = 8.33/50.; // ns / Ohms

  TH1F **pulse = new TH1F *[Npulses];
  TString pulId;

  TH1F *histPkUb = new TH1F("histPkUb","UB", 100, 0, 100); // In mV
  TH1F *histChUb = new TH1F("histChUb","UB", 300, 0, 300); // In pC

  TH1F *histPkUub = new TH1F("histPkUub","UUB", 100, 0, 100); // In mV
  TH1F *histChUub = new TH1F("histChUub","UUB", 300, 0, 300); // In pC

  TH1F *histPk = new TH1F("histPk","1. ns", 100, 0, 100); // In mV
  TH1F *histCh = new TH1F("histCh","1. ns", 300, 0, 300); // In pC

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
    histPkUb->Fill( tmpPk*fadc2mvUb ); // Storing in mV
    histChUb->Fill( tmpCh*(factorChUb*fadc2mvUb) ); // Sotoring in pC

    tmpPk = 0;
    tmpCh = 0;    
    for (double dt=0; dt<400; dt+=8.33)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      if( tmpPk < tmpsgnl )
        tmpPk = tmpsgnl;
      tmpCh += int(tmpsgnl/fadc2mvUub); // Sotoring in FADC
    }
    tmpPk /= fadc2mvUub; // Converting to FADC
    histPkUub->Fill( tmpPk*fadc2mvUub ); // Storing in mV
    histChUub->Fill( tmpCh*(factorChUub*fadc2mvUub) ); // Storing in pC

    tmppk1ns = 0;
    tmpch1ns = 0;    
    for (double dt=0; dt<400; dt+=1.)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      pulse[pl]->SetBinContent( int(dt), tmpsgnl );
      if( tmppk1ns < tmpsgnl )
        tmppk1ns = tmpsgnl;
      tmpch1ns += tmpsgnl; // Sotoring in FADC
    }
    histPk->Fill( tmppk1ns ); // Storing in mV
    histCh->Fill( tmpch1ns/50. ); // Storing in pC
  }

  TPaveStats *ptstats;

  c0->cd();
  pulse[0]->Draw();
  pulse[0]->SetStats(kFALSE);
  pulse[0]->SetLineColor(kBlack);
  pulse[0]->SetLineWidth(2);
  pulse[0]->SetTitle("");
  pulse[0]->GetXaxis()->SetTitle("Time [ns]");
  pulse[0]->GetYaxis()->SetTitle("Amplitude [mV]");
  pulse[0]->GetYaxis()->SetRangeUser(0, 80);
  for ( int pl=1; pl<Npulses/100.; pl++ )
  {
    pulse[pl]->SetStats(kFALSE);
    pulse[pl]->SetLineColor(kBlack);
    pulse[pl]->SetLineWidth(2);
    pulse[pl]->Draw("SAMES");
  }
  c0->Print("../plots/artifitialPulses.pdf");


  c1->cd();
  histPk->SetLineColor(kBlack);
  histPk->SetLineWidth(2);
  histPk->SetFillColor(kBlack);
  histPk->SetFillStyle(3001);
  histPk->SetTitle("");
  histPk->GetXaxis()->SetTitle("Peak [mV]");
  histPk->GetXaxis()->SetRangeUser(52, 85);
  histPk->GetYaxis()->SetTitle("Counts [au]");
  histPk->GetYaxis()->SetRangeUser(0, 1400);
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

  ptstats = new TPaveStats(0.68,0.72,0.89,0.89,"brNDC");
  histPk->SetName("1.0 ns sampling");
  ptstats->SetParent(histPk);
  ptstats->SetTextColor(kBlack);
  histPk->GetListOfFunctions()->Add(ptstats);
  cout << histPk->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.46,0.89,0.62,"brNDC");
  histPkUb->SetName("UB: 1FADC -> 1.95 mV");
  ptstats->SetParent(histPkUb);
  ptstats->SetTextColor(kBlue);
  histPkUb->GetListOfFunctions()->Add(ptstats);
  cout << histPkUb->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.20,0.89,0.36,"brNDC");
  histPkUub->SetName("UUB: 1FADC -> 0.49 mV");
  ptstats->SetTextColor(kRed);
  histPkUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histPkUub);
  cout << histPkUub->GetMean() << endl;

	c1->Print("../plots/samplingDigiPkHistos.pdf");


  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 900);
  c2->cd();

  histCh->SetLineColor(kBlack);
  histCh->SetLineWidth(2);
  histCh->SetFillColor(kBlack);
  histCh->SetFillStyle(3001);
  histCh->SetTitle("");
  histCh->GetXaxis()->SetTitle("Charge [pC]");
  histCh->GetXaxis()->SetRangeUser(82, 144);
  histCh->GetYaxis()->SetTitle("Counts [au]");
  histCh->GetYaxis()->SetRangeUser(0, 750);
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

  ptstats = new TPaveStats(0.68,0.71,0.89,0.89, "brNDC"); //0.67,0.71,0.89,0.89,"brNDC");
  histCh->SetName("1.0 ns sampling");
  ptstats->SetParent(histCh);
  ptstats->SetTextColor(kBlack);
  histCh->GetListOfFunctions()->Add(ptstats);
  cout << histCh->GetMean() << endl;

  ptstats = new TPaveStats(0.65,0.45,0.89,0.69,"brNDC");
  histChUb->SetName("UB: 1FADC -> (25 ns/50 #Omega)(1.95 mV)");
  ptstats->SetParent(histChUb);
  ptstats->SetTextColor(kBlue);
  histChUb->GetListOfFunctions()->Add(ptstats);
  cout << histChUb->GetMean() << endl;

  ptstats = new TPaveStats(0.65,0.22,0.89,0.43,"brNDC");
  histChUub->SetName("UUB: 1FADC -> (8.33 ns/50 #Omega)(0.49 mV)"); //ns/50 #Omega)(0.49 mV)");
  ptstats->SetTextColor(kRed);
  histChUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histChUub);
  cout << histChUub->GetMean() << endl;

	c2->Print("../plots/samplingDigiChHistos.pdf");

  f->Write();
}
