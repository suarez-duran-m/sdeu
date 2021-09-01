#include "Riostream.h"
void samplingDigiPulsesAmpli() 
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("samplingDigiPulsesAmpli.C","");
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
  meanp0 *=1.5;
  sigmp0 = myfunc->GetParError(0);
  meanp1 = 9.55307e+01;
  sigmp1 = myfunc->GetParError(1);
  meanp2 = 1.5537e+01;
  sigmp2 = myfunc->GetParError(2);
   
  int Npulses = 1e4;
  double fadc2mvUb = 1.95;//2000./1024.; //1.95; // mV
  double fadc2mvUub = 0.49; //2000./4096.; //0.49; // mV

  double factorChUb = 0.5;//25./50.; //25. ns / 50. Ohms
  double factorChUub = 0.17;//8.33/50.; //8.33 ns / 50. Ohms

  TH1F **pulse = new TH1F *[Npulses];
  TString pulId;

  TH1I *histPkUb = new TH1I("histPkUb","UB", 200, 0, 200); // In mV
  TH1I *histChUb = new TH1I("histChUb","UB", 200, 0, 200); // In pC

  TH1I *histPkUub = new TH1I("histPkUub","UUB", 200, 0, 200); // In mV
  TH1I *histChUub = new TH1I("histChUub","UUB", 200, 0, 200); // In pC

  TH1I *histPk = new TH1I("histPk","1. ns", 200, 0, 200); // In mV
  TH1I *histCh = new TH1I("histCh","1. ns", 200, 0, 200); // In pC

  double tmpsgnl = 0.;
  double tmpPk = 0;
  int tmpCh = 0;

  double tmppk1ns = 0.;
  double tmpch1ns = 0.;
  double tmp = 0.;

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
      tmpCh += int(tmpsgnl / fadc2mvUb); // Storing in FADC
    }
    tmpPk /= fadc2mvUb; // Converting to FADC
    histPkUb->Fill( tmpPk*fadc2mvUb ); // Storing in mV
    tmp = 100.*(factorChUb*fadc2mvUb)+0.5; // to avoid strange bin-peak 
    histChUb->Fill( tmpCh*(int(tmp)/100.) );//(tmp/10.) ); //(factorChUb*fadc2mvUb) ); // Storing in pC

    tmpPk = 0;
    tmpCh = 0;    
    for (double dt=0; dt<400; dt+=8.33)
    {
      tmpsgnl = tmpp0*TMath::Landau( dt, tmpp1, tmpp2 );
      if( tmpPk < tmpsgnl )
        tmpPk = tmpsgnl;
      tmpCh += int(tmpsgnl / fadc2mvUub); // Sotoring in FADC
    }
    tmpPk /= fadc2mvUub; // Converting to FADC
    histPkUub->Fill( tmpPk*fadc2mvUub ); // Storing in mV
    histChUub->Fill( tmpCh*(factorChUub*fadc2mvUub) ); //(factorChUub*fadc2mvUub) ); // Storing in pC

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
    histCh->Fill( tmpch1ns/50. ); //tmpch1ns/50. ); // Storing in pC
  }

  TPaveStats *ptstats;
  
/*
  TCanvas *c0 = new TCanvas("c0","c0",1600, 900);
  TH1F *h0 = new TH1F("h0","h0", 16, 0, 400);
  TH1F *h1 = new TH1F("h1","h1", 16, 0, 400);

  c0->cd();
  pulse[0]->Draw();
  pulse[0]->SetStats(kFALSE);
  pulse[0]->SetLineColor(kBlack);
  pulse[0]->SetLineWidth(2);
  pulse[0]->SetTitle("");
  pulse[0]->GetXaxis()->SetTitle("Time [ns]");
  pulse[0]->GetYaxis()->SetTitle("Amplitude [mV]");
  pulse[0]->GetYaxis()->SetRangeUser(0, 116);
  pulse[0]->Draw();

  for ( int pl=0; pl<Npulses/10.; pl++ )
  {
    pulse[pl]->SetStats(kFALSE);
    pulse[pl]->SetLineColor(kBlack);
    pulse[pl]->SetLineWidth(2);
    pulse[pl]->Draw("SAMES");
  }
  c0->Print("../plots/artifitialPulsesAmp.pdf");
  */

  TCanvas *c1 = new TCanvas("c1","c1",1600, 900);

  c1->cd();
  histPk->SetLineColor(kBlack);
  histPk->SetLineWidth(2);
  histPk->SetFillColor(kBlack);
  histPk->SetFillStyle(3001);
  histPk->SetTitle("");
  histPk->GetXaxis()->SetTitle("Max. Amplitude [mV]");
  histPk->GetXaxis()->SetRangeUser(82, 122);
  //histPk->GetXaxis()->SetTitleSize(28);
  histPk->GetYaxis()->SetTitle("Counts [au]");
  histPk->GetYaxis()->SetRangeUser(0, 1400);
  histPk->Draw("");

  histPkUb->SetLineColor(kBlue);
  histPkUb->SetLineWidth(2);
  histPkUb->SetFillColor(kBlue);
  histPkUb->SetFillStyle(3001);
  histPkUb->GetXaxis()->SetTitleSize(28);
  histPkUb->Draw("SAMES");

  histPkUub->SetLineColor(kRed);
  histPkUub->SetLineWidth(2);
  histPkUub->SetFillColor(kRed);
  histPkUub->SetFillStyle(3001);
  histPkUub->GetXaxis()->SetTitleSize(28);
  histPkUub->Draw("SAMES");

  ptstats = new TPaveStats(0.68,0.72,0.89,0.89,"brNDC");
  histPk->SetName("1.0 ns sampling");
  ptstats->SetParent(histPk);
  ptstats->SetTextColor(kBlack);
  histPk->GetListOfFunctions()->Add(ptstats);
  cout << "========\nPrinting averages\n" << endl;
  cout << histPk->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.46,0.89,0.62,"brNDC");
  histPkUb->SetName("UB: 1FADC -> 1.95 mV");
  ptstats->SetParent(histPkUb);
  ptstats->SetTextColor(kBlue);
  histPkUb->GetListOfFunctions()->Add(ptstats);
  cout << histPkUb->GetMean() << endl;

  ptstats = new TPaveStats(0.68,0.20,0.89,0.36,"brNDC");
  histPkUub->SetName("UUB:\n 1FADC -> 0.49 mV");
  ptstats->SetTextColor(kRed);
  histPkUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histPkUub);
  cout << histPkUub->GetMean() << endl;

	c1->Print("../plots/samplingiDigiPkHistosAmp.pdf");


  TCanvas *c2 = new TCanvas("c2", "c2", 1600, 900);
  c2->cd();

  histCh->SetLineColor(kBlack);
  histCh->SetLineWidth(2);
  histCh->SetFillColor(kBlack);
  histCh->SetFillStyle(3001);
  histCh->SetTitle("");
  histCh->GetXaxis()->SetTitle("Charge [pC]");
  histCh->GetXaxis()->SetRangeUser(135, 220);
  histCh->GetYaxis()->SetTitle("Counts [au]");
  histCh->GetYaxis()->SetRangeUser(0, 550);
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

  ptstats = new TPaveStats(0.67,0.71,0.89,0.89,"brNDC");
  histCh->SetName("1.0 ns sampling");
  ptstats->SetParent(histCh);
  ptstats->SetTextColor(kBlack);
  histCh->GetListOfFunctions()->Add(ptstats);
  cout << histCh->GetMean() << endl;

  ptstats = new TPaveStats(0.67,0.48,0.89,0.68,"brNDC");
  histChUb->SetName("UB: 1FADC -> (25 ns/50 #Omega)(1.95 mV)");
  ptstats->SetParent(histChUb);
  ptstats->SetTextColor(kBlue);
  histChUb->GetListOfFunctions()->Add(ptstats);
  cout << histChUb->GetMean() << endl;

  ptstats = new TPaveStats(0.67,0.25,0.89,0.45,"brNDC");
  histChUub->SetName("UUB: 1FADC -> (8.33 ns/50 #Omega)(0.49 mV)");
  ptstats->SetTextColor(kRed);
  histChUub->GetListOfFunctions()->Add(ptstats);
  ptstats->SetParent(histChUub);
  cout << histChUub->GetMean() << endl;

	c2->Print("../plots/samplingDigiChHistosAmp.pdf");


  /*
  ofstream rmdpulsesAmpli;
  rmdpulsesAmpli.open("rmdpulsesAmpli.dat");
  for (int x = 0; x<1e4; x++)
  {
    tmpp0 = gRandom->Gaus(1.5*meanp0, sigmp0/2.); // Twice the original amplitud
    tmpp1 = gRandom->Gaus(meanp1, sigmp1/2.);
    tmpp2 = gRandom->Gaus(meanp2, sigmp2/2.);
    for (int tmpbin=0; tmpbin<400;tmpbin++)
      rmdpulsesAmpli << tmpbin << " " << tmpp0*TMath::Landau( tmpbin, tmpp1, tmpp2 ) << "\n";
  }

  rmdpulsesAmpli.close();
  */
  
  /*
  ofstream rmdpulsesMean;
  rmdpulsesMean.open("rmdpulsesMean.dat");
  for (int x = 0; x<1e4; x++)
  {
    tmpp0 = gRandom->Gaus(meanp0, sigmp0/2.); // Twice the original amplitud
    tmpp1 = gRandom->Gaus(1.5*meanp1, sigmp1/2.);
    tmpp2 = gRandom->Gaus(meanp2, sigmp2/2.);
    for (int tmpbin=0; tmpbin<400;tmpbin++)
      rmdpulsesMean << tmpbin << " " << tmpp0*TMath::Landau( tmpbin, tmpp1, tmpp2 ) << "\n";
  }

  rmdpulsesMean.close();
  */

/*
  ofstream rmdpulsesSig;
  rmdpulsesSig.open("rmdpulsesSig.dat");
  for (int x = 0; x<1e4; x++)
  {
    tmpp0 = gRandom->Gaus(meanp0, sigmp0/2.); // Twice the original amplitud
    tmpp1 = gRandom->Gaus(meanp1, sigmp1/2.);
    tmpp2 = gRandom->Gaus(1.5*meanp2, sigmp2/2.);
    for (int tmpbin=0; tmpbin<400;tmpbin++)
      rmdpulsesSig << tmpbin << " " << tmpp0*TMath::Landau( tmpbin, tmpp1, tmpp2 ) << "\n";
  }

  rmdpulsesSig.close();
  */

  // ===============================
  // *** Doing Area vs Amplitude ***

  /*
  ofstream rmdpulsesAmp;
  for (double fc = 1.1; fc<1.5; fc+=0.05 )
  {
    rmdpulsesAmp.open("rmdpulsesAmp"+to_string(fc)+".dat");
    for (int x = 0; x<1e4; x++)
    {
      tmpp0 = gRandom->Gaus(fc*meanp0, sigmp0/2.); // Twice the original amplitud
      tmpp1 = gRandom->Gaus(meanp1, sigmp1/2.);
      tmpp2 = gRandom->Gaus(meanp2, sigmp2/2.);
      for (int tmpbin=0; tmpbin<400; tmpbin++)
        rmdpulsesAmp << tmpbin << " " << tmpp0*TMath::Landau( tmpbin, tmpp1, tmpp2 ) << "\n";
    }
    rmdpulsesAmp.close();
  }
  */

  //tmp->Draw();


  //RooRealVar xx("x","x", 0,400);
  //RooRealVar mean("mean","mean of gaussian", 91.17,0,400);
  //RooRealVar sigma("sigma","width of gaussian",0,17.13,400);

  //chFit->Draw();
  //h->Draw();

  f->Write();
}
