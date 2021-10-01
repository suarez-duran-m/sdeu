TCanvas *canvasStyle(TString name)
{
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

void histoStyle(TH1F *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.9);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void histoStyle(TGraphErrors *hist)
{
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.1);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitleSize(0.05);
}


void plottingAoPubUub(int st)
{
  TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
  dir.ReplaceAll("plottingAoPubUub.C","");
  dir.ReplaceAll("/./","/");
  TString stid;
  stid.Form("%d", st);

  // =====================================
  // *** *** *** Reading files *** *** ***

  // ============================
  // *** *** Open for UUB *** ***

  ifstream in;
  ifstream in2;
  in.open(Form("%s../../underHisto/uubaoptime"+stid+".dat", dir.Data()));
  in2.open(Form("%s../../nouub/underHistos/ubaoptime"+stid+".dat", dir.Data()));

  int nlinesub = 121; 
  int nlinesuub = 169; // 100; for 1740 // 169; for 863
  int nlines = nlinesub + nlinesuub;
  double tmpTime = 0;
  double tmpYaop = 0;
  double tmpErr = 0;
  double pmt = 0;
  int tmpcnt = 0;

  double uubtimePmt1[nlines];
  double uubaopPmt1[nlines];
  double uuberrPmt1[nlines];

  double uubtimePmt2[nlines];
  double uubaopPmt2[nlines];
  double uuberrPmt2[nlines];

  double uubtimePmt3[nlines];
  double uubaopPmt3[nlines];
  double uuberrPmt3[nlines];

  double ubtimePmt1[nlines];
  double ubaopPmt1[nlines];
  double uberrPmt1[nlines];

  double ubtimePmt2[nlines];
  double ubaopPmt2[nlines];
  double uberrPmt2[nlines];

  double ubtimePmt3[nlines];
  double ubaopPmt3[nlines];
  double uberrPmt3[nlines];


  while ( tmpcnt < 3*nlines ) 
  {
    if ( tmpcnt < 3*nlinesub )
    {
      in2 >> pmt >> tmpTime >> tmpYaop >> tmpErr;
      if ( pmt==1 )
      {
        ubtimePmt1[tmpcnt] = tmpTime;
        ubaopPmt1[tmpcnt] = tmpYaop*.5; // 25 ns / 50 ohm = .5 nF
        uberrPmt1[tmpcnt] = tmpErr;

        uubtimePmt1[tmpcnt] = tmpTime;
        uubaopPmt1[tmpcnt] = 0; 
        uuberrPmt1[tmpcnt] = 0; 
      }
      
      else if ( pmt==2 )
      {
        ubtimePmt2[tmpcnt-nlinesub] = tmpTime;
        ubaopPmt2[tmpcnt-nlinesub] = tmpYaop*.5;
        uberrPmt2[tmpcnt-nlinesub] = tmpErr;

        uubtimePmt2[tmpcnt-nlinesub] = tmpTime;
        uubaopPmt2[tmpcnt-nlinesub] = 0; 
        uuberrPmt2[tmpcnt-nlinesub] = 0; 
      }
      else if ( pmt==3 )
      {
        ubtimePmt3[tmpcnt-nlinesub*2] = tmpTime;
        ubaopPmt3[tmpcnt-nlinesub*2] = tmpYaop*.5;
        uberrPmt3[tmpcnt-nlinesub*2] = tmpErr;

        uubtimePmt3[tmpcnt-nlinesub*2] = tmpTime;
        uubaopPmt3[tmpcnt-nlinesub*2] = 0; 
        uuberrPmt3[tmpcnt-nlinesub*2] = 0; 
      }
    }
    else
    {
      in >> pmt >> tmpTime >> tmpYaop >> tmpErr;
      if ( pmt==1 )
      {
        uubtimePmt1[tmpcnt-2*nlinesub] = tmpTime;
        uubaopPmt1[tmpcnt-2*nlinesub] = tmpYaop*0.17; // 8.33 ns / 50 ohm = 0.17 nF
        uuberrPmt1[tmpcnt-2*nlinesub] = tmpErr;

        ubtimePmt1[tmpcnt-2*nlinesub] = tmpTime;
        ubaopPmt1[tmpcnt-2*nlinesub] = 0; 
        uberrPmt1[tmpcnt-2*nlinesub] = 0;
      }
      else if ( pmt==2 )
      {
        uubtimePmt2[tmpcnt-2*nlinesub-nlinesuub] = tmpTime;
        uubaopPmt2[tmpcnt-2*nlinesub-nlinesuub] = tmpYaop*.17;
        uuberrPmt2[tmpcnt-2*nlinesub-nlinesuub] = tmpErr;

        ubtimePmt2[tmpcnt-2*nlinesub-nlinesuub] = tmpTime;
        ubaopPmt2[tmpcnt-2*nlinesub-nlinesuub] = 0; 
        uberrPmt2[tmpcnt-2*nlinesub-nlinesuub] = 0;
      }
      else if ( pmt==3 )
      {
        uubtimePmt3[tmpcnt-2*nlinesub-nlinesuub*2] = tmpTime;
        uubaopPmt3[tmpcnt-2*nlinesub-nlinesuub*2] = tmpYaop*.17;
        uuberrPmt3[tmpcnt-2*nlinesub-nlinesuub*2] = tmpErr;

        ubtimePmt3[tmpcnt-2*nlinesub-nlinesuub*2] = tmpTime;
        ubaopPmt3[tmpcnt-2*nlinesub-nlinesuub*2] = 0; 
        uberrPmt3[tmpcnt-2*nlinesub-nlinesuub*2] = 0;
      }
    }
    tmpcnt++;
  }
  in.close();
  in2.close();

  TH1F *uubAoPpmt1 = new TH1F("uubAoPpmt1", "", nlines-1, uubtimePmt1);
  TH1F *uubAoPpmt2 = new TH1F("uubAoPpmt2", "", nlines-1, uubtimePmt2);
  TH1F *uubAoPpmt3 = new TH1F("uubAoPpmt3", "", nlines-1, uubtimePmt3);

  TH1F *ubAoPpmt1 = new TH1F("ubAoPpmt1", "", nlines-1, ubtimePmt1);
  TH1F *ubAoPpmt2 = new TH1F("ubAoPpmt2", "", nlines-1, ubtimePmt2);
  TH1F *ubAoPpmt3 = new TH1F("ubAoPpmt3", "", nlines-1, ubtimePmt3);


  for ( int kk=0; kk<nlines; kk++ )
  {
    uubAoPpmt1->SetBinContent( kk+1, uubaopPmt1[kk] );
    uubAoPpmt1->SetBinError( kk+1, uuberrPmt1[kk] );

    uubAoPpmt2->SetBinContent( kk+1, uubaopPmt2[kk] );
    uubAoPpmt2->SetBinError( kk+1, uuberrPmt2[kk] );

    uubAoPpmt3->SetBinContent( kk+1, uubaopPmt3[kk] );
    uubAoPpmt3->SetBinError( kk+1, uuberrPmt3[kk] );

    ubAoPpmt1->SetBinContent( kk, ubaopPmt1[kk] );
    ubAoPpmt1->SetBinError( kk, uberrPmt1[kk] );

    ubAoPpmt2->SetBinContent( kk, ubaopPmt2[kk] );
    ubAoPpmt2->SetBinError( kk, uberrPmt2[kk] );

    ubAoPpmt3->SetBinContent( kk, ubaopPmt3[kk] );
    ubAoPpmt3->SetBinError( kk, uberrPmt3[kk] );
  }

  // ====================================
  // *** *** *** Plotting AoP *** *** ***

  gStyle->SetOptStat(0);
  TLegend *leg;
  TCanvas *c1 = canvasStyle("c1");

  c1->cd();
  ubAoPpmt1->SetStats(0);
  ubAoPpmt1->SetMarkerColor(kCyan-2);
  ubAoPpmt1->SetLineColor(kCyan-2);
  ubAoPpmt1->GetXaxis()->SetTimeDisplay(1);
  ubAoPpmt1->GetXaxis()->SetTimeFormat("%d/%m");
  ubAoPpmt1->SetMarkerSize(1.5);
  ubAoPpmt1->SetMarkerStyle(8);
  ubAoPpmt1->GetYaxis()->SetRangeUser(1.1, 2.4);
  ubAoPpmt1->GetYaxis()->SetTitle("AoP [nF]");
  ubAoPpmt1->GetXaxis()->SetTitle("Time since August 1st, 2020 [day/month]");
  histoStyle(ubAoPpmt1);
  ubAoPpmt1->Draw("P");

  uubAoPpmt1->SetMarkerColor(kRed-3);
  uubAoPpmt1->SetLineColor(kRed-3);
  uubAoPpmt1->SetMarkerSize(1.5);
  uubAoPpmt1->SetMarkerStyle(8);
  uubAoPpmt1->Draw("P sames");

  leg = new TLegend(0.7, 0.66, 0.96, 0.94);
  leg->SetHeader("Station "+stid+" PMT1");
  leg->SetTextSize(0.06);
  leg->AddEntry(ubAoPpmt1, "UB ");
  leg->AddEntry(uubAoPpmt1, "UUB");
  leg->SetTextAlign(22);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->Draw();
  c1->Print("../../plots/uububAoPtimePmt1St"+stid+".pdf");

  TCanvas *c2 = canvasStyle("c2");

  c2->cd();
  ubAoPpmt2->SetStats(0);
  ubAoPpmt2->SetMarkerColor(kCyan-2);
  ubAoPpmt2->SetLineColor(kCyan-2);
  ubAoPpmt2->GetXaxis()->SetTimeDisplay(1);
  ubAoPpmt2->GetXaxis()->SetTimeFormat("%d/%m");
  ubAoPpmt2->SetMarkerSize(1.5);
  ubAoPpmt2->SetMarkerStyle(8);
  ubAoPpmt2->GetYaxis()->SetRangeUser(1.1, 2.4);
  ubAoPpmt2->GetYaxis()->SetTitle("AoP [nF]");
  ubAoPpmt2->GetXaxis()->SetTitle("Time since August 1st, 2020 [day/month]");
  histoStyle(ubAoPpmt2);
  ubAoPpmt2->Draw("P");

  uubAoPpmt2->SetMarkerColor(kRed-3);
  uubAoPpmt2->SetLineColor(kRed-3);
  uubAoPpmt2->SetMarkerSize(1.5);
  uubAoPpmt2->SetMarkerStyle(8);
  uubAoPpmt2->Draw("P sames");

  leg = new TLegend(0.7, 0.66, 0.96, 0.94);
  leg->SetHeader("Station "+stid+" PMT2");
  leg->SetTextSize(0.06);
  leg->AddEntry(ubAoPpmt2, "UB ");
  leg->AddEntry(uubAoPpmt2, "UUB");
  leg->SetTextAlign(22);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->Draw();
  c2->Print("../../plots/uububAoPtimePmt2St"+stid+".pdf");


  TCanvas *c3 = canvasStyle("c3");

  c3->cd();
  ubAoPpmt3->SetStats(0);
  ubAoPpmt3->SetMarkerColor(kCyan-2);
  ubAoPpmt3->SetLineColor(kCyan-2);
  ubAoPpmt3->GetXaxis()->SetTimeDisplay(1);
  ubAoPpmt3->GetXaxis()->SetTimeFormat("%d/%m");
  ubAoPpmt3->SetMarkerSize(1.5);
  ubAoPpmt3->SetMarkerStyle(8);
  ubAoPpmt3->GetYaxis()->SetRangeUser(1.1, 2.4);
  ubAoPpmt3->GetYaxis()->SetTitle("AoP [nF]");
  ubAoPpmt3->GetXaxis()->SetTitle("Time since August 1st, 2020 [day/month]");
  histoStyle(ubAoPpmt3);
  ubAoPpmt3->Draw("P");

  uubAoPpmt3->SetMarkerColor(kRed-3);
  uubAoPpmt3->SetLineColor(kRed-3);
  uubAoPpmt3->SetMarkerSize(1.5);
  uubAoPpmt3->SetMarkerStyle(8);
  uubAoPpmt3->Draw("P sames");

  leg = new TLegend(0.7, 0.66, 0.96, 0.94);
  leg->SetHeader("Station "+stid+" PMT3");
  leg->SetTextSize(0.06);
  leg->AddEntry(ubAoPpmt3, "UB ");
  leg->AddEntry(uubAoPpmt3, "UUB");
  leg->SetTextAlign(22);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);
  leg->Draw();

  c3->Print("../../plots/uububAoPtimePmt3St"+stid+".pdf");
}
