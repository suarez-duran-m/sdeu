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

void fillQpkTimeVals(bool ifIsUub, int pmt, int st_id, 
    vector<double> &retQpkVect, vector<double> &retTimeVect) {
  TString bnCdas = (ifIsUub) ?
    "~/2021/sdeu/underHisto/results/uubChPkPMT" :
    "~/2021/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  TString monthUub[3] = {"Aug", "Sep", "Oct"};
  int nMonths = (ifIsUub) ?
    sizeof(monthUub)/sizeof(*monthUub) : 
    sizeof(monthUub)/sizeof(*monthUub)-1;
  
  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 0;
  int lstYear = (ifIsUub) ? nYears-1 : 1;
  
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<nMonths; month_i++ ) {
      fname = bnCdas + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      chargeInfo->SetBranchAddress("timeEvnt", &fetchTime);

      fetchQpkVals = 0.;
      fetchTime = 0;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals <= 0 )
          continue;
        retQpkVect.push_back( fetchQpkVals );
        retTimeVect.push_back( fetchTime );
      }
      f->Clear();
      f->Close();
    }
  }
}


void fetchingQpksPerStationVsTime(bool ifIsUub, int st_id) {

  TLegend *leg;
  TString strStId;
  strStId.Form("%d", st_id);
  TString strIfUub = (ifIsUub) ? "UUB" : "UB";

  vector < double > qpkValsPmt1;
  vector < double > evtTimePmt1;
  vector < double > qpkValsPmt2;
  vector < double > evtTimePmt2;
  vector < double > qpkValsPmt3;
  vector < double > evtTimePmt3;
  
  fillQpkTimeVals(ifIsUub, 1, st_id, qpkValsPmt1, evtTimePmt1);
  fillQpkTimeVals(ifIsUub, 2, st_id, qpkValsPmt2, evtTimePmt2);
  fillQpkTimeVals(ifIsUub, 3, st_id, qpkValsPmt3, evtTimePmt3);

  TGraph *grpPmt1 = new TGraph(qpkValsPmt1.size(), &evtTimePmt1[0], &qpkValsPmt1[0]);
  TGraph *grpPmt2 = new TGraph(qpkValsPmt2.size(), &evtTimePmt2[0], &qpkValsPmt2[0]);
  TGraph *grpPmt3 = new TGraph(qpkValsPmt3.size(), &evtTimePmt3[0], &qpkValsPmt3[0]);
  
  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  grpPmt1->SetTitle("");
  grpPmt1->GetXaxis()->SetTitle("Time since September 1st, 2021 [day/month]");
  grpPmt1->GetXaxis()->SetTimeFormat("%m/%d");
  grpPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpPmt1->GetYaxis()->SetRangeUser(900, 1650);
  grpPmt1->GetYaxis()->SetTitle("Q^{Pk}_{VEM} [FADC]");
  grpPmt1->SetMarkerStyle(72);
  grpPmt1->SetMarkerColor(kBlack);
  grpPmt1->SetMarkerSize(2);
  grpPmt1->Draw("AP");

  grpPmt2->SetMarkerStyle(73);
  grpPmt2->SetMarkerColor(kGreen+3);
  grpPmt2->SetMarkerSize(2);
  grpPmt2->Draw("P");

  grpPmt3->SetMarkerStyle(74);
  grpPmt3->SetMarkerColor(kBlue);
  grpPmt3->SetMarkerSize(2);
  grpPmt3->Draw("P");

  leg = new TLegend(0.82,0.75,0.95,0.95);
  leg->SetHeader("Station 827, "+strIfUub);
  leg->AddEntry(grpPmt1, "PMT1", "p");
  leg->AddEntry(grpPmt2, "PMT2", "p");
  leg->AddEntry(grpPmt3, "PMT3", "p");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();

  c1->Print("../plots/qpksVsTimeSt"+strStId+strIfUub+".pdf");
  
  exit(0);
} 
