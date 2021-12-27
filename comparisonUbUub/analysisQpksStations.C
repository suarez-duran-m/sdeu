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

  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = (ifIsUub) ?
    sizeof(monthUub)/sizeof(*monthUub) : 
    sizeof(monthUub)/sizeof(*monthUub)-1;
  
  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
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

void doAvePerTime(bool ifIsUub, vector<double> qpkVect, vector<double> timeVect,
    vector<double> &retAveQpk, vector<double> &retRmsQpk, vector<double> &retDayQpk,
    vector<double> &retNQpkDay, TH1D *retDistQpkRel, TH1D *retDistErrMean ) {

  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;//1248652818; //1217116818; // August 1st
  int oneDay = 86400;
  double aveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;
  vector < double > qpks1day;
 
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    timeDiff = timeVect[qpk_i] - currentDay;
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= qpkInDay;
        if ( qpkInDay > 1 ) {
          retAveQpk.push_back( aveDay );
          rms = sqrt(qpk2/qpkInDay - aveDay*aveDay);
          // Using the error of the mean
          retRmsQpk.push_back( rms/sqrt(qpkInDay) ); //sqrt(qpk2/qpkInDay - aveDay*aveDay) );
          retDistErrMean->Fill( rms/sqrt(qpkInDay) );
          retNQpkDay.push_back( qpkInDay );
          // Filling Qpk/<Qpk>
          for ( auto & qpk : qpks1day )
            retDistQpkRel->Fill( qpk/aveDay );
        }
        else {
          retAveQpk.push_back( aveDay );
          retRmsQpk.push_back( 1.0 );
        }
        retDayQpk.push_back( currentDay+oneDay/2 );
      }
      aveDay = 0.;
      qpk2 = 0.;
      qpkInDay = 0;      
      currentDay += oneDay;
      qpks1day.clear();
    }
    if ( qpkVect[qpk_i] > 0. ) {
      aveDay += qpkVect[qpk_i];
      qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
      qpkInDay++;
      qpks1day.push_back( qpkVect[qpk_i] );
    }
    if ( timeDiff > 2*oneDay )
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
  }
}

double getMean(vector<double> vecValues) {
  double mean = 0;
  for( auto & val_i : vecValues )
    mean += val_i;

  return mean/vecValues.size();
}

double getRms(vector<double> vecValues) {
  double sumval2 = 0.;
  double sumval = 0.;

  for( auto & val_i : vecValues ) {
    sumval2 += (val_i*val_i);
    sumval += val_i;
  }
  sumval /= vecValues.size();
  sumval2 /= vecValues.size();
  return sqrt( sumval2 - sumval*sumval );
}


void analysisQpksStations() {

  //TFile *outFile = new TFile("analysisQpksStations.root", "RECREATE");
  TLegend *leg;
  TString strMean;
  TString strErr;

  // Reading list of stations
  string stListPath = "/home/msd/2021/sdeu/fullUubStationsListVert.txt";
  ifstream fileStList;
  double St_i;
  vector < double > stListId;
  fileStList.open(stListPath);
  while( fileStList.good() ) {
    fileStList >> St_i;
    stListId.push_back(St_i);
  }
  stListId.pop_back();
  fileStList.close();

  vector < vector < double > > qpkValsPmt1;
  vector < vector < double > > evtTimePmt1;
  vector < vector < double > > qpkValsPmt2;
  vector < vector < double > > evtTimePmt2;
  vector < vector < double > > qpkValsPmt3;
  vector < vector < double > > evtTimePmt3;

  qpkValsPmt1.resize(2);
  evtTimePmt1.resize(2);
  qpkValsPmt2.resize(2);
  evtTimePmt2.resize(2);
  qpkValsPmt3.resize(2);
  evtTimePmt3.resize(2);

  vector < vector < double > > aveDayPmt1;
  vector < vector < double > > rmsDayPmt1;
  vector < vector < double > > aveChi2NdfPmt1;
  vector < vector < double > > rmsChi2NdfPmt1;
  vector < vector < double > > timeAvePmt1;
  vector < vector < double > > nQpkDayPmt1;

  aveDayPmt1.resize(2);
  rmsDayPmt1.resize(2);
  aveChi2NdfPmt1.resize(2);
  rmsChi2NdfPmt1.resize(2);
  timeAvePmt1.resize(2);
  nQpkDayPmt1.resize(2);

  vector < vector < double > > aveDayPmt2;
  vector < vector < double > > rmsDayPmt2;
  vector < vector < double > > timeAvePmt2;
  vector < vector < double > > aveChi2NdfPmt2;
  vector < vector < double > > rmsChi2NdfPmt2;
  vector < vector < double > > nQpkDayPmt2;

  aveDayPmt2.resize(2);
  rmsDayPmt2.resize(2);
  timeAvePmt2.resize(2);
  aveChi2NdfPmt2.resize(2);
  rmsChi2NdfPmt2.resize(2);
  nQpkDayPmt2.resize(2);

  vector < vector < double > > aveDayPmt3;
  vector < vector < double > > rmsDayPmt3;
  vector < vector < double > > timeAvePmt3;
  vector < vector < double > > aveChi2NdfPmt3;
  vector < vector < double > > rmsChi2NdfPmt3;
  vector < vector < double > > nQpkDayPmt3;

  aveDayPmt3.resize(2);
  rmsDayPmt3.resize(2);
  timeAvePmt3.resize(2);
  aveChi2NdfPmt3.resize(2);
  rmsChi2NdfPmt3.resize(2);
  nQpkDayPmt3.resize(2);

  int nbins = 300;
  double frtBin = 0.;
  double lstBin = 3.;
  TH1D *distQpkRelPmt1 = new TH1D("distQpkRelPmt1", "", nbins, frtBin, lstBin);
  TH1D *distQpkRelPmt2 = new TH1D("distQpkRelPmt2", "", nbins, frtBin, lstBin);
  TH1D *distQpkRelPmt3 = new TH1D("distQpkRelPmt3", "", nbins, frtBin, lstBin);

  nbins = 150;
  frtBin = 0.;
  lstBin = 30.;
  TH1D *distErrMeanPmt1 = new TH1D("distErrMeanPmt1", "", nbins, frtBin, lstBin);
  TH1D *distErrMeanPmt2 = new TH1D("distErrMeanPmt2", "", nbins, frtBin, lstBin);
  TH1D *distErrMeanPmt3 = new TH1D("distErrMeanPmt3", "", nbins, frtBin, lstBin);

  nbins = 100;
  frtBin = -50;
  lstBin = 50;
  TH1D *distDiffAveQpkDayPmt1 = new TH1D("distDiffAveQpkDayPmt1", "", nbins, frtBin, lstBin);
  TH1D *distDiffAveQpkDayPmt2 = new TH1D("distDiffAveQpkDayPmt2", "", nbins, frtBin, lstBin);
  TH1D *distDiffAveQpkDayPmt3 = new TH1D("distDiffAveQpkDayPmt3", "", nbins, frtBin, lstBin);

  int st_id = 0;
  double aveNqpkDayUb = 0.;
  double aveNqpkDayUub = 0.;
  vector < double > diffAveQpkDayPmt1;
  vector < double > diffAveQpkDayPmt2;
  vector < double > diffAveQpkDayPmt3;
  vector < vector < double > > aveQpkPerDayPmt1;
  vector < vector < double > > aveQpkPerDayPmt2;
  vector < vector < double > > aveQpkPerDayPmt3;

  aveQpkPerDayPmt1.resize(2);
  aveQpkPerDayPmt2.resize(2);
  aveQpkPerDayPmt3.resize(2);

  for ( auto & st_i : stListId ) {
    st_id = (int)st_i;
  
    fillQpkTimeVals(false, 1, st_id, qpkValsPmt1[0], evtTimePmt1[0]);
    doAvePerTime(false, qpkValsPmt1[0], evtTimePmt1[0], aveDayPmt1[0], rmsDayPmt1[0], 
        timeAvePmt1[0], nQpkDayPmt1[0], distQpkRelPmt1, distErrMeanPmt1);

    fillQpkTimeVals(true, 1, st_id, qpkValsPmt1[1], evtTimePmt1[1]);
    doAvePerTime(true, qpkValsPmt1[1], evtTimePmt1[1], aveDayPmt1[1], rmsDayPmt1[1],
        timeAvePmt1[1], nQpkDayPmt1[1], distQpkRelPmt1, distErrMeanPmt1);
  
    fillQpkTimeVals(false, 2, st_id, qpkValsPmt2[0], evtTimePmt2[0]);
    doAvePerTime(false, qpkValsPmt2[0], evtTimePmt2[0], aveDayPmt2[0], rmsDayPmt2[0],
        timeAvePmt2[0], nQpkDayPmt2[0], distQpkRelPmt2, distErrMeanPmt2);

    fillQpkTimeVals(true, 2, st_id, qpkValsPmt2[1], evtTimePmt2[1]);
    doAvePerTime(true, qpkValsPmt2[1], evtTimePmt2[1], aveDayPmt2[1], rmsDayPmt2[1],
        timeAvePmt2[1], nQpkDayPmt2[1], distQpkRelPmt2, distErrMeanPmt2);

    fillQpkTimeVals(false, 3, st_id, qpkValsPmt3[0], evtTimePmt3[0]);
    doAvePerTime(false, qpkValsPmt3[0], evtTimePmt3[0], aveDayPmt3[0], rmsDayPmt3[0],
        timeAvePmt3[0], nQpkDayPmt3[0], distQpkRelPmt3, distErrMeanPmt3);

    fillQpkTimeVals(true, 3, st_id, qpkValsPmt3[1], evtTimePmt3[1]);
    doAvePerTime(true, qpkValsPmt3[1], evtTimePmt3[1], aveDayPmt3[1], rmsDayPmt3[1],
        timeAvePmt3[1], nQpkDayPmt3[1], distQpkRelPmt3, distErrMeanPmt3);

    aveNqpkDayUb = getMean(nQpkDayPmt1[0]);
    aveNqpkDayUub = getMean(nQpkDayPmt1[1]);
    if ( aveNqpkDayUb > 0 && aveNqpkDayUub > 0 ) {
      diffAveQpkDayPmt1.push_back( aveNqpkDayUb - aveNqpkDayUub );
      nQpkDayPmt1.clear();
      nQpkDayPmt1.clear();
      nQpkDayPmt1.resize(2);
      aveQpkPerDayPmt1[0].push_back( aveNqpkDayUb );
      aveQpkPerDayPmt1[1].push_back( aveNqpkDayUub );
      if ( fabs(aveNqpkDayUb - aveNqpkDayUub) > (3.96+5.20) )
        cout << st_id << " PMT1" << endl;
    }

    aveNqpkDayUb = getMean(nQpkDayPmt2[0]);
    aveNqpkDayUub = getMean(nQpkDayPmt2[1]);
    if ( aveNqpkDayUb > 0 && aveNqpkDayUub > 0 ) {
      diffAveQpkDayPmt2.push_back( aveNqpkDayUb - aveNqpkDayUub );
      nQpkDayPmt2.clear();
      nQpkDayPmt2.clear();
      nQpkDayPmt2.resize(2);
      aveQpkPerDayPmt2[0].push_back( aveNqpkDayUb );
      aveQpkPerDayPmt2[1].push_back( aveNqpkDayUub );
      if ( fabs(aveNqpkDayUb - aveNqpkDayUub) > (3.87+4.87) )
        cout << st_id << " PMT2" << endl;
    }

    aveNqpkDayUb = getMean(nQpkDayPmt3[0]);
    aveNqpkDayUub = getMean(nQpkDayPmt3[1]);
    if ( aveNqpkDayUb > 0 && aveNqpkDayUub > 0 ) {
      diffAveQpkDayPmt3.push_back( aveNqpkDayUb - aveNqpkDayUub );
      nQpkDayPmt3.clear();
      nQpkDayPmt3.clear();
      nQpkDayPmt3.resize(2);
      aveQpkPerDayPmt3[0].push_back( aveNqpkDayUb );
      aveQpkPerDayPmt3[1].push_back( aveNqpkDayUub );
      if ( fabs(aveNqpkDayUb - aveNqpkDayUub) > (4.02+4.92) )
        cout << st_id << " PMT3" << endl;
    }
    
    qpkValsPmt1.clear();
    qpkValsPmt1.resize(2);
    evtTimePmt1.clear();
    evtTimePmt1.resize(2);
    qpkValsPmt2.clear();
    qpkValsPmt2.resize(2);
    evtTimePmt2.clear();
    evtTimePmt2.resize(2);
    qpkValsPmt3.clear();
    qpkValsPmt3.resize(2);
    evtTimePmt3.clear();
    evtTimePmt3.resize(2);

    cout << "Finish for station " << st_id << endl;
    
    //if ( st_id == 830 )
      //break;
  }

  for ( auto val_i : diffAveQpkDayPmt1 )
    distDiffAveQpkDayPmt1->Fill( val_i );
  
  for ( auto val_i : diffAveQpkDayPmt2 )
    distDiffAveQpkDayPmt2->Fill( val_i );
 
  for ( auto val_i : diffAveQpkDayPmt3 )
    distDiffAveQpkDayPmt3->Fill( val_i );

  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  distDiffAveQpkDayPmt1->SetTitle("");
  distDiffAveQpkDayPmt1->SetStats(kFALSE);
  distDiffAveQpkDayPmt1->GetXaxis()->SetTitle("#LT N_{Q^{pk}_{VEM}} #GT^{UB} - #LT N_{Q^{pk}_{VEM}} #GT^{UUB}");
  distDiffAveQpkDayPmt1->GetYaxis()->SetTitle("Entries [ua]");
  distDiffAveQpkDayPmt1->SetMarkerStyle(72);
  distDiffAveQpkDayPmt1->SetMarkerColor(kRed);
  distDiffAveQpkDayPmt1->SetLineColor(kRed);
  distDiffAveQpkDayPmt1->SetMarkerSize(1);
  distDiffAveQpkDayPmt1->Draw();
  
  distDiffAveQpkDayPmt2->SetMarkerStyle(73);
  distDiffAveQpkDayPmt2->SetMarkerColor(kBlue);
  distDiffAveQpkDayPmt2->SetLineColor(kBlue);
  distDiffAveQpkDayPmt2->SetMarkerSize(1);
  distDiffAveQpkDayPmt2->Draw("same");
  
  distDiffAveQpkDayPmt3->SetMarkerStyle(74);
  distDiffAveQpkDayPmt3->SetMarkerColor(kGreen+3);
  distDiffAveQpkDayPmt3->SetLineColor(kGreen+3);
  distDiffAveQpkDayPmt3->SetMarkerSize(1);
  distDiffAveQpkDayPmt3->Draw("same");

  leg = new TLegend(0.65,0.4,0.95,0.95);
  //leg->SetHeader("Station "+strStId+", "+strIfUub+" Selected");
  strMean.Form("%.2f", distDiffAveQpkDayPmt1->GetMean());
  strErr.Form("%.2f", distDiffAveQpkDayPmt1->GetMeanError());
  leg->AddEntry(distDiffAveQpkDayPmt1, "PMT1", "l");
  leg->AddEntry(distDiffAveQpkDayPmt1, "#mu = "+strMean+" #pm "+strErr,"");
  strErr.Form("%.2f", distDiffAveQpkDayPmt1->GetRMS());
  leg->AddEntry(distDiffAveQpkDayPmt1, "RMS = "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt1[0]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt1[0])/sqrt(aveQpkPerDayPmt1.size()));
  leg->AddEntry(distDiffAveQpkDayPmt1, "Ave. Day UB: "+strMean+" #pm "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt1[1]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt1[1])/sqrt(aveQpkPerDayPmt1.size()));
  leg->AddEntry(distDiffAveQpkDayPmt1, "Ave. Day UUB: "+strMean+" #pm "+strErr,"");
  
  strMean.Form("%.2f", distDiffAveQpkDayPmt2->GetMean());
  strErr.Form("%.2f", distDiffAveQpkDayPmt2->GetMeanError());
  leg->AddEntry(distDiffAveQpkDayPmt2, "PMT2", "l");
  leg->AddEntry(distDiffAveQpkDayPmt2, "#mu = "+strMean+" #pm "+strErr,"");
  strErr.Form("%.2f", distDiffAveQpkDayPmt2->GetRMS());  
  leg->AddEntry(distDiffAveQpkDayPmt2, "RMS = "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt2[0]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt2[0])/sqrt(aveQpkPerDayPmt2.size()));
  leg->AddEntry(distDiffAveQpkDayPmt2, "Ave. Day UB: "+strMean+" #pm "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt2[1]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt2[1])/sqrt(aveQpkPerDayPmt2.size()));
  leg->AddEntry(distDiffAveQpkDayPmt2, "Ave. Day UUB: "+strMean+" #pm "+strErr,"");

  strMean.Form("%.2f", distDiffAveQpkDayPmt3->GetMean());
  strErr.Form("%.2f", distDiffAveQpkDayPmt3->GetMeanError());
  leg->AddEntry(distDiffAveQpkDayPmt3, "PMT3", "l");
  leg->AddEntry(distDiffAveQpkDayPmt3, "#mu = "+strMean+" #pm "+strErr,"");
  strErr.Form("%.2f", distDiffAveQpkDayPmt3->GetRMS());
  leg->AddEntry(distDiffAveQpkDayPmt3, "RMS = "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt3[0]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt3[0])/sqrt(aveQpkPerDayPmt3.size()));
  leg->AddEntry(distDiffAveQpkDayPmt3, "Ave. Day UB: "+strMean+" #pm "+strErr,"");
  strMean.Form("%.2f", getMean(aveQpkPerDayPmt3[1]));
  strErr.Form("%.2f", getRms(aveQpkPerDayPmt3[1])/sqrt(aveQpkPerDayPmt3.size()));
  leg->AddEntry(distDiffAveQpkDayPmt3, "Ave. Day UUB: "+strMean+" #pm "+strErr,"");

  leg->SetTextSize(0.03);
  leg->SetBorderSize(0); 
  leg->SetFillStyle(0);
  leg->Draw();
  //leg->SetName("legend");
  //leg->Write();

  c1->Print("../plots2/distDiffAveQpkDayPmts.pdf");

  //outFile->Write();
  //outFile->Close();

  //exit(0);
} 
