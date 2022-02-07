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
    vector<double> &retQpkVect, vector<double> &retTimeVect, 
    vector<double> &retChi2, vector<int> &retNdf) {
  TString bnCdas = (ifIsUub) ?
    "~/postdoc/sdeu/underHisto/results/uubChPkPMT" :
    "~/postdoc/sdeu/nouub/underHistos/results/ubChPkPMT"; 

  TString pmtId;
  TString strStId;
  TString fname;
  pmtId.Form("%d", pmt);
  strStId.Form("St%d", (int)st_id);

  TString monthUub[4] = {"Aug", "Sep", "Oct", "Nov"};
  int nMonths = sizeof(monthUub)/sizeof(*monthUub);
  
  int nYears = 5;
  TString strYear[5] = {"2016", "2018", "2019", "2020", "2021"};
  int stYear = (ifIsUub) ? nYears-1 : 1;
  int lstYear = (ifIsUub) ? nYears-1 : 1;
  
  TString strChargeData = "ChargeData";
  TFile *f;
  TTree *chargeInfo;
  double fetchQpkVals = 0.;
  int fetchTime = 0;
  double fetchChi2 = 0.;
  int fetchNdf = 0;

  for ( int year=stYear; year<=lstYear; year++ ) {
    for ( int month_i=0; month_i<nMonths; month_i++ ) {
      fname = bnCdas + pmtId + strStId + "lrb35" + monthUub[month_i] + strYear[year];
      f = TFile::Open(fname+".root");
      chargeInfo = (TTree*)f->Get(strChargeData);
      chargeInfo->SetBranchAddress("chargeVal", &fetchQpkVals);
      chargeInfo->SetBranchAddress("timeEvnt", &fetchTime);
      chargeInfo->SetBranchAddress("chi2", &fetchChi2);
      chargeInfo->SetBranchAddress("ndf", &fetchNdf);
      fetchQpkVals = 0.;
      fetchTime = 0;
      for( int etry=0; etry<chargeInfo->GetEntries(); etry++) {
        chargeInfo->GetEntry(etry);
        if ( fetchQpkVals <= 0 )
          continue;
        retQpkVect.push_back( fetchQpkVals );
        retTimeVect.push_back( fetchTime );
        retChi2.push_back( fetchChi2 );
        retNdf.push_back( fetchNdf );
        //if ( pmt == 1 )
          //cout << fetchTime << " " << fetchQpkVals << endl;
      }
      f->Clear();
      f->Close();
    }
  }
}

void doAvePerHour(bool ifIsUub, vector<double> qpkVect, vector<double> timeVect,
    vector<double> &retAveQpk, vector<double> &retRmsQpk, vector<double> &retDayQpk) {

  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;//1248652818; //1217116818; // August 1st
  int oneDay = 86400;
  double aveDay = 0.;
  double qpk2 = 0.;
  double rms = 0.;
  int qpkInDay = 0;
  int timeDiff = 0;
 
  for ( int qpk_i=0; qpk_i < qpkVect.size(); qpk_i++ ) {
    timeDiff = timeVect[qpk_i] - currentDay;
    if ( int(timeVect[qpk_i]) > 1311839895 && int(timeVect[qpk_i]) < 1312422592 )
      cout << "MSD " << int(timeVect[qpk_i]) << " " << qpkVect[qpk_i] << endl;
    //cout << (int)timeVect[qpk_i] << " " << qpkVect[qpk_i] << endl;
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= qpkInDay;
        if ( qpkInDay > 1 ) {
          retAveQpk.push_back( aveDay );
          rms = sqrt(qpk2/qpkInDay - aveDay*aveDay);
          // Using the error of the mean
          retRmsQpk.push_back( rms/sqrt(qpkInDay) );//rms/sqrt(qpkInDay) ); //sqrt(qpk2/qpkInDay - aveDay*aveDay) );
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
    }
    if ( qpkVect[qpk_i] > 0. ) {
      aveDay += qpkVect[qpk_i];
      qpk2 += qpkVect[qpk_i]*qpkVect[qpk_i];
      qpkInDay++;
    }
    if ( timeDiff > 2*oneDay )
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
  }
}

void doAvePerHour(bool ifIsUub, vector<double> chi2Vect, vector<int> ndfVect,
    vector<double> timeVect, vector<double> &retAveChi2Ndf, 
    vector<double> &retRmsChi2Ndf) { // Return the average Chi2Ndf as function of time

  int currentDay = (ifIsUub) ? 1311811218 : 1217116818;//1248652818; //1217116818; // August 1st
  int oneDay = 86400;
  double aveDay = 0.;
  double aveDay2 = 0.;
  int chi2InDay = 0;
  int timeDiff = 0;
  double rms = 0.;
 
  for ( int chi2_i=0; chi2_i < chi2Vect.size(); chi2_i++ ) {
    timeDiff = timeVect[chi2_i] - currentDay;
    if ( timeDiff > oneDay ) {
      if ( aveDay > 0 ) {
        aveDay /= chi2InDay;
        retAveChi2Ndf.push_back( aveDay );
        rms = sqrt(aveDay2/chi2InDay - aveDay*aveDay);
        // Using the error of the mean
        retRmsChi2Ndf.push_back( rms/sqrt(chi2InDay) ); //sqrt(aveDay2/chi2InDay - aveDay*aveDay) );
      }
      aveDay = 0.;
      aveDay2 = 0.;
      chi2InDay = 0;      
      currentDay += oneDay;
    }
    aveDay += chi2Vect[chi2_i]/ndfVect[chi2_i];
    aveDay2 += (chi2Vect[chi2_i]/ndfVect[chi2_i])*(chi2Vect[chi2_i]/ndfVect[chi2_i]);
    chi2InDay++;
    
    if ( timeDiff > 2*oneDay )
      currentDay += oneDay*(timeDiff/oneDay)-oneDay;
  }
}

void fetchingQpksPerStationVsTime(bool ifIsUub, int st_id, bool ifFit) {

  TLegend *leg;
  TString strStId;
  strStId.Form("%d", st_id);
  TString strIfUub = (ifIsUub) ? "UUB" : "UB";
  TString strChi2;
  TString strNdf;

  vector < double > qpkValsPmt1;
  vector < double > evtTimePmt1;
  vector < double > qpkValsPmt2;
  vector < double > evtTimePmt2;
  vector < double > qpkValsPmt3;
  vector < double > evtTimePmt3;

  vector < double > chi2Pmt1;
  vector < int > ndfPmt1;
  vector < double > chi2Pmt2;
  vector < int > ndfPmt2;
  vector < double > chi2Pmt3;
  vector < int > ndfPmt3;

  vector < double > aveDayPmt1;
  vector < double > rmsDayPmt1;
  vector < double > aveChi2NdfPmt1;
  vector < double > rmsChi2NdfPmt1;
  vector < double > timeAvePmt1;

  vector < double > aveDayPmt2;
  vector < double > rmsDayPmt2;
  vector < double > timeAvePmt2;
  vector < double > aveChi2NdfPmt2;
  vector < double > rmsChi2NdfPmt2;

  vector < double > aveDayPmt3;
  vector < double > rmsDayPmt3;
  vector < double > timeAvePmt3;
  vector < double > aveChi2NdfPmt3;
  vector < double > rmsChi2NdfPmt3;
  
  fillQpkTimeVals(ifIsUub, 1, st_id, qpkValsPmt1, evtTimePmt1, chi2Pmt1, ndfPmt1);
  doAvePerHour(ifIsUub, qpkValsPmt1, evtTimePmt1, aveDayPmt1, rmsDayPmt1, timeAvePmt1);
  //cout << endl << "MSD" << endl;
  doAvePerHour(ifIsUub, chi2Pmt1, ndfPmt1, evtTimePmt1, aveChi2NdfPmt1, rmsChi2NdfPmt1);

  fillQpkTimeVals(ifIsUub, 2, st_id, qpkValsPmt2, evtTimePmt2, chi2Pmt2, ndfPmt2);
  doAvePerHour(ifIsUub, qpkValsPmt2, evtTimePmt2, aveDayPmt2, rmsDayPmt2, timeAvePmt2);
  doAvePerHour(ifIsUub, chi2Pmt2, ndfPmt2, evtTimePmt2, aveChi2NdfPmt2, rmsChi2NdfPmt2);

  fillQpkTimeVals(ifIsUub, 3, st_id, qpkValsPmt3, evtTimePmt3, chi2Pmt3, ndfPmt3);
  doAvePerHour(ifIsUub, qpkValsPmt3, evtTimePmt3, aveDayPmt3, rmsDayPmt3, timeAvePmt3);
  doAvePerHour(ifIsUub, chi2Pmt3, ndfPmt3, evtTimePmt3, aveChi2NdfPmt3, rmsChi2NdfPmt3);

  TGraphErrors *grpPmt1 = new TGraphErrors(aveDayPmt1.size(), &timeAvePmt1[0], &aveDayPmt1[0], 0, &rmsDayPmt1[0]);
  TGraphErrors *grpPmt2 = new TGraphErrors(aveDayPmt2.size(), &timeAvePmt2[0], &aveDayPmt2[0], 0, &rmsDayPmt2[0]);
  //TGraph *grpPmt1 = new TGraph(qpkValsPmt2.size(), &evtTimePmt2[0], &qpkValsPmt2[0]);
  TGraphErrors *grpPmt3 = new TGraphErrors(aveDayPmt3.size(), &timeAvePmt3[0], &aveDayPmt3[0], 0, &rmsDayPmt3[0]);

  TGraphErrors *grpChi2NdfPmt1
    = new TGraphErrors(aveChi2NdfPmt1.size(),&timeAvePmt1[0],&aveChi2NdfPmt1[0], 0, &rmsChi2NdfPmt1[0]);
 
  TGraphErrors *grpChi2NdfPmt2 
    = new TGraphErrors(aveChi2NdfPmt2.size(),&timeAvePmt2[0],&aveChi2NdfPmt2[0], 0, &rmsChi2NdfPmt2[0]);
  TGraphErrors *grpChi2NdfPmt3 
    = new TGraphErrors(aveChi2NdfPmt3.size(),&timeAvePmt3[0],&aveChi2NdfPmt3[0], 0, &rmsChi2NdfPmt3[0]);


  int xMinTime = (ifIsUub) ? 1311724818 : 1217030418; //1248566418; //1217030418;  
  int xMaxTime = (ifIsUub) ? 1322352018 : 1227657618; //1256515218; //1224979218;

  double yMin = 0.;
  double yMax = 0.;

  double yMinMaxPm1 = grpPmt1->GetYaxis()->GetXmax();
  double yMinMaxPm2 = grpPmt2->GetYaxis()->GetXmax();
  yMax = (yMinMaxPm1 > yMinMaxPm2) ? yMinMaxPm1 : yMinMaxPm2;
  double yMinMaxPm3 = grpPmt3->GetYaxis()->GetXmax();
  yMax = (yMax > yMinMaxPm3) ? yMax : yMinMaxPm3;
  
  yMinMaxPm1 = grpPmt1->GetYaxis()->GetXmin();
  yMinMaxPm2 = grpPmt2->GetYaxis()->GetXmin();
  yMin = (yMinMaxPm1 < yMinMaxPm2) ? yMinMaxPm1 : yMinMaxPm2;
  yMinMaxPm3 = grpPmt3->GetYaxis()->GetXmin();
  yMin = (yMin < yMinMaxPm3) ? yMin : yMinMaxPm3;

  int binMinFit = 0;
  int binMaxFit = 0;
  
  if ( ifFit ) {
    binMinFit = 1316333394;
    binMaxFit = 1317022866;
    grpPmt1->Fit("pol0","", "", binMinFit,binMaxFit);
    binMinFit = 1320273234;
    binMaxFit = 1321356690;
    grpPmt2->Fit("pol0","", "", binMinFit,binMaxFit);
    binMinFit = 1112519698;
    binMaxFit = 1114116370;
    grpPmt3->Fit("pol0","", "", binMinFit,binMaxFit);
  } 
  
  //cout << "Chi2: " << grpPmt1->GetFunction("pol0")->GetChisquare()
    //<< " Chi2/ndf: " << grpPmt1->GetFunction("pol0")->GetChisquare() / grpPmt1->GetFunction("pol0")->GetNDF()
    //<< endl;
 
  /*
  TCanvas *c0 = canvasStyle("c0");
  c0->cd();

  grpChi2NdfPmt1->SetTitle("");
  grpChi2NdfPmt1->GetXaxis()->SetTitle("Time [Year/Day/Month]");
  grpChi2NdfPmt1->GetXaxis()->SetTimeFormat("%y/%m/%d");
  grpChi2NdfPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpChi2NdfPmt1->GetXaxis()->SetRangeUser(xMinTime, xMaxTime);
  grpChi2NdfPmt1->GetYaxis()->SetTitle("Chi2/ndf");
  grpChi2NdfPmt1->SetMarkerStyle(72);
  grpChi2NdfPmt1->SetMarkerColor(kRed);
  grpChi2NdfPmt1->SetLineColor(kRed);
  grpChi2NdfPmt1->SetMarkerSize(1);
  grpChi2NdfPmt1->Draw("AP");

  grpChi2NdfPmt2->SetMarkerStyle(73);
  grpChi2NdfPmt2->SetMarkerColor(kBlue);
  grpChi2NdfPmt2->SetMarkerSize(1);
  grpChi2NdfPmt2->SetLineColor(kBlue);
  grpChi2NdfPmt2->Draw("P");
  
  grpChi2NdfPmt3->SetMarkerStyle(74);
  grpChi2NdfPmt3->SetMarkerColor(kGreen+3);
  grpChi2NdfPmt3->SetMarkerSize(1);
  grpChi2NdfPmt3->Draw("P");

  leg = new TLegend(0.82,0.75,0.95,0.95);
  leg->SetHeader("Station "+strStId+", "+strIfUub);
  leg->AddEntry(grpPmt1, "PMT1", "p");
  leg->AddEntry(grpPmt2, "PMT2", "p");
  leg->AddEntry(grpPmt3, "PMT3", "p");
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0); 
  leg->SetFillStyle(0);
  leg->Draw();
  //c0->Print("../plots/chi2NdfQpksVsTimeSt"+strStId+strIfUub+".pdf");
  */


  TCanvas *c1 = canvasStyle("c1");
  c1->cd();

  grpPmt1->SetTitle("");
  grpPmt1->GetXaxis()->SetTitle("Time [Year/Day/Month]");
  grpPmt1->GetXaxis()->SetTimeFormat("%y/%m/%d");
  grpPmt1->GetXaxis()->SetTimeOffset(315964782,"gmt");
  grpPmt1->GetXaxis()->SetRangeUser(xMinTime, xMaxTime);
  //grpPmt1->GetYaxis()->SetRangeUser(1200, 1700);
  grpPmt1->GetYaxis()->SetRangeUser(yMin, yMax);
  grpPmt1->GetYaxis()->SetTitle("Q^{Pk}_{VEM} [FADC]");
  grpPmt1->SetMarkerStyle(72);
  grpPmt1->SetMarkerColor(kRed);
  grpPmt1->SetLineColor(kRed);
  grpPmt1->SetMarkerSize(1);
  if ( ifFit )
    grpPmt1->GetFunction("pol0")->SetLineColor(kBlack);
  grpPmt1->Draw("AP");
  
  grpPmt2->SetMarkerStyle(73);
  grpPmt2->SetMarkerColor(kBlue);
  grpPmt2->SetLineColor(kBlue);
  grpPmt2->SetMarkerSize(1);
  if ( ifFit )
    grpPmt2->GetFunction("pol0")->SetLineColor(kBlack);
  grpPmt2->Draw("P");
  
  grpPmt3->SetMarkerStyle(74);
  grpPmt3->SetMarkerColor(kGreen+3);
  grpPmt3->SetLineColor(kGreen+3);
  grpPmt3->SetMarkerSize(1);
  //if ( ifFit )
    //grpPmt3->GetFunction("pol0")->SetLineColor(kBlack);
  grpPmt3->Draw("P");

  leg = new TLegend(0.75,0.7,0.95,0.95);
  leg->SetHeader("Station "+strStId+", "+strIfUub+" Selected");
  if ( ifFit ) {
    strChi2.Form("%.2f",grpPmt1->GetFunction("pol0")->GetChisquare());
    strNdf.Form("%d",grpPmt1->GetFunction("pol0")->GetNDF());
    leg->AddEntry(grpPmt1, "PMT1", "p");
    leg->AddEntry(grpPmt1, "Chi2/ndf: "+strChi2+"/"+strNdf, "");

    strChi2.Form("%.2f",grpPmt2->GetFunction("pol0")->GetChisquare());
    strNdf.Form("%d",grpPmt2->GetFunction("pol0")->GetNDF());
    leg->AddEntry(grpPmt2, "PMT2", "p");
    leg->AddEntry(grpPmt2, "Chi2/ndf: "+strChi2+"/"+strNdf, "");

    //strChi2.Form("%.2f",grpPmt3->GetFunction("pol0")->GetChisquare());
    //strNdf.Form("%d",grpPmt3->GetFunction("pol0")->GetNDF());
    leg->AddEntry(grpPmt3, "PMT3", "p");
    //leg->AddEntry(grpPmt3, "Chi2/ndf: "+strChi2+"/"+strNdf, "");
  }

  leg->SetTextSize(0.03);
  leg->SetBorderSize(0); 
  leg->SetFillStyle(0);
  leg->Draw();

  //c1->Print("../plots/qpksVsTimeHandSt"+strStId+strIfUub+".pdf"); 
  //exit(0);
} 
