void readDeltasAscii() {
    auto dataTree = new TTree("dataTree", "dataTree");

    dataTree->ReadFile("../deltas_python.dat", "gps/D:stId/I:pmtId/I:cQpk/D:cQpkErr/D:cQpkOff/D:cQpkErrOff/D:Qpk/D:QpkErr/D:QpkOff/D:QpkOffErr/D:delta/D:deltaErr/D:vh/D"); 

    //int nVals = dataTree->Draw("gps:stId:pmtId", "cQpkOff>0 && (cQpk-cQpkOff)/sqrt(cQpkErr**2+cQpkErrOff**2)>10.", "goff");
    //int nVals = dataTree->Draw("gps:stId:pmtId", "cQpkOff>0 && (cQpk-cQpkOff)/sqrt(cQpkErr**2+cQpkErrOff**2)<-10.", "goff");
    int nVals = dataTree->Draw("gps:stId:pmtId:delta", "pmtId==1 && 100.*delta<10.0 && vh<0.8", "goff");

    double *gpsVals = dataTree->GetVal(0);
    double *ids = dataTree->GetVal(1);
    double *pmts = dataTree->GetVal(2);
    double *deltas = dataTree->GetVal(3);

    //cout << "Ids > 10 %" << endl;
    //cout << "Ids < -10 %" << endl;
    auto distAwkward = new TH1D("distAwkward", "", 200, -5.0, 10.0);
    auto distAwkward3 = new TH1D("distAwkward3", "", 200, -5.0, 10.0);
    for( int i=0; i<nVals; i++ ) {
        //cout << ids[i] << " " << pmts[i] << " " << int(gpsVals[i]) << endl;
        distAwkward->Fill(100.*deltas[i]);
        //cout << ids[i] << ", ";
    }
    //dataTree->ReadFile("../deltas_python.dat", "gps/D:stId/I:pmtId/I:cQpk/D:cQpkErr/D:cQpkOff/D:cQpkErrOff/D:Qpk/D:QpkErr/D:QpkOff/D:QpkOffErr/D:delta/D:deltaErr/D:vh/D"); 
    int nVals3 = dataTree->Draw("gps:stId:pmtId:delta", "pmtId==3 && 100.*delta<10.0 && vh<0.8", "goff");
    double *deltas3 = dataTree->GetVal(3);
    for( int i=0; i<nVals3; i++ )
        distAwkward3->Fill(100.*deltas3[i]);
    cout << endl;
    //
    /*
    cout << "GPS for < -10 %" << endl;
    for( int i=0; i<nVals; i++)
        cout << int(gpsVals[i]) << ", ";
    cout << endl;
    */
    
    TCanvas *c = new TCanvas("c", "c", 1600, 900);
    c->cd();
    //dataTree->Draw("stId", "cQpkOff>0 && (cQpk-cQpkOff)/sqrt(cQpkErr**2+cQpkErrOff**2)>10.");
    //dataTree->Draw("(cQpk-cQpkOff)/sqrt(cQpkErr**2+cQpkErrOff**2)", "cQpkOff>0");
    distAwkward->Draw();
    distAwkward3->SetLineColor(kGreen+3);
    distAwkward3->Draw("same");
    cout << "MSD " << distAwkward->GetMean() << " " << distAwkward3->GetMean() << endl;
    cout << "MSD " << 200.*(distAwkward->GetMean() - distAwkward3->GetMean())/(distAwkward->GetMean() + distAwkward3->GetMean()) << endl;
}