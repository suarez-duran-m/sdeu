void test() {
  auto fnew = TFile::Open("new.root", "recreate");
  fnew->cd();
  auto newtree = new TTree("newtree", "newtree");
  double gps;
  newtree->Branch("gps", &gps, "gps/D");
  
  auto f = TFile::Open("fittedHisto_delta_1344112041_1726_3.root"); 
  f->cd();
  auto tree = (TTree*)f->Get("T");

  f->cd();
  tree->SetBranchAddress("gpsTime", &gps); 
  tree->GetEntry(0); 

  fnew->cd();
  newtree->Fill(); 
  newtree->Write();
  
  f->cd();
  f->Close(); 
  fnew->cd();
  fnew->Write(); 
  fnew->Close();
}
