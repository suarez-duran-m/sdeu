#include "Riostream.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"

Double_t ChiSquareDistr(Double_t *x, Double_t *par) {
  // Chisquare density distribution for nrFree degrees of freedom
  Double_t nrFree = par[0];
  Double_t chi2 = x[0];
  
  if (chi2 > 0) {
    Double_t lambda = nrFree/2.;
    Double_t norm = TMath::Gamma(lambda)*TMath::Power(2.,lambda);
    return TMath::Power(chi2,lambda-1)*TMath::Exp(-0.5*chi2)/norm;
  } 
  else
    return 0.0;
}

void chi2Dist(Int_t nrFree)
{
  TCanvas *c1 = new TCanvas("c1","");
  c1->cd();

  //TF1 *f1 = new TF1("chi-square distribution", ChiSquareDistr, 0.01, 10, 1);
  TF1 *f1 = new TF1("Chi2 distribution n=4", "x/(4*TMath::Exp(x/2.))", 0, 20);
  TF1 *f2 = new TF1("Chi2 distribution n=5", "x^(1.5)/(2^(2.5)*(2*1.7724)*TMath::Exp(x/2.))", 0, 20);

  //f1->SetParameter(0, Double_t(nrFree));
  f1->SetTitle("Chi2 distribution f(z,n)");
  f1->SetLineColor(49);
  f1->Draw();

  f1->GetHistogram()->GetXaxis()->SetTitle("z");
  f1->GetHistogram()->SetYTitle("f(z;n)");
  
  f2->SetTitle("Chi2 distribution n=5");
  f2->SetLineColor(40);
  f2->Draw("same");

  TLegend *leg = new TLegend(0.75,0.78,0.98,0.88);
  leg->AddEntry(f1, "n=4", "l");
  leg->AddEntry(f2, "n=5", "l");
  leg->SetTextSize(0.05);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->Draw();
  c1->Modified();
  c1->Print("../plots2/chi2Distribution.pdf");
  
  exit(0);
}
