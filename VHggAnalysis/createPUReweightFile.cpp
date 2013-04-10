#include <iostream>
#include "TFile.h"
#include "TH1F.h"


int main() {

  TFile* file_lept = TFile::Open("finalizedTrees_THq_feasibility/THq_tHqLeptonic_mH125_8TeV_testtest_presel_CSV.root ");
  TFile* file_hadr = TFile::Open("finalizedTrees_THq_feasibility/THq_tHqHadronic_mH125_8TeV_testtest_presel_CSV.root ");

  TH1F* h1_lept = (TH1F*)file_lept->Get("nvertex");
  TH1F* h1_hadr = (TH1F*)file_hadr->Get("nvertex");

  h1_lept->Scale(1./h1_lept->Integral());
  h1_hadr->Scale(1./h1_hadr->Integral());

  TFile* file_weights = TFile::Open("thq_hadronic_puweights.root", "recreate");
  TH1F* h1_weights = new TH1F(*h1_lept);
  h1_weights->SetName("puweights");
  h1_weights->Divide(h1_hadr);

  file_weights->cd();
  h1_weights->Write(); 
  file_weights->Close(); 

  return 0;

}
