#include <iostream>
#include "CommonTools/DrawBase.h"
#include "TH1F.h"


void scanCut( DrawBase* db, TTree* tree_thq, TTree* tree_tth, const std::string& var );


int main() {

  DrawBase* db = new DrawBase("scanTTH");
  db->set_outputdir("prova");

  TFile* file_tth = TFile::Open("finalizedTrees_THq_feasibility/THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_sel4_CSV.root");
  TFile* file_thq = TFile::Open("finalizedTrees_THq_feasibility/THq_tHq_mH125_8TeV_testtest_sel4_CSV.root");

  TTree* tree_tth = (TTree*)file_tth->Get("tree_passedEvents");
  TTree* tree_thq = (TTree*)file_thq->Get("tree_passedEvents");

  scanCut( db, tree_thq, tree_tth, "BDT_lept");
  scanCut( db, tree_thq, tree_tth, "LD_lept");

  return 0;

}



void scanCut( DrawBase* db, TTree* tree_thq, TTree* tree_tth, const std::string& var ) {

  int nsteps = 20;
  float xMin = 0.; 
  if( var=="BDT_lept" ) xMin = -1.;
  float xMax = 1.;
  float stepsize = (xMax - xMin)/((float)nsteps);


  TGraph* gr_contam = new TGraph(0);

  for( unsigned istep=0; istep<nsteps; ++istep ) {

    float cut = xMin + (float)istep*stepsize;

    TH1D* mgg_tth = new TH1D("mgg_tth", "", 4, 115., 135.);
    TH1D* mgg_thq = new TH1D("mgg_thq", "", 4, 115., 135.);

    char selection[500];
    sprintf( selection, "eventWeight*(isLeptonic && %s > %f)", var.c_str(), cut );

    tree_tth->Project( "mgg_tth", "mgg", selection );
    tree_thq->Project( "mgg_thq", "mgg", selection );

    float tth = mgg_tth->Integral( 2, 3 ); //120-130
    float thq = mgg_thq->Integral( 2, 3 ); //120-130

    float tth_contam = 100.*tth/(thq*34.);

    gr_contam->SetPoint( istep, cut, tth_contam );

    delete mgg_tth;
    delete mgg_thq;

  }

  
  gr_contam->SetMarkerSize(1.3);
  gr_contam->SetMarkerColor(46);
  gr_contam->SetMarkerStyle(21);


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();


  TH2D* h2_axes = new TH2D("axes", "", 10, xMin, xMax, 10, 0., 40);
  if( var == "LD_lept" )
    h2_axes->SetXTitle( "LD Threshold" );
  else
    h2_axes->SetXTitle( "BDT Threshold" );
  h2_axes->SetYTitle( "ttH Contamination [%]");

  h2_axes->Draw();
  gr_contam->Draw("Psame");

  TLine* tenpercent = new TLine(xMin, 10., xMax, 10.);
  tenpercent->Draw("same");

  TPaveText* labeltop = db->get_labelTop();
  labeltop->Draw("same");

  gPad->RedrawAxis();
  
  std::string canvasName = db->get_outputdir() + "/tthContam_" + var + ".eps";
  c1->SaveAs(canvasName.c_str());

  delete h2_axes;
  delete c1;

}
