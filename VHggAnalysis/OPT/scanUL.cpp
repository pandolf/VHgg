#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TH2D.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include <iostream>
#include <fstream>

#include "DrawBase.h"

#include "../cl95cms.C"




int main( int argc, char* argv[] ) {

  std::string selectionType = "prova";
  if( argc>1 ) {
    std::string selectionType_str(argv[1]);
    selectionType = selectionType_str;
  }

  int nbtags = -1;
  if( argc>2 ) {
    std::string nbtags_str(argv[2]);
    nbtags = atoi(nbtags_str.c_str());
  }

  if( nbtags!=0 && nbtags!=1 && nbtags!=2 ) {
    std::cout << "nbtags must be 0, 1, or 2. Exiting." << std::endl;
    exit(11);
  }


  // this sets the style:
  DrawBase* db = new DrawBase("OPT_UL");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(30000.);

  TPaveText* label_sqrt = db->get_labelSqrt();
    

  std::string optcutsdir = "optcuts_" + selectionType;
  std::string ULFileName = optcutsdir + "/ULscan.txt";

  ofstream ofs_UL(ULFileName.c_str());
  ofs_UL << "Expected for 30 fb-1:" << std::endl;
  ofs_UL << "Seff   \tS     \tB +- s(B)\t" << std::endl;

  TGraphErrors* gr_UL = new TGraphErrors(0);
  float UL_max = 0.;
  float UL_min = 999.;
  float effS_UL_min = 0.;
  float effmax = 0.;

//  TFile* signalFile = TFile::Open("../VHgg_HToGG_M-125_8TeV-pythia6_presel_JP.root");
  // optimized working point chosen when looking only at VH signal:
  TFile* signalFile = TFile::Open("../VHgg_WH_ZH_HToGG_M-125_8TeV-pythia6_presel_JP.root");
  TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");


  TChain* backgroundTree = new TChain("tree_passedEvents");
  backgroundTree->Add("../VHgg_DiPhoton_8TeV-pythia6_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../VHgg_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../VHgg_VV_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../VHgg_VGG_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../VHgg_TT_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../VHgg_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6_presel_JP.root/tree_passedEvents");




  for( unsigned iEff=1; iEff<10; ++iEff ) {

    char infileName[300];
    sprintf( infileName, "%s/cuts_%dbtag_Seff%d.txt", optcutsdir.c_str(), nbtags, iEff*10);
    ifstream ifs(infileName);
    std::cout << "-> Opening Seff file: " << infileName << std::endl;
  
    std::vector<std::string> varNames;
    std::vector<float> cutsMin;
    std::vector<float> cutsmin;

    while( ifs.good() && !ifs.eof() ) {

      std::string varName;
      float cutMin, cutmin;

      ifs >> varName >> cutMin >> cutmin;

      varNames.push_back( varName );
      cutsMin.push_back( cutMin );
      cutsmin.push_back( cutmin );

    } //while file is good
  
    ifs.close();

    // eliminate last element (last line is read and is empty):
    varNames.pop_back();
    cutsMin.pop_back();
    cutsmin.pop_back();


    std::string selection = "eventWeight*( ";
    for( unsigned ivar=0; ivar<varNames.size(); ++ivar ) {
      if( ivar!=0 ) selection += " && ";
      char thisCut[200];
      sprintf( thisCut, "%s >= %f && %s < %f", varNames[ivar].c_str(), cutsMin[ivar], varNames[ivar].c_str(), cutsmin[ivar] );
      std::string thisCut_str(thisCut);
      selection += thisCut_str;
    }

//  char btagCut[200];
//  if( nbtags==2 )
//    sprintf( btagCut, "nbjets_loose >= %d", nbtags );
//  else
//    sprintf( btagCut, "nbjets_loose == %d", nbtags );
//  std::string btagCut_str(btagCut);
//  selection += btagCut_str;
    selection += " )";

std::cout << "selection: " << selection << std::endl;

    TH1F* h1_bg = new TH1F("bg", "", 80, 100., 180.);
    h1_bg->Sumw2();
    TH1F* h1_signal = new TH1F("signal", "", 80, 100., 180.);
    h1_signal->Sumw2();
   
    signalTree->Project( "signal", "mgg", selection.c_str() );
    backgroundTree->Project( "bg", "mgg", selection.c_str() );

    double lumi = 30000.;
    
    int mggbinmin = h1_signal->FindBin(120.);
    int mggbinmax = h1_signal->FindBin(130.);

    double signal = h1_signal->Integral(mggbinmin,mggbinmax);
    double background_error;
    double background = h1_bg->IntegralAndError( mggbinmin, mggbinmax, background_error );
   
    signal *= lumi;
    background *= lumi;
    background_error *= lumi;

    //float signal_xsec = 2.28E-03*(19.37 + 1.573 + 0.6966 + 0.3943 + 0.1302); 
    // VH cross section only:
    float signal_xsec = 2.28E-03*(0.6966); 
    float total_signal = signal_xsec*db->get_lumi();
    float effS = signal/total_signal;
    //float UL = StatTools::computeUL( signal+background, background, background_error );
    float UL = CLA( db->get_lumi(), 0., effS, 0., background, 0. );
    float ULSM = UL/signal_xsec;


std::cout << "signal: " << signal << " bg: " << background << " +- " << background_error << std::endl;


    if( effS > effmax )
      effmax = effS;

    gr_UL->SetPoint( iEff-1, 100.*effS, ULSM );

    if( ULSM > UL_max ) {
      UL_max = ULSM;
    }

    if( ULSM < UL_min ) {
      UL_min = ULSM;
      effS_UL_min = effS;
    }

    float ymax = h1_signal->GetMaximum() + h1_bg->GetMaximum();
    ymax*=1.5;

//  THStack* stack = new THStack();
//  stack->Add( h1_bg );
//  stack->Add( h1_signal );

//  TH2D* h2_axes = new TH2D("axes", "", 3, -0.5, 2.5, 10, 0., ymin);
//  h2_axes->GetXaxis()->SetLabelSize(0.085);
//  h2_axes->GetXaxis()->SetBinLabel(1, "#mu#mu");
//  h2_axes->GetXaxis()->SetBinLabel(2, "e#mu");
//  h2_axes->GetXaxis()->SetBinLabel(3, "ee");
//  h2_axes->SetYTitle("Events");


//  TLegend* legend = new TLegend(0.6, 0.75, 0.88, 0.88);
//  legend->SetFillColor(0);
//  legend->SetTextSize(0.035);
//  legend->AddEntry( h1_signal, "Signal", "F");
//  legend->AddEntry( h1_bg, "Background", "F");

//  char canvasName[250];
//  sprintf( canvasName, "%s/yieldPlot_Seff%d.eps", optcutsdir.c_str(), iEff*10);

    //TPaveText* label = new TPaveText( 0.15, 0.65, 0.45, 0.85, "brNDC");
    //label->SetFillColor(0);
    //label->SetTextSize(0.035);
    //label->AddText("L = 5 fb^{-1}");
    //char signalLabel[100];
    //sprintf( signalLabel, "s = %.2f (%d%%)", s, (int)(((float)h1_signal->GetEntries()/nTotal_s)*100) );
    //label->AddText( signalLabel );
    //char bgLabel[100];
    //sprintf( bgLabel, "b = %.2f", b_pred);
    //label->AddText( bgLabel );
    //char signifLabel[100];
    //sprintf( signifLabel, "UL = %.2f", UL);
    //label->AddText( signifLabel );

//  char b_text[500];
//  sprintf( b_text, "BG = %.2f #pm %.2f", b_pred, b_pred_err );
//  char obs_text[500];
//  sprintf( obs_text, "OBS = %.2f", obs );
//  char UL_text[100];
//  sprintf( UL_text, "UL = %.3f", UL );
//  TPaveText* label_UL = new TPaveText( 0.23, 0.73, 0.5, 0.88, "brNDC" );
//  label_UL->SetFillColor(0);
//  label_UL->SetTextSize(0.035);
//  label_UL->AddText(b_text);
//  label_UL->AddText(obs_text);
//  label_UL->AddText(UL_text);

//  TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
//  c1->cd();
//  h2_axes->Draw();
//  //h1_signal->Draw("same");
//  stack->Draw("histo same");
//  h1_bg->Draw("0 E2 same");
//  legend->Draw("same");
//  //label->Draw("same");
//  label_UL->Draw("same");
//  label_sqrt->Draw("same");
//  gPad->RedrawAxis();
//  c1->SaveAs(canvasName);

//  delete c1;
//  delete legend;
//  delete h2_axes;
    //delete stack;
    

    ofs_UL << effS << "\t" << signal << "\t" << background << " +- " << background_error << "\t" << UL << std::endl;

    delete h1_signal;
    delete h1_bg;

std::cout << "### " << iEff << "   UL: " << UL << "  UL/SM: " << UL/signal_xsec << std::endl;
  } // for iEff

  std::cout << "> > >   BEST UL: " << UL_min << std::endl;
  std::cout << "> > >   signal eff: " << effS_UL_min << std::endl;

  ofs_UL.close();

  db->resetStyle();

  gr_UL->SetMarkerSize(2.);
  gr_UL->SetMarkerStyle(21);
  gr_UL->SetMarkerColor(kRed+3);


  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effmax*100., 10, 0., 1.6*UL_max ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("UL / SM (30 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  gr_UL->Draw("P same");
  label_sqrt->Draw("same");

  char UL_vs_Seff_name[250];
  sprintf(UL_vs_Seff_name, "%s/UL_vs_Seff_%dbtag.eps", optcutsdir.c_str(), nbtags );
  c_gr->SaveAs(UL_vs_Seff_name);

}




