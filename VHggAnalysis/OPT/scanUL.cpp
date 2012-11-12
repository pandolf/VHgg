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


float data_lumi=9300.;



int main( int argc, char* argv[] ) {


  if( argc!= 3 ) {
    std::cout << "USAGE: ./scanUL [optType] [category]" << std::endl;
    exit(11);
  }

  std::string selectionType = "prova";
  std::string selectionType_str(argv[1]);
  selectionType = selectionType_str;


  int category = -1;
  std::string category_str(argv[2]);
  category = atoi(category_str.c_str());


  //if( category<0 || category>5 ) {
  //  std::cout << "category must be between 0 and 5. Exiting." << std::endl;
  //  exit(11);
  //}


  // this sets the style:
  DrawBase* db = new DrawBase("OPT_UL");
  db->set_lumiOnRightSide();
  db->set_lumiNormalization(30000.);

  TPaveText* label_sqrt = db->get_labelSqrt();
    

  std::string optcutsdir = "optcuts_" + selectionType;
  char ULFileName_char[500];
  sprintf( ULFileName_char, "%s/ULscan_%d.txt", optcutsdir.c_str(), category);
  std::string ULFileName(ULFileName_char);

  ofstream ofs_UL(ULFileName.c_str());
  ofs_UL << "Expected for 30 fb-1:" << std::endl;
  ofs_UL << "Seff   \tS     \tB +- s(B)\t" << std::endl;

  TGraphErrors* gr_UL = new TGraphErrors(0);
  TGraphErrors* gr_UL_bgDATAave = new TGraphErrors(0);
  TGraphErrors* gr_UL_bgMCave = new TGraphErrors(0);
  float UL_max = 0.;
  float UL_min = 999.;
  float effS_UL_min = 0.;
  float effmax = 0.;

//  TFile* signalFile = TFile::Open("../VHgg_HToGG_M-125_8TeV-pythia6_presel_JP.root");
  //TFile* signalFile = TFile::Open("../VHgg_WH_ZH_HToGG_M-125_8TeV-pythia6_presel_JP.root");
  //TTree* signalTree = (TTree*)signalFile->Get("tree_passedEvents");

  // optimized working point chosen when looking only at VH (and ttH) signal:
  TChain* signalTree = new TChain("tree_passedEvents");
  if( category>1 ) // VH
    signalTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_JP.root");
  else // ttH
    signalTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_JP.root");


  TChain* backgroundTree = new TChain("tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_DiPhoton_8TeV-pythia6_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_VV_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_VGG_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_TT_8TeV_presel_JP.root/tree_passedEvents");
  backgroundTree->Add("../finalizedTrees_micheli_noPUID/TTVHgg_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6_presel_JP.root/tree_passedEvents");


  TChain* backgroundTree_data = new TChain("tree_passedEvents");
  backgroundTree_data->Add("../finalizedTrees_micheli_noPUID/TTVHgg_DATA_Run2012ABC_presel_JP.root");
  //backgroundTree_data->Add("/cmsrm/pc23/micheli/finalizedTrees_micheli_noPUID/TTVHgg_DATA_Run2012ABC_presel_invertedPhotID_JP.root");


  std::string categoryText;
  if( category==0 )
    categoryText = "ttH Leptonic";
  else if( category==1 )
    categoryText = "ttH Hadronic";
  else if( category==2 )
    categoryText = "VH 2 b-tag";
  else if( category==3 )
    categoryText = "VH 1 b-tag";
  else if( category==4 )
    categoryText = "VH 0 b-tag";



  for( unsigned iEff=1; iEff<=10; ++iEff ) {

    // use category 5 only for VH 2 and 1 tags
    int category_forFile = (category==2 || category==3 || category>4) ? 5 : category;

    char infileName[300];
    sprintf( infileName, "%s/cuts_cat%d_Seff%d.txt", optcutsdir.c_str(), category_forFile, iEff*10);
    ifstream ifs(infileName);
    std::cout << "-> Opening Seff file: " << infileName << std::endl;
  
    std::vector<std::string> varNames;
    std::vector<float> cutsMin;
    std::vector<float> cutsmin;

    while( ifs.good() && !ifs.eof() ) {

      std::string varName;
      float cutMin, cutmin;

      ifs >> varName >> cutMin >> cutmin;
    
      if( varName=="category" ) continue;

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

    // category:
    char categoryCut[300];
    sprintf( categoryCut, " && category==%d", category );
    std::string categoryCut_str(categoryCut);

    if( category<5 ) {
      selection += categoryCut_str;
    } else { //tried also these other btag combinations:
      if( category==6 ) selection += " && nbjets_medium>=1 && nbjets_loose>=2 ";
      if( category==7 ) selection += " && nbjets_medium>=2 ";
      if( category==8 ) selection += " && nbjets_medium>=1 ";
      if( category==9 ) selection += " && nbjets_loose>=1 ";
    }

    // just to be sure, add mgg cut by hand:
    selection += " && mgg>100. && mgg<180. ";

    std::string selection_bg(selection);
    std::string selection_sig(selection);

    // for signal: consider only signal region
    selection_sig += " && (mgg>120. && mgg<130.) ";

    // for BG: blind signal region (will then extrapolate from sidebands)
    selection_bg += " && (mgg<120. || mgg>130.) ";

//  char btagCut[200];
//  if( category==2 )
//    sprintf( btagCut, "nbjets_loose >= %d", category );
//  else
//    sprintf( btagCut, "nbjets_loose == %d", category );
//  std::string btagCut_str(btagCut);
//  selection += btagCut_str;

    // close parenthesis, you motherfucka
    selection += " )";
    selection_bg += " )";
    selection_sig += " )";

std::cout << "selection: " << selection_bg << std::endl;

    int nbinsx = 80.;
    TH1F* h1_bgDATA = new TH1F("bgDATA", "", nbinsx, 100., 180.);
    h1_bgDATA->Sumw2();
    TH1F* h1_bgMC = new TH1F("bgMC", "", nbinsx, 100., 180.);
    h1_bgMC->Sumw2();
    TH1F* h1_signal = new TH1F("signal", "", nbinsx, 100., 180.);
    h1_signal->Sumw2();
   
    signalTree         ->Project( "signal", "mgg", selection_sig.c_str() );
    backgroundTree     ->Project( "bgMC", "mgg", selection.c_str() );
    backgroundTree_data->Project( "bgDATA", "mgg", selection_bg.c_str() );

std::cout << "entries in BG(data): " << h1_bgDATA->GetEntries() << std::endl;

////TFile* file = TFile::Open("prova.root", "recreate");
////file->cd();
////h1_bg->Write();
////h1_signal->Write();
////h1_bgMC->Write();
////file->Close();
////exit(9);
    double lumi = 30000.;
    
    //int mggbinmin = h1_signal->FindBin(120.);
    //int mggbinmax = h1_signal->FindBin(130.);

    double signal = h1_signal->Integral();


    //double background_error;
    //double background = h1_bg->IntegralAndError( 1, nbinsx, background_error );

    // background from data: scaled from sidebands, assuming flat distribution:
    double backgroundDATA_ave_error;
    double backgroundDATA_ave = h1_bgDATA->IntegralAndError( 1, nbinsx, backgroundDATA_ave_error );
    backgroundDATA_ave *= (10./70.); 
    backgroundDATA_ave_error *= (10./70.);

    // background from MC: integrate over full range, assume flat, average:
    double backgroundMC_ave_error;
    double backgroundMC_ave = h1_bgMC->IntegralAndError( 1, nbinsx, backgroundMC_ave_error );
    backgroundMC_ave *= (10./80.); 
    backgroundMC_ave_error *= (10./80.);

    signal *= lumi;
    //background *= lumi;
    //background_error *= lumi;

    // BG from data corresponds to 12 fb-1
    backgroundDATA_ave *= lumi/data_lumi;
    backgroundDATA_ave_error *= lumi/data_lumi;

    backgroundMC_ave *= lumi;
    backgroundMC_ave_error *= lumi;



//std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$         bg: " << background << " bg_ave: " << background_ave << std::endl;
std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$         bg_ave: " << backgroundDATA_ave << " bgMC_ave: " << backgroundMC_ave << std::endl;
if( backgroundDATA_ave==0. && backgroundMC_ave==0. ) continue;
   
    //float signal_xsec = 2.28E-03*(19.37 + 1.573 + 0.6966 + 0.3943 + 0.1302); 
    // VH cross section only:
    float signal_xsec;
    if( category>1 ) // VH:
      signal_xsec = 2.28E-03*(0.6966 + 0.3943); 
    else // ttH
      signal_xsec = 2.28E-03*(0.1302);
    float total_signal = signal_xsec*db->get_lumi();
    float effS = signal/total_signal;
    //float UL = StatTools::computeUL( signal+background, background, background_error );
    //float UL = CLA( db->get_lumi(), 0., effS, 0., background, 0. );
    //float ULSM = UL/signal_xsec;

    float UL_bgDATAave = CLA( db->get_lumi(), 0., effS, 0., backgroundDATA_ave, 0. );
    float ULSM_bgDATAave = UL_bgDATAave/signal_xsec;

    float UL_bgMCave = CLA( db->get_lumi(), 0., effS, 0., backgroundMC_ave, 0. );
    float ULSM_bgMCave = UL_bgMCave/signal_xsec;


//std::cout << "signal: " << signal << " bg: " << background << " +- " << background_error << std::endl;
std::cout << "signal: " << signal << " bg: " << backgroundDATA_ave << " +- " << backgroundDATA_ave_error << std::endl;


    if( effS > effmax )
      effmax = effS;

    //gr_UL->SetPoint( iEff-1, 100.*effS, ULSM );
    gr_UL_bgDATAave->SetPoint( iEff-1, 100.*effS, ULSM_bgDATAave );
    gr_UL_bgMCave->SetPoint( iEff-1, 100.*effS, ULSM_bgMCave );

    if( ULSM_bgDATAave > UL_max ) {
      UL_max = ULSM_bgDATAave;
    }

    if( ULSM_bgDATAave < UL_min ) {
      UL_min = ULSM_bgDATAave;
      effS_UL_min = effS;
    }

    float ymax = h1_signal->GetMaximum() + h1_bgDATA->GetMaximum();
    ymax*=1.5;


    h1_bgDATA->Scale(30000./data_lumi);
    h1_bgMC->Scale(30000.);
    h1_signal->Scale(30000.);

    h1_bgMC->SetFillColor(38);
    h1_signal->SetFillColor(46);
    h1_bgDATA->SetMarkerStyle(20);
    h1_bgDATA->SetMarkerSize(1.5);

    h1_bgMC->Rebin(5.);
    h1_signal->Rebin(5.);
    h1_bgDATA->Rebin(5.);

    THStack* stack = new THStack();
    stack->Add( h1_bgMC );
    stack->Add( h1_signal );

    TH2D* h2_axes = new TH2D("axes", "", 10, 100., 180., 10, 0., h1_bgMC->GetMaximum()*1.5);
    h2_axes->SetXTitle("DiPhoton Invariant Mass [GeV]");
    h2_axes->SetYTitle("Events");


    TLegend* legend = new TLegend(0.6, 0.65, 0.9, 0.88, categoryText.c_str());
    legend->SetFillColor(0);
    legend->SetTextSize(0.035);
    legend->AddEntry( h1_signal, "Signal", "F");
    legend->AddEntry( h1_bgMC, "BG (MC)", "F");
    legend->AddEntry( h1_bgDATA, "BG (Data)", "P");

    char canvasName[250];
    sprintf( canvasName, "%s/yieldPlot_cat%d_Seff%d.eps", optcutsdir.c_str(), category, iEff*10);

    TPaveText* label = db->get_labelTop();

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

    TCanvas* c1 = new TCanvas("c1", "c1", 600., 600.);
    c1->cd();
    h2_axes->Draw();
    //h1_signal->Draw("same");
    stack->Draw("histo same");
    h1_bgDATA->Draw("E same");
    legend->Draw("same");
    label->Draw("same");
    //label_UL->Draw("same");
    //label_sqrt->Draw("same");
    gPad->RedrawAxis();
    c1->SaveAs(canvasName);

    delete c1;
    delete legend;
    delete h2_axes;
    //delete stack;
    

    ofs_UL << effS << "\t" << signal << "\t" << backgroundDATA_ave << " +- " << backgroundDATA_ave_error << "\t" << UL_bgDATAave << std::endl;

    delete h1_signal;
    delete h1_bgDATA;
    delete h1_bgMC;

std::cout << "### " << iEff << "   UL: " << UL_bgDATAave << "  UL/SM: " << UL_bgDATAave/signal_xsec << std::endl;
  } // for iEff

  std::cout << "> > >   BEST UL: " << UL_min << std::endl;
  std::cout << "> > >   signal eff: " << effS_UL_min << std::endl;

  ofs_UL.close();

  db->resetStyle();

//gr_UL->SetMarkerSize(2.);
//gr_UL->SetMarkerStyle(21);
//gr_UL->SetMarkerColor(kRed+3);

  gr_UL_bgDATAave->SetMarkerSize(2.);
  gr_UL_bgDATAave->SetMarkerStyle(20);
  gr_UL_bgDATAave->SetMarkerColor(kRed+2);

  gr_UL_bgMCave->SetMarkerSize(2.);
  gr_UL_bgMCave->SetMarkerStyle(29);
  gr_UL_bgMCave->SetMarkerColor(kOrange+1);

  TLegend* legendUL = new TLegend(0.5, 0.65, 0.9, 0.9, categoryText.c_str());
  legendUL->SetFillColor(0);
  legendUL->SetTextSize(0.04);
  legendUL->AddEntry(gr_UL_bgDATAave, "BG from Data", "P");
  legendUL->AddEntry(gr_UL_bgMCave, "BG from MC", "P");

  TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1.3*effmax*100., 10, 0., 1.6*UL_max ); 
  //TH2D* h2_axes_gr = new TH2D("axes_gr", "", 10, 0., 1., 10, 0., 5.);
  h2_axes_gr->SetYTitle("UL / SM (30 fb^{-1})");
  h2_axes_gr->SetXTitle("Signal Efficiency [%]");


  TCanvas* c_gr = new TCanvas("c_gr", "c_gr", 600., 600.);
  c_gr->cd();

  
  h2_axes_gr->Draw();
  //gr_UL->Draw("P same");
  gr_UL_bgDATAave->Draw("P same");
  gr_UL_bgMCave->Draw("P same");
  label_sqrt->Draw("same");
  legendUL->Draw("same");

  char UL_vs_Seff_name[250];
  sprintf(UL_vs_Seff_name, "%s/UL_vs_Seff_cat%d.eps", optcutsdir.c_str(), category );
  c_gr->SaveAs(UL_vs_Seff_name);

}




