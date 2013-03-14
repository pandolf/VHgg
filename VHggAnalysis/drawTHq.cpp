#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"

bool separate_signals = true;
bool use_ttHSignal = true;


void printYields( DrawBase* db, const std::string& suffix, bool doUL=false );


int main(int argc, char* argv[]) {

  if(  argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawTHq [(string)redntpProdVersion][(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string  redntpProdVersion(argv[1]);
  std::string selType(argv[2]);


  std::string bTaggerType = "CSV";
  if( argc>=4 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }


  // no stack is used for shape comparisons between signal
  // and the main backgrounds (diphoton, gammajet and ttbar)
  DrawBase* db_nostack = new DrawBase("THq_nostack");

  // stack is used for the MC stack after full selection
  DrawBase* db_stack = new DrawBase("THq_stack");
  DrawBase* db_stack_UL = new DrawBase("THq_stack_UL");




  db_nostack->set_lumiOnRightSide();
  db_nostack->set_shapeNormalization();

  db_stack->set_lumiOnRightSide();
  db_stack->set_lumiNormalization(30000.);
  db_stack->set_noStack(false);

  db_stack_UL->set_lumiOnRightSide();
  db_stack_UL->set_lumiNormalization(30000.);
  db_stack_UL->set_noStack(false);            




  std::string outputdir_str = "THq_plots/"+redntpProdVersion+"/THqPlots_MConly_" + selType;
  outputdir_str += "_" + bTaggerType;
  db_nostack->set_outputdir(outputdir_str);
  db_stack->set_outputdir(outputdir_str);
  db_stack_UL->set_outputdir(outputdir_str+"/UL");



  int signalFillColor = 46;

  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";




  
  std::string tHqFileName = "finalizedTrees_THq_feasibility/THq_tHqLeptonic_mH125_8TeV_testtest_00";
  tHqFileName += "_" + selType;
  tHqFileName += "_" + bTaggerType;
  tHqFileName += ".root";
  TFile* tHqFile = TFile::Open(tHqFileName.c_str());
  db_nostack->add_mcFile( tHqFile, "tHq", "tHq (125)", kBlack, 0);
  db_stack->add_mcFile( tHqFile, "tHq", "tHq (125)", kWhite, 0);


  std::string ttHFileName =inputDir +  "THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  ttHFileName += "_" + selType;
  ttHFileName += "_" + bTaggerType;
  ttHFileName += ".root";
  TFile* ttHFile = TFile::Open(ttHFileName.c_str());
  db_nostack->add_mcFile( ttHFile, "TTH", "ttH (125)", kRed+1, 0);
  db_stack->add_mcFile( ttHFile, "TTH", "ttH (125)", kRed+3, 0);



  std::string VHFileName = inputDir +  "THq_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  VHFileName += "_" + selType;
  VHFileName += "_" + bTaggerType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  db_nostack->add_mcFile( VHFile, "VH", "VH (125)", kRed+3, 0);
  db_stack->add_mcFile( VHFile, "VH", "VH (125)", kRed+3, 0);

  
  std::string GluGluHFileName = inputDir +  "THq_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  GluGluHFileName += "_" + selType;
  GluGluHFileName += "_" + bTaggerType;
  GluGluHFileName += ".root";
  TFile* GluGluHFile = TFile::Open(GluGluHFileName.c_str());
  //if( !use_ttHSignal )
    db_nostack->add_mcFile( GluGluHFile, "GluGluH", "ggF H (125) ", signalFillColor+1, 0);
  db_stack->add_mcFile( GluGluHFile, "GluGluH", "ggF H (125)", signalFillColor+1, 0);

  std::string VBFHFileName = inputDir +  "THq_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  VBFHFileName += "_" + selType;
  VBFHFileName += "_" + bTaggerType;
  VBFHFileName += ".root";
  TFile* VBFHFile = TFile::Open(VBFHFileName.c_str());
  db_stack->add_mcFile( VBFHFile, "VBFH", "VBF H ", signalFillColor, 0);
  //    db_nostack->add_mcFile( VBFHFile, "VBFH", "VBF H", signalFillColor, 0);


  // inclusive signal file for stack plot
  std::string HToGGFileName = inputDir +  "THq_HToGG_M-125_8TeV-pythia6";
  HToGGFileName += "_" + selType;
  HToGGFileName += "_" + bTaggerType;
  HToGGFileName += ".root";
  TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());
  db_stack_UL->add_mcFile( HToGGFile, "HToGG", "H (125)", signalFillColor, 0);







  std::string DiPhotonFileName = inputDir +  "THq_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db_nostack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 38);
  db_stack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);
  db_stack_UL->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 38);

  std::string GammaJetFileName = inputDir +  "THq_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
 //  db_nostack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 29);
  db_stack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);
  db_stack_UL->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 29);

  std::string DiBosonFileName = inputDir +  "THq_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db_stack->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);
  db_stack_UL->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);


  std::string TriBosonFileName = inputDir +  "THq_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db_stack->add_mcFile( TriBosonFile, "Vgg", "V#gamma#gamma", 40);
  db_stack_UL->add_mcFile( TriBosonFile, "Vgg", "V#gamma#gamma", 40);


  std::string TTFileName = inputDir +  "THq_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db_stack->add_mcFile( TTFile, "TT", "Top", 44);
  db_stack_UL->add_mcFile( TTFile, "TT", "Top", 44);
  //if( use_ttHSignal )
  db_nostack->add_mcFile( TTFile, "TT", "Top", 29);


  std::string QCDFileName = inputDir +  "THq_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  db_stack->add_mcFile( QCDFile, "QCD", "QCD", 41);
  db_stack_UL->add_mcFile( QCDFile, "QCD", "QCD", 41);








  bool log = true;


  db_nostack->drawHisto("nvertex", "Number of Primary Vertexes");
  db_nostack->drawHisto("nvertex_PUW", "Number of Primary Vertexes");
  db_nostack->reset();
  if( !use_ttHSignal )
    db_nostack->set_xAxisMax(6.5);
  else 
    db_nostack->set_xAxisMax(8.5);
  db_nostack->drawHisto("njets", "Number of Jets", "", "Events");
  db_nostack->drawHisto("nbjets_loose", "Number of b-Jets (Loose)", "", "Events");
  db_nostack->drawHisto("nbjets_medium", "Number of b-Jets (Medium)", "", "Events");
  db_nostack->reset();


  db_nostack->drawHisto_fromTree("tree_passedEvents", "epfMet",   "eventWeight*isLeptonic", 50, 0., 250., "pfMet", "Particle Flow #slash{E}_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_lept",  "eventWeight", 50, 0., 200., "pt_lept", "Lepton p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_lept", "eventWeight", 50, -2.5, 2.5, "eta_lept", "Lepton #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_qJet",  "eventWeight", 50, 0., 500., "pt_qJet", "q-Jet Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight", 50, -5, 5, "eta_qJet", "q-Jet Candidate #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "qgl_qJet", "eventWeight", 50, 0., 1.0001, "qgl_qJet", "q-Jet Candidate QGL", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_bJet",  "eventWeight", 50, 0., 500., "pt_bJet", "b-Jet Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet", "eventWeight", 50, -5, 5, "eta_bJet", "b-Jet Candidate #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_top",  "eventWeight", 50, 0., 500., "pt_top", "Top Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_top", "eventWeight", 50, -5, 5, "eta_top", "Top Candidate #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "nCentralJets", "eventWeight", 7, -0.5, 6.5, "nCentralJets", "Number of Additional Jets");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "hardestCentralJetPt", "eventWeight", 50, 0., 200., "hardestCentralJetPt", "Hardest Additional Jet p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight", 50, 0., 3.15, "deltaPhi_top_higgs", "#Delta#Phi (top-diphoton)", "rad");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet-etagg", "eventWeight", 50, -5., 5., "deltaEta_bJet_higgs", "#Delta#eta (bJet-diphoton)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_lept-etagg", "eventWeight*isLeptonic", 50, -5., 5., "deltaEta_lept_higgs", "#Delta#eta (lepton-diphoton)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_top -etagg", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaEta_top_higgs",  "#Delta#eta (top-diphoton)");




  //bool doUL = (selType != "presel" );
  bool doUL = true;

  db_stack->drawHisto("mgg", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "incl", doUL );



  //Upper limits
  db_stack_UL->set_rebin(5);
  db_stack_UL->drawHisto("mgg", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack_UL, "all", doUL );


  return 0;

}




void printYields( DrawBase* db, const std::string& suffix, bool doUL ) {

  std::string yieldsFileName = db->get_outputdir() + "/yields_" + suffix + ".txt";
  ofstream yieldsFile(yieldsFileName.c_str());


  //float xMin = 122.;
  //float xMax = 128.;
  float xMin = 120.;
  float xMax = 130.;

  std::vector<TH1D*> histos = db->get_lastHistos_mc();

  int binXmin = histos[0]->FindBin(xMin);
  int binXmax = histos[0]->FindBin(xMax) -1;

  bool foundSignal = false;
  float totalBG = 0.;
  float totalBG_ave = 0.;
  float signal = 0.;

  yieldsFile << std::endl << "Yields (@ 30 fb-1): " << std::endl;
  for( unsigned int ii=0; ii<histos.size(); ++ii ) {
    yieldsFile << db->get_mcFile(ii).datasetName << " " << histos[ii]->Integral(binXmin, binXmax) << std::endl;
    if( db->get_mcFile(ii).datasetName != "HToGG" ) {
      totalBG += histos[ii]->Integral(binXmin, binXmax);
      totalBG_ave += histos[ii]->Integral(1, histos[ii]->GetNbinsX());
    } else {
      foundSignal = true;
      signal = histos[ii]->Integral(binXmin, binXmax);
    }
  }

  totalBG_ave *= (10.)/(histos[0]->GetXaxis()->GetXmax()-histos[0]->GetXaxis()->GetXmin());

  yieldsFile << "Total BG: " << totalBG << " (averaged: " << totalBG_ave << ")" << std::endl;

  float signal_xsec = 2.28E-03*(19.37 + 1.573 + 0.6966 + 0.3943 + 0.1302); 
  float total_signal = signal_xsec*db->get_lumi();
  float effS = signal/total_signal;
  yieldsFile << "Signal efficiency: " << effS << std::endl;

  if( !foundSignal ) std::cout << "WARNING!!! DIDN'T FIND SIGNAL HToGG!" << std::endl; 

  
  if( doUL && foundSignal ) {

    float ul = CLA( db->get_lumi(), 0., effS, 0., totalBG, 0. );
    float ul_bgave = CLA( db->get_lumi(), 0., effS, 0., totalBG_ave, 0. );
    yieldsFile << std::endl << "No error on BG:" << std::endl;
    yieldsFile << "UL: " << ul << "    (average BG): " << ul_bgave << std::endl;
    yieldsFile << "UL/SM: " << ul/signal_xsec << "    (average BG): " << ul_bgave/signal_xsec << std::endl;
    float ul_bgerr = CLA( db->get_lumi(), 0., effS, 0., totalBG, 0.05*totalBG );
    yieldsFile << std::endl << "5\% error on BG:" << std::endl;
    yieldsFile << "UL: " << ul_bgerr << std::endl;
    yieldsFile << "UL/SM: " << ul_bgerr/signal_xsec << std::endl;

  }



}

