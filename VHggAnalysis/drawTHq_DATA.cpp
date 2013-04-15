#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawTHq [(string)redntpProdVersion][(string)selType] [bTaggerType=\"CSV\"]" << std::endl;
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
  DrawBase* db = new DrawBase("THq");



  db->set_lumiOnRightSide();
  //db->set_lumiNormalization();
  db->set_lumi(19500.);
  db->set_shapeNormalization();
  db->set_drawZeros(false);




  std::string outputdir_str = "THq_plots/"+redntpProdVersion+"/THqPlots_DATA_" + selType;
  outputdir_str += "_" + bTaggerType;
  db->set_outputdir(outputdir_str);



  int signalFillColor = 46;

  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";



  
  
  std::string dataFileName = "HGG_translated_tHq_loose.root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, "Data" );


  std::string THqFileName = inputDir + "/THq_tHqLeptonic_mH125_8TeV_testtest";
  THqFileName += "_" + selType;
  THqFileName += "_" + bTaggerType;
  THqFileName += ".root";
  TFile* THqFile = TFile::Open(THqFileName.c_str());
  db->add_mcFile( THqFile, "THq", "THq (125)", kBlack, 0);

  std::string ttHFileName = inputDir + "/THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  ttHFileName += "_" + selType;
  ttHFileName += "_" + bTaggerType;
  ttHFileName += ".root";
  TFile* ttHFile = TFile::Open(ttHFileName.c_str());
  db->add_mcFile( ttHFile,  "TTH", "ttH (125)", kRed+1, 0);


  //db->drawHisto_fromTree("tree_passedEvents", "pt_bJet", "eventWeight*(run>10000 || isLeptonic)", 20, 0., 250., "pt_bJet", "b-Jet Candidate p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight*(isLeptonic)", 20, -5., 5.0001, "eta_qJet", "q-Jet Candidate #eta");
  db->drawHisto_fromTree("tree_passedEvents", "charge_lept", "eventWeight*(isLeptonic)", 3, -1.5, 1.5, "charge_lept", "Lepton Charge");
  db->drawHisto_fromTree("tree_passedEvents", "mt_top", "eventWeight*(isLeptonic)", 20, 0., 500., "mt_top", "Top Transverse Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight*(isLeptonic)", 20, 0., 3.15, "deltaPhi_top_higgs", "#Delta#Phi(top-diphoton)", "rad");
  db->drawHisto_fromTree("tree_passedEvents", "deltaEta_lept_qJet", "eventWeight*(isLeptonic)", 20, 0., 5., "deltaEta_lept_qJet", "#Delta#eta(lepton-qJet)", "");
  db->drawHisto_fromTree("tree_passedEvents", "BDT_lept", "eventWeight*(isLeptonic)", 20, -1., 1.0001, "BDT_lept", "Leptonic Channel BDT");


  return 0;

}


