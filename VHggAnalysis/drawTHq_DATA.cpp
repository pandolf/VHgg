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
  DrawBase* db_stack = new DrawBase("THq_stack");



  db->set_lumiOnRightSide();
  //db->set_lumiNormalization();
  db->set_lumi(19500.);
  db->set_shapeNormalization();
  db->set_drawZeros(false);

  db_stack->set_drawZeros(false);



  std::string outputdir_str = "THq_plots/"+redntpProdVersion+"/THqPlots_DATA_" + selType;
  outputdir_str += "_" + bTaggerType;
  db->set_outputdir(outputdir_str);
  db_stack->set_outputdir(outputdir_str);



  int signalFillColor = 46;

  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";



  
  
  std::string dataFileName = "HGG_translated_tHq_loose.root";
  //std::string dataFileName = inputDir + "/THq_DATA_Run2012";
  //dataFileName += "_" + selType;
  //dataFileName += "_" + bTaggerType;
  //dataFileName += ".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db->add_dataFile( dataFile, "data", "Sidebands Data" );
  db_stack->add_dataFile( dataFile, "data", "Data" );


  std::string THqFileName = inputDir + "/THq_tHq_mH125_8TeV_testtest";
  THqFileName += "_" + selType;
  THqFileName += "_" + bTaggerType;
  THqFileName += ".root";
  TFile* THqFile = TFile::Open(THqFileName.c_str());
  db->add_mcFile( THqFile, "THq", "THq (125)", kWhite, 0, -1, kBlack, 2);
  db_stack->add_mcFile( THqFile, "THq", "THq (Ct=-1)", kWhite, 0, -1, kBlack, 2);

  std::string ttHFileName = inputDir + "/THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  ttHFileName += "_" + selType;
  ttHFileName += "_" + bTaggerType;
  ttHFileName += ".root";
  TFile* ttHFile = TFile::Open(ttHFileName.c_str());
  db->add_mcFile( ttHFile,  "TTH", "ttH (125)", kRed+1, 0);
  db_stack->add_mcFile( ttHFile,  "TTH", "ttH (Ct=-1)", kGray, 0);


  //db->drawHisto_fromTree("tree_passedEvents", "pt_bJet", "eventWeight*(run>10000 || isLeptonic)", 20, 0., 250., "pt_bJet", "b-Jet Candidate p_{T}", "GeV");
  db->set_legendTitle("Leptonic Channel");
  db->drawHisto_fromTree("tree_passedEvents", "njets", "eventWeight*(isLeptonic)", 7, 1.5, 8.5, "nJets_lept", "Number of Jets");
  db->drawHisto_fromTree("tree_passedEvents", "nbjets_loose", "eventWeight*(isLeptonic)", 5, 0.5, 5.5, "nBJetsLoose_lept", "Number of Loose b-Jets");
  db->drawHisto_fromTree("tree_passedEvents", "nbjets_medium", "eventWeight*(isLeptonic)", 5, 0.5, 5.5, "nBJetsMedium_lept", "Number of Medium b-Jets");
  db->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight*(isLeptonic)", 20, -5., 5.0001, "eta_qJet_lept", "q-Jet Candidate #eta");
  db->drawHisto_fromTree("tree_passedEvents", "charge_lept", "eventWeight*(isLeptonic)", 3, -1.5, 1.5, "charge_lept_lept", "Lepton Charge");
  db->drawHisto_fromTree("tree_passedEvents", "mt_top", "eventWeight*(isLeptonic)", 20, 0., 500., "mt_top_lept", "Top Transverse Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight*(isLeptonic)", 20, 0., 3.15, "deltaPhi_top_higgs_lept", "#Delta#Phi(top-diphoton)", "rad");
  db->drawHisto_fromTree("tree_passedEvents", "deltaEta_lept_qJet", "eventWeight*(isLeptonic)", 20, 0., 5., "deltaEta_lept_qJet_lept", "#Delta#eta(lepton-qJet)", "");
  db->drawHisto_fromTree("tree_passedEvents", "BDT_lept", "eventWeight*(isLeptonic)", 40, -1., 1.0001, "BDT_lept", "Leptonic Channel BDT");
  db->drawHisto_fromTree("tree_passedEvents", "LD_lept", "eventWeight*(isLeptonic)", 40, 0., 1.0001, "LD_lept", "Leptonic Channel LD");


  db->set_legendTitle("Hadronic Channel");
  db->drawHisto_fromTree("tree_passedEvents", "njets", "eventWeight*(!isLeptonic)", 5, 3.5, 8.5, "nJets_hadr", "Number of Jets");
  db->drawHisto_fromTree("tree_passedEvents", "nCentralJets", "eventWeight*(!isLeptonic)", 7, -0.5, 6.5, "nCentralJets_hadr", "Number of Additional Jets (|#eta| < 2)");
  db->drawHisto_fromTree("tree_passedEvents", "nFwdJets", "eventWeight*(!isLeptonic)", 7, -0.5, 6.5, "nFwdJets_hadr", "Number of Additional Jets (|#eta| > 2)");
  db->drawHisto_fromTree("tree_passedEvents", "nbjets_loose", "eventWeight*(!isLeptonic)", 5, 0.5, 5.5, "nBJetsLoose_hadr", "Number of Loose b-Jets");
  db->drawHisto_fromTree("tree_passedEvents", "nbjets_medium", "eventWeight*(!isLeptonic)", 5, 0.5, 5.5, "nBJetsMedium_hadr", "Number of Medium b-Jets");
  db->drawHisto_fromTree("tree_passedEvents", "m_top", "eventWeight*(!isLeptonic)", 10, 130., 210, "m_top_hadr", "Top Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "m_W", "eventWeight*(!isLeptonic)", 10, 20., 220., "m_W_hadr", "W Mass", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "pt_qJet", "eventWeight*(!isLeptonic)", 20, 45., 245., "pt_qJet_hadr", "p_{T} (qJet)", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight*(!isLeptonic)", 20, -5., 5., "eta_qJet_hadr", "#eta (qJet)");
  db->drawHisto_fromTree("tree_passedEvents", "epfMet", "eventWeight*(!isLeptonic)", 20, 0., 200., "epfMet_hadr", "Particle Flow #slash{E}_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "Ht", "eventWeight*(!isLeptonic)", 20, 0., 500., "Ht_hadr", "H_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mHt", "eventWeight*(!isLeptonic)", 25, 0., 500., "mHt_hadr", "#slash{H}_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "mSt", "eventWeight*(!isLeptonic)", 25, 0., 500., "mSt_hadr", "#slash{S}_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight*(!isLeptonic)", 20, 0., 3.15, "deltaPhi_top_higgs_hadr", "#Delta#Phi(top-diphoton)", "rad");
  db->drawHisto_fromTree("tree_passedEvents", "deltaEta_bJet_qJet", "eventWeight*(!isLeptonic)", 10, 0., 5., "deltaEta_bJet_qJet_hadr", "#Delta#eta(bJet-qJet)", "");
  db->drawHisto_fromTree("tree_passedEvents", "ptPhot1", "eventWeight*(!isLeptonic)", 50, 35., 235., "ptPhot1_hadr", "Leading Photon p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptPhot1*120./mgg", "eventWeight*(!isLeptonic)", 50, 35., 235., "ptRunPhot1_hadr", "Running Leading Photon p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptPhot2", "eventWeight*(!isLeptonic)", 50, 20., 220., "ptPhot2_hadr", "Subleading Photon p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptPhot2*120./mgg", "eventWeight*(!isLeptonic)", 50, 20., 220., "ptRunPhot2_hadr", "Running Subleading Photon p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "ptgg", "eventWeight*(!isLeptonic)", 50, 0., 250., "ptgg_hadr", "DiPhoton p_{T}", "GeV");
  db->drawHisto_fromTree("tree_passedEvents", "deltaR_top_phot_min_hadr", "eventWeight*(!isLeptonic)", 20, 0., 5., "deltaR_top_phot_min_hadr",  "Minimal #DeltaR(top, #gamma)", "");
  db->drawHisto_fromTree("tree_passedEvents", "deltaR_qJet_phot_min_hadr", "eventWeight*(!isLeptonic)", 20, 0., 5., "deltaR_qJet_phot_min_hadr",  "Minimal #DeltaR(qJet, #gamma)", "");


  db_stack->set_lumiNormalization(19500.);
  db_stack->set_mcWeight( "THq", 34. );
  db_stack->set_mcWeight( "TTH", 2.4 );
  db_stack->set_yAxisMax( 6 );
  db_stack->set_legendTitle("Hadronic Channel");
  db_stack->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(!isLeptonic)", 80, 100., 180., "mgg_hadr", "Diphoton Mass", "GeV");
  db_stack->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(!isLeptonic && nCentralJets<=1 )", 80, 100., 180., "mgg_hadr_CJV", "Diphoton Mass", "GeV");
  db_stack->set_yAxisMax( 3.6 );
  db_stack->set_legendTitle("Leptonic Channel");
  db_stack->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(isLeptonic)", 80, 100., 180., "mgg_lept", "Diphoton Mass", "GeV");
  db_stack->drawHisto_fromTree("tree_passedEvents", "mgg", "eventWeight*(isLeptonic && LD_lept>0.25)", 80, 100., 180., "mgg_lept_LD", "Diphoton Mass", "GeV");


  // add all other processes:

  std::string VHFileName = inputDir +  "THq_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  VHFileName += "_" + selType;
  VHFileName += "_" + bTaggerType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  db->add_mcFile( VHFile, "VH", "VH", 42, 0);

  
  std::string GluGluHFileName = inputDir +  "THq_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  GluGluHFileName += "_" + selType;
  GluGluHFileName += "_" + bTaggerType;
  GluGluHFileName += ".root";
  TFile* GluGluHFile = TFile::Open(GluGluHFileName.c_str());
  db->add_mcFile( GluGluHFile, "GluGluH", "ggF H", signalFillColor+1, 0);


  std::string VBFHFileName = inputDir +  "THq_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  VBFHFileName += "_" + selType;
  VBFHFileName += "_" + bTaggerType;
  VBFHFileName += ".root";
  TFile* VBFHFile = TFile::Open(VBFHFileName.c_str());
  db->add_mcFile( VBFHFile, "VBFH", "VBF H ", signalFillColor, 0);


  db->reset();
//// inclusive signal file for stack plot
//std::string HToGGFileName = inputDir +  "THq_HToGG_M-125_8TeV-pythia6";
//HToGGFileName += "_" + selType;
//HToGGFileName += "_" + bTaggerType;
//HToGGFileName += ".root";
//TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());






  DrawBase* db2 = new DrawBase("db2");
  db2->set_lumiOnRightSide();
  db2->set_lumiNormalization(19500.);
  db2->set_outputdir(outputdir_str);

  db2->add_dataFile( dataFile, "data", "Sidebands Data" );

  std::string DiPhotonFileName = inputDir +  "THq_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db2->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);

  std::string GammaJetFileName = inputDir +  "THq_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
  db2->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);

  std::string DiBosonFileName = inputDir +  "THq_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db2->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);


  std::string TriBosonFileName = inputDir +  "THq_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db2->add_mcFile( TriBosonFile, "Vgg", "V#gamma#gamma", 40);


  std::string TTFileName = inputDir +  "THq_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db2->add_mcFile( TTFile, "TT", "Top", 44);



  std::string QCDFileName = inputDir +  "THq_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  //db->add_mcFile( QCDFile, "QCD", "QCD", 41);



  // and now set the weights 
  // - inflate the higgs processes
  // - set the sideband/(sideband+SR) scaling in order to compare to data (all)
  float SF = 60./80.;

  //db2->set_mcWeight( "THq", 34.*SF );
  //db2->set_mcWeight( "TTH", 2.4*SF );
  //db2->set_mcWeight( "VH", 2.4*SF );
  //db2->set_mcWeight( "GluGluH", 2.4*SF );
  //db2->set_mcWeight( "VBFH", 2.4*SF );
  db2->set_mcWeight( "DiPhoton", SF );
  db2->set_mcWeight( "GammaJet", SF );
  db2->set_mcWeight( "DiBoson", SF );
  db2->set_mcWeight( "Vgg", SF );
  db2->set_mcWeight( "TT", SF );

  db2->set_drawZeros(true);


  db2->set_getBinLabels(true);
  db2->drawHisto("cutFlow_lept", "", "", "Events", true);
  db2->set_getBinLabels(false);



  return 0;

}


