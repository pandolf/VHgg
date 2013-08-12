#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"

bool separate_signals = true;
bool Ct_minus1 = true;

void printYields( DrawBase* db, const std::string& suffix, bool doUL=false, float massWindow=5. );
void drawBDTRoc( DrawBase* db );


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
  DrawBase* db_nostack_hadr = new DrawBase("THq_nostack_hadr");

  // stack is used for the MC stack after full selection
  DrawBase* db_stack = new DrawBase("THq_stack");




  db_nostack->set_lumiOnRightSide();
  db_nostack->set_shapeNormalization();

  db_nostack_hadr->set_lumiOnRightSide();
  db_nostack_hadr->set_shapeNormalization();

  db_stack->set_lumiOnRightSide();
  db_stack->set_lumiNormalization(20000.);
  db_stack->set_noStack(false);





  std::string outputdir_str = "THq_plots/"+redntpProdVersion+"/THqPlots_MConly_" + selType;
  outputdir_str += "_" + bTaggerType;
  if( Ct_minus1 ) outputdir_str += "_CtMinusOne";
  db_nostack->set_outputdir(outputdir_str);
  db_nostack_hadr->set_outputdir(outputdir_str);
  db_stack->set_outputdir(outputdir_str);



  int signalFillColor = 46;

  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";




  
  //std::string tHqFileName = "finalizedTrees_THq_feasibility/THq_tHqLeptonic_mH125_8TeV_testtest";
  std::string tHqFileName = inputDir + "/THq_tHq_mH125_8TeV_testtest";
  tHqFileName += "_" + selType;
  tHqFileName += "_" + bTaggerType;
  tHqFileName += ".root";
  TFile* tHqFile = TFile::Open(tHqFileName.c_str());
  if( !Ct_minus1 ) {
    db_nostack->add_mcFile( tHqFile, "tHq", "tHq", kBlack, 0);
    db_nostack_hadr->add_mcFile( tHqFile, "tHq", "tHq", kBlack, 0);
    db_stack->add_mcFile( tHqFile, "tHq", "tHq", kWhite, 0);
  } else {
    db_nostack->add_mcFile( tHqFile, "tHq", "tHq", kBlack, 0);
    db_nostack_hadr->add_mcFile( tHqFile, "tHq", "tHq", kBlack, 0);
    db_stack->add_mcFile( tHqFile, "tHq", "tHq (C_{t} = -1)", kWhite, 0);
    db_stack->set_mcWeight( "tHq", 34. );
  }


  std::string ttHFileName =inputDir +  "THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  ttHFileName += "_" + selType;
  ttHFileName += "_" + bTaggerType;
  ttHFileName += ".root";
  TFile* ttHFile = TFile::Open(ttHFileName.c_str());
  db_nostack->add_mcFile( ttHFile, "TTH", "ttH", kRed+1, 0);
  db_nostack_hadr->add_mcFile( ttHFile, "TTH", "ttH", kRed+1, 0);
  db_stack->add_mcFile( ttHFile, "TTH", "ttH", kGray, 0);
  if( Ct_minus1 )
    db_stack->set_mcWeight( "TTH", 2.4 );




  std::string VHFileName = inputDir +  "THq_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  VHFileName += "_" + selType;
  VHFileName += "_" + bTaggerType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  //db_nostack->add_mcFile( VHFile, "VH", "VH", signalFillColor+1, 0);
  db_stack->add_mcFile( VHFile, "VH", "VH", 42, 0);
  if( Ct_minus1 )
    db_stack->set_mcWeight( "VH", 2.4 );

  
  std::string GluGluHFileName = inputDir +  "THq_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  GluGluHFileName += "_" + selType;
  GluGluHFileName += "_" + bTaggerType;
  GluGluHFileName += ".root";
  TFile* GluGluHFile = TFile::Open(GluGluHFileName.c_str());
  //db_nostack->add_mcFile( GluGluHFile, "GluGluH", "ggF H ", signalFillColor+1, 0);
  db_stack->add_mcFile( GluGluHFile, "GluGluH", "ggF H", signalFillColor+1, 0);
  if( Ct_minus1 )
    db_stack->set_mcWeight( "GluGluH", 2.4 );


  std::string VBFHFileName = inputDir +  "THq_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  VBFHFileName += "_" + selType;
  VBFHFileName += "_" + bTaggerType;
  VBFHFileName += ".root";
  TFile* VBFHFile = TFile::Open(VBFHFileName.c_str());
  db_stack->add_mcFile( VBFHFile, "VBFH", "VBF H ", signalFillColor, 0);
  if( Ct_minus1 )
    db_stack->set_mcWeight( "VBFH", 2.4 );
  //    db_nostack->add_mcFile( VBFHFile, "VBFH", "VBF H", signalFillColor, 0);


  // inclusive signal file for stack plot
  std::string HToGGFileName = inputDir +  "THq_HToGG_M-125_8TeV-pythia6";
  HToGGFileName += "_" + selType;
  HToGGFileName += "_" + bTaggerType;
  HToGGFileName += ".root";
  TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());







  std::string DiPhotonFileName = inputDir +  "THq_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db_nostack_hadr->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 38);
  db_stack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);

  std::string GammaJetFileName = inputDir +  "THq_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
 //  db_nostack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 29);
  db_stack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);

  std::string DiBosonFileName = inputDir +  "THq_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db_stack->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);


  std::string TriBosonFileName = inputDir +  "THq_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db_stack->add_mcFile( TriBosonFile, "Vgg", "V#gamma#gamma", 40);


  std::string TTFileName = inputDir +  "THq_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db_stack->add_mcFile( TTFile, "TT", "Top", 44);
  //db_nostack->add_mcFile( TTFile, "TT", "Top", 39);

  std::string TTGGFileName = inputDir +  "THq_TTbarGG_0Jet_Summer12-PU_S7_START52_V9-v1";
  TTGGFileName += "_" + selType;
  TTGGFileName += "_" + bTaggerType;
  TTGGFileName += ".root";
  TFile* TTGGFile = TFile::Open(TTGGFileName.c_str());
  db_nostack->add_mcFile( TTGGFile, "TTGG", "tt+#gamma#gamma", 29);
  db_nostack_hadr->add_mcFile( TTGGFile, "TTGG", "tt+#gamma#gamma", 29);


  std::string QCDFileName = inputDir +  "THq_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  db_stack->add_mcFile( QCDFile, "QCD", "QCD", 41);








  bool log = true;


  db_nostack->set_widenLegend(false);

  db_nostack->drawHisto("nvertex", "Number of Primary Vertexes");
  db_nostack->drawHisto("nvertex_PUW", "Number of Primary Vertexes");
  db_nostack->reset();
  db_nostack->drawHisto("njets", "Number of Jets", "", "Events");
  db_nostack->drawHisto("nbjets_loose", "Number of b-Jets (Loose)", "", "Events");
  db_nostack->drawHisto("nbjets_medium", "Number of b-Jets (Medium)", "", "Events");
  db_nostack->reset();


  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptPhot1",   "eventWeight", 50, 0., 300.,  "ptPhot1", "Leading Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptPhot2",   "eventWeight", 50, 0., 150.,  "ptPhot2", "Subleading Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptRunPhot1",   "eventWeight", 50, 0., 300.,  "ptRunPhot1", "Running Leading Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptRunPhot2",   "eventWeight", 50, 0., 150.,  "ptRunPhot2", "Running Subleading Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptgg",   "eventWeight", 50, 0., 300.,  "ptgg", "Diphoton p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "ptRungg",   "eventWeight", 50, 0., 300.,  "ptRungg", "Diphoton Running p_{T}", "GeV");


  // leptonic channel plots:
  db_nostack->set_legendTitle("Leptonic Channel");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "isMu",   "eventWeight*isLeptonic", 2, 0., 2., "isMu", "Muonic Event");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "epfMet",   "eventWeight*isLeptonic", 50, 0., 250., "pfMet", "Particle Flow #slash{E}_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "BDT_lept",   "eventWeight*isLeptonic", 50, -1., 1.0001, "BDT_lept", "Leptonic BDT");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "BDT2_lept",   "eventWeight*isLeptonic", 50, -1., 1.0001, "BDT2_lept", "Leptonic BDT");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "BDTold_lept",   "eventWeight*isLeptonic", 50, -1., 1.0001, "BDTold_lept", "Leptonic BDT");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "LD_lept",   "eventWeight*isLeptonic", 50, 0., 1.0001, "LD_lept", "Leptonic LD");
//  db_nostack->drawHisto_fromTree("tree_passedEvents", "BDTBDT_lept",   "eventWeight*isLeptonic", 50, -1., 1.0001, "BDTBDT_lept", "BDT Leptonic BDT");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "njets",   "eventWeight*isLeptonic", 9, -0.5, 8.5, "njets_lept", "Number of Jets (p_{T} > 20 GeV)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "nbjets_loose",   "eventWeight*isLeptonic", 5, -0.5, 4.5,  "nbjets_loose_lept", "Number of b-Jets (Loose)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "nbjets_medium",   "eventWeight*isLeptonic", 5, -0.5, 4.5, "nbjets_medium_lept", "Number of b-Jets (Medium)");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_lept",  "eventWeight*isLeptonic", 50, 0., 200., "pt_lept", "Lepton p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_lept", "eventWeight*isLeptonic", 50, -2.5, 2.5, "eta_lept", "Lepton #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "charge_lept", "eventWeight*isLeptonic", 3, -1.5, 1.5, "charge_lept", "Lepton charge", "", "Events", false, 2);

  db_nostack->drawHisto_fromTree("tree_passedEvents", "nCentralJets", "eventWeight*isLeptonic", 7, -0.5, 6.5, "nCentralJets_lept", "Number of Additional Jets");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "hardestCentralJetPt", "eventWeight*isLeptonic", 50, 0., 200., "hardestCentralJetPt_lept", "Hardest Additional Jet p_{T}", "GeV", "Events", true);

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_qJet",  "eventWeight*isLeptonic", 50, 0., 250., "pt_qJet_lept", "q-Jet Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight*isLeptonic", 50, -5, 5, "eta_qJet_lept", "q-Jet Candidate #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "qgl_qJet", "eventWeight*isLeptonic", 50, 0., 1.0001, "qgl_qJet_lept", "q-Jet Candidate QGL", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight*isLeptonic", 50, 0., 3.15, "deltaPhi_top_higgs_lept", "#Delta#Phi (top-diphoton)", "rad");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaR_lept_phot", "eventWeight*isLeptonic", 50, 0., 6., "deltaR_lept_phot", "Minimal #DeltaR (lepton-photon)", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_top", "eventWeight*isLeptonic", 50, -5, 5, "eta_top_lept", "Top Candidate #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_lept-etagg", "eventWeight*isLeptonic", 50, -5., 5., "deltaEta_lept_higgs", "#Delta#eta (lepton-diphoton)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_lept-eta_qJet", "eventWeight*isLeptonic", 50, -5., 5., "deltaEta_lept_qJet_lept", "#Delta#eta (lepton-qJet)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet-eta_qJet", "eventWeight*isLeptonic", 50, -5., 5., "deltaEta_bJet_qJet_lept", "#Delta#eta (bJet-qJet)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet-eta_qJet", "eventWeight*isLeptonic*(abs(eta_qJet)>2.)", 50, -5., 5., "deltaEta_bJet_qJet_lept_cutEtaQJet", "#Delta#eta (bJet-qJet)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "mt_top", "eventWeight*(isLeptonic)", 50, 0., 500., "mt_top_lept",  "Top Transverse Mass", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "m_top", "eventWeight*(isLeptonic)", 50, 0., 500., "m_top_lept",  "Top Mass", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "mt_W", "eventWeight*(isLeptonic)", 50, 0., 250., "mt_W_lept",  "W Transverse Mass", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "cosThetaStar", "eventWeight*(isLeptonic)", 50, -1., 1.0001, "cosThetaStar_lept",  "cos(#theta*)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "cosThetaStar_tH", "eventWeight*(isLeptonic)", 50, -1., 1.0001, "cosThetaStar_tH_lept",  "cos(#theta*) (tH boost)");

  drawBDTRoc( db_nostack );


  // hadronic channel plots:
  db_nostack_hadr->set_legendTitle("Hadronic Channel");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "njets",   "eventWeight*(!isLeptonic)", 9, -0.5, 8.5, "njets_hadr", "Number of Jets (p_{T} > 20 GeV)");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "nbjets_loose",   "eventWeight*(!isLeptonic)", 5, -0.5, 4.5,  "nbjets_loose_hadr", "Number of b-Jets (Loose)");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "nbjets_medium",   "eventWeight*(!isLeptonic)", 5, -0.5, 4.5, "nbjets_medium_hadr", "Number of b-Jets (Medium)");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "nCentralJets", "eventWeight*(!isLeptonic)", 7, -0.5, 6.5, "nCentralJets_hadr", "Number of Additional Jets");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "hardestCentralJetPt", "eventWeight*(!isLeptonic)", 50, 0., 200., "hardestCentralJetPt_hadr", "Hardest Additional Jet p_{T}", "GeV", "Events", true);

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "pt_qJet",  "eventWeight*(!isLeptonic)", 50, 20., 120., "pt_qJet_hadr", "q-Jet Candidate p_{T}", "GeV");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight*(!isLeptonic)", 50, -5, 5, "eta_qJet_hadr", "q-Jet Candidate #eta", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "qgl_qJet", "eventWeight*(!isLeptonic)", 50, 0., 1.0001, "qgl_qJet_hadr", "q-Jet Candidate QGL", "");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "qgl_jetW1", "eventWeight*(!isLeptonic)", 50, 0., 1.0001, "qgl_jetW1_hadr", "Leading W-Jet QGL", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "qgl_jetW2", "eventWeight*(!isLeptonic)", 50, 0., 1.0001, "qgl_jetW2_hadr", "Subleading W-Jet QGL", "");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight*(!isLeptonic)", 50, 0., 3.15, "deltaPhi_top_higgs_hadr", "#Delta#Phi (top-diphoton)", "rad");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "eta_top-etagg", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaEta_top_higgs_hadr", "#Delta#eta (top-diphoton)", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "abs(eta_top)-abs(etagg)", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaFabsEta_top_higgs_hadr", "#Delta|#eta| (top-diphoton)", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "zeppen", "eventWeight*(!isLeptonic)", 50, -5., 5., "zeppen",  "Zeppendfeld Variable", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "eta_bJet-eta_qJet", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaEta_bJet_qJet_hadr", "#Delta#eta (bJet-qJet)");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "abs(eta_bJet)-abs(eta_qJet)", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaFabsEta_bJet_qJet_hadr", "#Delta|#eta| (bJet-qJet)");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "eta_top", "eventWeight*(!isLeptonic)", 50, -5, 5, "eta_top_hadr", "Top Candidate #eta", "");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "pt_top", "eventWeight*(!isLeptonic)", 50, 0., 250., "pt_top_hadr", "Top Candidate p_{T}", "GeV");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "m_top", "eventWeight*(!isLeptonic)", 50, 50., 500., "m_top_hadr",  "Top Mass", "GeV");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "m_W", "eventWeight*(!isLeptonic)", 50, 20., 250., "m_W_hadr",  "W Mass", "GeV");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "m_W", "eventWeight*(!isLeptonic && abs(m_top-172.5)<50.)", 50, 20., 250., "m_W_mTopCut_hadr",  "W Mass", "GeV");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "deltaR_top_phot_min_hadr", "eventWeight*(!isLeptonic)", 50, 0., 5., "deltaR_top_phot_min_hadr",  "Minimal #DeltaR(top, #gamma)", "");

  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "cosThetaStar", "eventWeight*(!isLeptonic)", 50, -1., 1.0001, "cosThetaStar_hadr",  "cos(#theta*)");
  db_nostack_hadr->drawHisto_fromTree("tree_passedEvents", "cosThetaStar_tH", "eventWeight*(!isLeptonic)", 50, -1., 1.0001, "cosThetaStar_tH_hadr",  "cos(#theta*) (tH boost)");



  // inclusive plots:

  db_nostack->set_legendTitle("");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_qJet",  "eventWeight", 50, 0., 500., "pt_qJet", "q-Jet Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_qJet", "eventWeight", 50, -5, 5, "eta_qJet", "q-Jet Candidate #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "qgl_qJet", "eventWeight", 50, 0., 1.0001, "qgl_qJet", "q-Jet Candidate QGL", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_bJet",  "eventWeight", 50, 0., 500., "pt_bJet", "b-Jet Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet", "eventWeight", 50, -5, 5, "eta_bJet", "b-Jet Candidate #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "pt_top",  "eventWeight", 50, 0., 500., "pt_top", "Top Candidate p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_top", "eventWeight", 50, -5, 5, "eta_top", "Top Candidate #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents", "deltaPhi_top_higgs", "eventWeight", 50, 0., 3.15, "deltaPhi_top_higgs", "#Delta#Phi (top-diphoton)", "rad");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_bJet-etagg", "eventWeight", 50, -5., 5., "deltaEta_bJet_higgs", "#Delta#eta (bJet-diphoton)");
  db_nostack->drawHisto_fromTree("tree_passedEvents", "eta_top -etagg", "eventWeight*(!isLeptonic)", 50, -5., 5., "deltaEta_top_higgs",  "#Delta#eta (top-diphoton)");




  //bool doUL = (selType != "presel" );
  bool doUL = true;
  //db_stack->set_rebin(5);
  db_stack->set_legendTitle("Leptonic Channel");
  db_stack->drawHisto("mgg_lept", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "lept", doUL );
  db_stack->drawHisto("mgg_lept_BDT", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "lept_BDT", doUL );
  db_stack->drawHisto("mgg_lept_LD", "DiPhoton Invariant Mass", "GeV");
  float massWindow = 0.5;
  for( ; massWindow<8.; massWindow+=0.5 ) 
    printYields( db_stack, "lept_LD", doUL, massWindow );


  db_stack->set_legendTitle("Hadronic Channel");
  db_stack->drawHisto("mgg_hadr", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "hadr", doUL );
  db_stack->drawHisto("mgg_hadr_centralJetVeto", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "hadr_centralJetVeto", doUL );
  db_stack->reset();



  ////bool doUL = (selType != "presel" );
  //bool doUL = true;

  //db_stack->drawHisto("mgg", "DiPhoton Invariant Mass", "GeV");
  //printYields( db_stack, "incl", doUL );





  return 0;

}




void printYields( DrawBase* db, const std::string& suffix, bool doUL, float massWindow ) {

  float xMin = 125.-massWindow;
  float xMax = 125.+massWindow;
  //float xMin = 120.;
  //float xMax = 130.;


  char yieldsFileName[1000];
  if( massWindow<1. )
    sprintf( yieldsFileName, "%s/yields_%s_0%.0f.txt", db->get_outputdir().c_str(), suffix.c_str(), massWindow*10. );
  else
    sprintf( yieldsFileName, "%s/yields_%s_%.0f.txt", db->get_outputdir().c_str(), suffix.c_str(), massWindow*10. );
  ofstream yieldsFile(yieldsFileName);


  std::vector<TH1D*> histos = db->get_lastHistos_mc();

  int binXmin = histos[0]->FindBin(xMin);
  int binXmax = histos[0]->FindBin(xMax) -1;

  bool foundSignal = false;
  float totalBG = 0.;
  float totalBG_ave = 0.;
  float totalBG_allSMH = 0.;
  float signal = 0.;
  float signal_noSMH = 0.;

  float massRange = (histos[0]->GetXaxis()->GetXmax()-histos[0]->GetXaxis()->GetXmin());

  yieldsFile << std::endl << "Yields (@ 20 fb-1): " << std::endl;

  for( unsigned int ii=0; ii<histos.size(); ++ii ) {

    std::string dataset = db->get_mcFile(ii).datasetName;

    bool isSMH = (dataset=="VBFH"||dataset=="TTH"||dataset=="VH"||dataset=="GluGluH");

    if( db->get_mcFile(ii).datasetName != "tHq" ) {

      if( isSMH ) { 
    
        // also SM H BR is enhanced by 2.4, but the enhancement over SM expectation is signal for us:
        float BG_SF = 1.;
        float S_SF = 0.;
        if( Ct_minus1 ) {
          BG_SF = 1./2.4;
          S_SF = 1.4/2.4;
        }

        float BG = histos[ii]->Integral(binXmin, binXmax)*BG_SF;

        float s = histos[ii]->Integral(binXmin, binXmax)*S_SF;

        totalBG += BG;
        totalBG_ave += BG;
        totalBG_allSMH += BG;

        signal += s;

        yieldsFile << db->get_mcFile(ii).datasetName << " (SM) " << BG << std::endl;
        yieldsFile << db->get_mcFile(ii).datasetName << " (extra) " << s << std::endl;

      } else {

        totalBG += histos[ii]->Integral(binXmin, binXmax);

      
        float SF = (xMax-xMin)/massRange; 
        float BG_ave = histos[ii]->Integral(1, histos[ii]->GetNbinsX())*SF;
        totalBG_ave += BG_ave;
        yieldsFile << db->get_mcFile(ii).datasetName << " " << histos[ii]->Integral(1, histos[ii]->GetNbinsX())*SF << std::endl;

      }

    } else {

      foundSignal = true;
      float thq = histos[ii]->Integral(binXmin, binXmax);
      signal += thq;
      signal_noSMH += thq;
      yieldsFile << db->get_mcFile(ii).datasetName << " " << thq << std::endl;
      //if( Ct_minus1 ) signal/=34.;

    }

  } //for datasets

  if( !foundSignal ) std::cout << "WARNING!!! DIDN'T FIND SIGNAL tHq!" << std::endl; 


  yieldsFile << "Total signal: " << signal << "  w/o SMH: " << signal_noSMH << std::endl;
  //yieldsFile << "Total BG / GeV: " << totalBG_ave/(xMax-xMin) << std::endl;
  yieldsFile << "Total BG (averaged): " << totalBG_ave << std::endl;
  yieldsFile << "Total SMH BG under peak: " << totalBG_allSMH << std::endl;
  yieldsFile << "Total extra SMH under peak: " << signal-signal_noSMH << std::endl;

  float signal_xsec = 2.28E-03*(0.0152*34. + 1.4*(19.37 + 1.573 + 0.6966 + 0.3943 + 0.1302));
  float signal_xsec_noSMH = 2.28E-03*(0.0152*34.);
  float total_signal = signal_xsec*db->get_lumi();
  float total_signal_noSMH = signal_xsec_noSMH*db->get_lumi();
  float effS = signal/total_signal;
  float effS_noSMH = signal_noSMH/total_signal_noSMH;
  yieldsFile << "Signal efficiency: " << effS << " w/o SMH: " << effS_noSMH << std::endl;


  
  if( doUL && foundSignal ) {

    float ul_bgave = CLA( db->get_lumi(), 0., effS, 0., totalBG_ave, 0. );
    float ul_bgave_noSMH = CLA( db->get_lumi(), 0., effS_noSMH, 0., totalBG_allSMH, 0. ); //fix this
    yieldsFile << std::endl << "No error on BG:" << std::endl;
    yieldsFile << "UL: " << ul_bgave << "  all SMH in BG: " << ul_bgave_noSMH << std::endl;
    yieldsFile << "UL/SM: " << ul_bgave/signal_xsec << "  all SMH in BG: " << ul_bgave_noSMH/signal_xsec_noSMH << std::endl;//fix this

    float ul_bgerrave = CLA( db->get_lumi(), 0., effS, 0., totalBG_ave, 0.30*totalBG_ave );
    float ul_bgerrave_noSMH = CLA( db->get_lumi(), 0., effS_noSMH, 0., totalBG_allSMH, 0.30*totalBG_allSMH );//fix this
    yieldsFile << std::endl << "30\% error on BG:" << std::endl;
    yieldsFile << "UL: " << ul_bgerrave << "  all SMH in BG: " << ul_bgerrave_noSMH << std::endl;
    yieldsFile << "UL/SM: " << ul_bgerrave/signal_xsec << "  all SMH in BG: " << ul_bgerrave_noSMH/signal_xsec_noSMH << std::endl;

  }


  std::cout << "-> Saved yields and ULs in " << yieldsFileName << std::endl; 

}




void drawBDTRoc( DrawBase* db ) {

  TH1F::AddDirectory(kTRUE);

  TFile* file_thq = db->get_mcFile("tHq").file;
  TFile* file_tth = db->get_mcFile("TTH").file;

  TTree* tree_thq = (TTree*)file_thq->Get("tree_passedEvents");
  TTree* tree_tth = (TTree*)file_tth->Get("tree_passedEvents");

  TH1F* h1_BDT_signal = new TH1F("BDT_signal", "", 100, -1., 1.0001); 
  TH1F* h1_BDT_bg = new TH1F("BDT_bg", "", 100, -1., 1.0001); 

  TH1F* h1_LD_signal = new TH1F("LD_signal", "", 100, 0., 1.0001); 
  TH1F* h1_LD_bg = new TH1F("LD_bg", "", 100, 0., 1.0001); 

  //TH1F* h1_LD_signal = new TH1F("LD_signal", "", 100, 0., 1.0001); 
  //TH1F* h1_LD_bg = new TH1F("LD_bg", "", 100, 0., 1.0001); 

  tree_thq->Project( "BDT_signal", "BDT_lept", "isLeptonic" );
  tree_thq->Project( "LD_signal", "LD_lept", "isLeptonic" );
  //tree_thq->Project( "LD_signal", "LD_lept", "isLeptonic" );

  tree_tth->Project( "BDT_bg", "BDT_lept", "isLeptonic" );
  tree_tth->Project( "LD_bg", "LD_lept", "isLeptonic" );
  //tree_tth->Project( "LD_bg", "LD_lept", "isLeptonic" );



  TGraph* gr_RoC_LD = new TGraph(0);
  TGraph* gr_RoC_BDT = new TGraph(0);

  int nbins = h1_LD_signal->GetNbinsX();

  for( unsigned int ibin=1; ibin<nbins+1; ++ibin ) {

  
    float eff_s_BDT = h1_BDT_signal->Integral( nbins-ibin, nbins )/h1_BDT_signal->Integral( 1, nbins );
    float eff_b_BDT = h1_BDT_bg->Integral( nbins-ibin, nbins )/h1_BDT_bg->Integral( 1, nbins );

    gr_RoC_BDT->SetPoint( ibin-1, 1.-eff_b_BDT, eff_s_BDT );
  
    float eff_s_LD = h1_LD_signal->Integral( nbins-ibin, nbins )/h1_LD_signal->Integral( 1, nbins );
    float eff_b_LD = h1_LD_bg->Integral( nbins-ibin, nbins )/h1_LD_bg->Integral( 1, nbins );
  
    gr_RoC_LD->SetPoint( ibin-1, 1.-eff_b_LD, eff_s_LD );

  }


  gr_RoC_BDT->SetMarkerSize(1.3);
  gr_RoC_BDT->SetMarkerStyle(24);
  gr_RoC_BDT->SetMarkerColor(kRed+3);

  gr_RoC_LD->SetMarkerSize(1.3);
  gr_RoC_LD->SetMarkerStyle(20);
  gr_RoC_LD->SetMarkerColor(kOrange+1);


  TCanvas* c1 = new TCanvas("c1_roc", "", 600, 600);
  c1->cd();

  TH2D* h2_axes = new TH2D("axes_roc", "", 10, 0., 1.0001, 10, 0., 1.0001);
  h2_axes->SetXTitle( "t#bar{t}H Rejection" );
  h2_axes->SetYTitle( "tHq Efficiency" );

  h2_axes->Draw();

  TLine* diag = new TLine(0., 1., 1., 0.);
  diag->Draw("same");


  TLegend* legend = new TLegend( 0.2, 0.2, 0.45, 0.45, "Leptonic Channel" );
  legend->SetFillColor(0);
  legend->SetTextSize(0.04);
  if( h1_BDT_signal!=0 && h1_BDT_bg!=0 )
    legend->AddEntry( gr_RoC_LD, "LD", "P");
  legend->AddEntry( gr_RoC_BDT, "BDT", "P");
  legend->Draw("same");

  TPaveText* labelTop = db->get_labelTop();
  labelTop->Draw("same");


  
  if( h1_BDT_signal!=0 && h1_BDT_bg!=0 ) 
    gr_RoC_LD->Draw("p same");
  gr_RoC_BDT->Draw("p same");

  gPad->RedrawAxis();

  char canvasName[500];
  sprintf( canvasName, "%s/RoC_BDT_lept.eps", db->get_outputdir().c_str() );
  c1->SaveAs(canvasName);

  delete c1;
  delete h2_axes;
  delete legend;
}

