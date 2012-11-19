#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"

bool separate_signals = false;


void printYields( DrawBase* db, const std::string& suffix, bool doUL=false, DrawBase*db_separateSignals=0 );


int main(int argc, char* argv[]) {

  if(  argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawTTVHgg [(string)redntpVersion] [(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }



  std::string dataDataset = "DATA_Run2012ABC";

  std::string redntpVersion(argv[1]);
  std::string selType(argv[2]);


  std::string bTaggerType = "JP";
  if( argc>=4 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }



  std::string filedir = "./";


  DrawBase* db_stack = new DrawBase("TTVHgg_stack");
  DrawBase* db_separateSignals = new DrawBase("TTVHgg_separateSignals");

  db_stack->set_lumiOnRightSide();

  db_stack->set_lumiNormalization(9300.);

  std::string outputdir_str = "TTVHggPlots_DATA_" + dataDataset + "_" + redntpVersion + "_" + selType + "_" + bTaggerType;
  db_stack->set_outputdir(outputdir_str);



  int signalFillColor = 46;


  std::string dataFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_" + dataDataset;
  dataFileName += "_" + selType;
  dataFileName += "_" + bTaggerType;
  dataFileName += ".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db_stack->add_dataFile( dataFile, "data" );



  // inclusive signal file for stack plot
  std::string HToGGFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_HToGG_M-125_8TeV-pythia6";
  HToGGFileName += "_" + selType;
  HToGGFileName += "_" + bTaggerType;
  HToGGFileName += ".root";
  TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());
  db_stack->add_mcFile( HToGGFile, "HToGG", "H (125)", signalFillColor, 0);

  // separate signals needed later, put them in db_separateSignals:
  std::string GluGluHFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  GluGluHFileName += "_" + selType;
  GluGluHFileName += "_" + bTaggerType;
  GluGluHFileName += ".root";
  TFile* GluGluHFile = TFile::Open(GluGluHFileName.c_str());
  db_separateSignals->add_mcFile( GluGluHFile, "GluGluH", "ggF H (125)", signalFillColor, 0);

  std::string VBFHFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1";
  VBFHFileName += "_" + selType;
  VBFHFileName += "_" + bTaggerType;
  VBFHFileName += ".root";
  TFile* VBFHFile = TFile::Open(VBFHFileName.c_str());
  db_separateSignals->add_mcFile( VBFHFile, "VBFH", "VBF H (125)", signalFillColor, 0);

  std::string VHFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  VHFileName += "_" + selType;
  VHFileName += "_" + bTaggerType;
  VHFileName += ".root";
  TFile* VHFile = TFile::Open(VHFileName.c_str());
  db_separateSignals->add_mcFile( VHFile, "VH", "VH (125)", signalFillColor, 0);

  std::string TTHFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2";
  TTHFileName += "_" + selType;
  TTHFileName += "_" + bTaggerType;
  TTHFileName += ".root";
  TFile* TTHFile = TFile::Open(TTHFileName.c_str());
  db_separateSignals->add_mcFile( TTHFile, "TTH", "ttH (125)", signalFillColor, 0);



  // and now BG's:

  std::string DiPhotonFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db_stack->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);

  std::string GammaJetFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
  db_stack->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", 38);


  std::string DiBosonFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db_stack->add_mcFile( DiBosonFile, "DiBoson", "Diboson", 39);

  std::string TriBosonFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db_stack->add_mcFile( TriBosonFile, "TriBoson", "V + #gamma#gamma", 40);

  std::string TTFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db_stack->add_mcFile( TTFile, "TT", "Top", 44);

  std::string QCDFileName = filedir + "finalizedTrees_" + redntpVersion + "/TTVHgg_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  db_stack->add_mcFile( QCDFile, "QCD", "QCD", 41);




  //db_stack->set_displayEmptyDatasets(false);



  bool log = true;


  db_stack->set_yAxisMaxScale(1.1);
  db_stack->drawHisto("nvertex", "Number of Reconstructed Vertexes");
  db_stack->drawHisto("nvertex_PUW", "Number of Reconstructed Vertexes");
  db_stack->reset();

  db_stack->drawHisto("njets", "Number of Jets", "", "Events", log);
  db_stack->drawHisto("nbjets_loose", "Number of b-Jets (Loose)", "", "Events", log);
  db_stack->drawHisto("nbjets_medium", "Number of b-Jets (Medium)", "", "Events", log);


  db_stack->set_rebin(8);
  db_stack->drawHisto("mjj", "Dijet Mass", "GeV");
  db_stack->set_legendTitle( "VH, no tag");
  db_stack->drawHisto("mjj_VHnotag", "Dijet Mass", "GeV");
  //db_stack->drawHisto("mjj_cutMgg_VHnotag", "Dijet Mass (122 < m_{#gamma#gamma} < 128 GeV)", "GeV");

  db_stack->set_legendTitle( "VH, b-tagged");
  db_stack->drawHisto("mjj_VHbtag", "Dijet Mass", "GeV");
  //db_stack->drawHisto("mjj_cutMgg_VHbtag", "Dijet Mass (122 < m_{#gamma#gamma} < 128 GeV)", "GeV");


  //db_stack->set_legendTitle( "2 b-tag Med Category" );
  //db_stack->drawHisto("mjj_VH2btagmed", "Dijet Mass", "GeV");
  //db_stack->drawHisto("mjj_cutMgg_VH2btagmed", "Dijet Mass (122 < m_{#gamma#gamma} < 128 GeV)", "GeV");


  db_stack->set_rebin(2);
  db_stack->drawHisto("qgljet0", "Lead Jet Q-G LD");
  db_stack->drawHisto("qgljet1", "Sublead Jet Q-G LD");
  //db_stack->set_legendTitle( "0 b-tag Category" );
  //db_stack->drawHisto("qgljet0_VHnotag", "Lead Jet Q-G LD");
  //db_stack->drawHisto("qgljet1_VHnotag", "Sublead Jet Q-G LD");
  //db_stack->set_legendTitle( "1 b-tag Category" );
  //db_stack->drawHisto("qgljet_VHbtag", "Non b-Tagged Jet Q-G LD");
  //db_stack->set_legendTitle( "" );
  //db_stack->set_rebin();

  db_stack->set_rebin(5);
  db_stack->drawHisto("ptphot0", "Lead Photon p_{T}", "GeV");
  db_stack->drawHisto("ptphot1", "Sublead Photon p_{T}", "GeV");
  db_stack->drawHisto("ptjet0", "Lead Jet p_{T}", "GeV");
  db_stack->drawHisto("ptjet1", "Sublead Jet p_{T}", "GeV");

  db_stack->set_rebin(2);

  db_stack->drawHisto("ptDiphot", "Diphoton p_{T}", "GeV");

  db_stack->drawHisto("deltaPhi", "#Delta#Phi(dijet-diphoton)", "rad");
  db_stack->drawHisto("ptDijet", "Dijet p_{T}", "GeV");
  db_stack->drawHisto("ptRatio", "Dijet p_{T} / Diphoton p_{T}");
  db_stack->drawHisto("ptDifference", "Dijet p_{T} - Diphoton p_{T}", "GeV");

  db_stack->drawHisto("deltaEtaJets", "Jet-Jet #Delta#eta");
  db_stack->drawHisto("deltaFabsEtaJets", "Jet-Jet #Delta|#eta|");
  db_stack->drawHisto("zeppen", "Zeppenfeld Variable");

  db_stack->drawHisto("deltaPhi_kinfit", "#Delta#Phi(dijet-diphoton)", "rad");
  db_stack->drawHisto("ptDijet_kinfit", "Dijet p_{T}", "GeV");
  db_stack->drawHisto("ptRatio_kinfit", "Dijet p_{T} / Diphoton p_{T}");
  db_stack->drawHisto("ptDifference_kinfit", "Dijet p_{T} - Diphoton p_{T}", "GeV");

  db_stack->drawHisto("deltaEtaJets_kinfit", "Jet-Jet #Delta#eta");
  db_stack->drawHisto("deltaFabsEtaJets_kinfit", "Jet-Jet #Delta|#eta|");
  db_stack->drawHisto("zeppen_kinfit", "Zeppenfeld Variable");

  db_stack->drawHisto("cosTheta1");
  db_stack->drawHisto("cosTheta2");
  db_stack->drawHisto("cosThetaStar");
  db_stack->drawHisto("helphi");
  db_stack->drawHisto("helphi1");

  db_stack->drawHisto("cosThetaStar_jets");
  db_stack->drawHisto("helicityAngle_V", "sin(#Phi*)");

  db_stack->set_rebin(20);
  db_stack->drawHisto("mVstar", "V* Mass", "GeV");
  db_stack->set_rebin(5);
  db_stack->set_xAxisMax(200);
  db_stack->drawHisto("ptVstar", "V* p_{T}", "GeV");
  db_stack->set_xAxisMax();
  db_stack->drawHisto("etaVstar", "V* #eta");
  db_stack->drawHisto("phiVstar", "V* #phi", "rad");

  db_stack->drawHisto("mVstar_kinfit", "V* Mass", "GeV");
  db_stack->drawHisto("ptVstar_kinfit", "V* p_{T}", "GeV");
  db_stack->drawHisto("etaVstar_kinfit", "V* #eta");
  db_stack->drawHisto("phiVstar_kinfit", "V* #phi", "rad");
  
  db_stack->set_rebin(40);
  db_stack->drawHisto("kinfit_chiSquareProbMax", "KinFit Max #chi^{2} Prob");
  db_stack->set_legendTitle( "60 < m_{jj} < 120 GeV");
  db_stack->drawHisto("kinfit_chiSquareProbMax_mjjWindow", "KinFit Max #chi^{2} Prob");
  db_stack->set_xAxisMax(0.1);
  db_stack->set_rebin(20);
  db_stack->set_flags("zoom");
  db_stack->set_legendTitle( "");
  db_stack->drawHisto("kinfit_chiSquareProbMax", "KinFit Max #chi^{2} Prob");
  db_stack->set_legendTitle( "60 < m_{jj} < 120 GeV");
  db_stack->drawHisto("kinfit_chiSquareProbMax_mjjWindow", "KinFit Max #chi^{2} Prob");
  db_stack->reset();
  //db_stack->drawHisto_fromTree("tree_passedEvents", "chiSquareProbMax", "eventWeight", 20, 0., 1, "kinfit_chiSquareProbMax_prova", "KinFit Max #chi^{2} Prob");
  //db_stack->drawHisto_fromTree("tree_passedEvents", "chiSquareProbMax", "eventWeight*(mjj>60. && mjj<120.)", 20, 0., 0.1, "kinfit_chiSquareProbMax_mjjwindow", "KinFit Max #chi^{2} Prob");


  db_stack->set_rebin(5);

  db_stack->drawHisto("mgg_prepresel", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "prepresel" );
  db_stack->drawHisto("mgg_presel", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "presel");

  bool doUL = (selType != "presel" );

  db_stack->drawHisto("mgg", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "incl", doUL, db_separateSignals );

  db_stack->set_legendTitle( "ttH Leptonic Category");
  db_stack->drawHisto("mgg_ttH_leptonic", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "ttH_leptonic", doUL, db_separateSignals );
  db_stack->set_legendTitle( "VH Hadronic");
  db_stack->drawHisto("mgg_ttH_hadronic", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "ttH_hadronic", doUL, db_separateSignals );
  db_stack->set_legendTitle( "VH, no tag");
  db_stack->drawHisto("mgg_VHnotag", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "VHnotag", doUL, db_separateSignals );
  db_stack->set_legendTitle( "VH, b-tagged");
  db_stack->drawHisto("mgg_VHbtag", "DiPhoton Invariant Mass", "GeV");
  printYields( db_stack, "VHbtag", doUL, db_separateSignals );
  //db_stack->set_legendTitle( "2 b-tag Med Category" );
  //db_stack->drawHisto("mgg_VH2btagmed", "DiPhoton Invariant Mass", "GeV");
  //printYields( db_stack, "2tagmed", doUL );



  return 0;

}




void printYields( DrawBase* db, const std::string& suffix, bool doUL, DrawBase* db_separateSignals ) {

  float lumi_ul = 30000.;


  std::string yieldsFileName = db->get_outputdir() + "/yields_" + suffix + ".txt";
  ofstream yieldsFile(yieldsFileName.c_str());


  //float xMin = 122.;
  //float xMax = 128.;
  float xMin = 120.;
  float xMax = 130.;

  std::vector<TH1D*> histos = db->get_lastHistos_mc();

  std::string histoName = histos[0]->GetName();


  int binXmin = histos[0]->FindBin(xMin);
  int binXmax = histos[0]->FindBin(xMax) -1;
  //std::cout <<  binXmin << std::endl;
  //std::cout <<  binXmax << std::endl;

  bool foundSignal = false;
  float totalBG = 0.;
  float totalBG_ave = 0.;
  float signal = 0.;

  yieldsFile << std::endl << "Yields (@ 30 fb-1): " << std::endl;
  for( unsigned int ii=0; ii<histos.size(); ++ii ) {
    yieldsFile << db->get_mcFile(ii).datasetName << " " << histos[ii]->Integral(binXmin, binXmax)*lumi_ul/db->get_lumi() << std::endl;
    if( db->get_mcFile(ii).datasetName != "HToGG" ) {
      totalBG += histos[ii]->Integral(binXmin, binXmax)*lumi_ul/db->get_lumi();
      totalBG_ave += histos[ii]->Integral(1, histos[ii]->GetNbinsX())*lumi_ul/db->get_lumi();
    } else {
      foundSignal = true;
      signal = histos[ii]->Integral(binXmin, binXmax)*lumi_ul/db->get_lumi();
      // print different signal contributions:
      if( db_separateSignals!=0 ) {
        for( unsigned int jj=0; jj<db_separateSignals->get_mcFiles().size(); ++jj ) {
          TFile* thisFile = db_separateSignals->get_mcFile(jj).file;
          TH1D* thisHisto = (TH1D*)thisFile->Get(histoName.c_str());
          int this_binXmin = thisHisto->FindBin(xMin);
          int this_binXmax = thisHisto->FindBin(xMax);
          float thisYield = thisHisto->Integral(this_binXmin, this_binXmax)*lumi_ul;
          yieldsFile << "     " << db_separateSignals->get_mcFile(jj).datasetName << " " << thisYield << " (" << 100.*thisYield/signal << "%%)" <<std::endl;
        } //for separate signals
      } //if db_separatesignals
    } //if is signal
  }

  totalBG_ave *= (10.)/(histos[0]->GetXaxis()->GetXmax()-histos[0]->GetXaxis()->GetXmin());

  yieldsFile << "Total BG: " << totalBG << " (averaged: " << totalBG_ave << ")" << std::endl;

  float signal_xsec = 2.28E-03*(19.37 + 1.573 + 0.6966 + 0.3943 + 0.1302); 
  float total_signal = signal_xsec*lumi_ul;
  float effS = signal/total_signal;
  yieldsFile << "Signal efficiency: " << effS << std::endl;

  if( !foundSignal ) std::cout << "WARNING!!! DIDN'T FIND SIGNAL HToGG!" << std::endl; 

  
  if( doUL && foundSignal ) {

    float ul = CLA( lumi_ul, 0., effS, 0., totalBG, 0. );
    float ul_bgave = CLA( lumi_ul, 0., effS, 0., totalBG_ave, 0. );
    yieldsFile << std::endl << "No error on BG:" << std::endl;
    yieldsFile << "UL: " << ul << "    (average BG): " << ul_bgave << std::endl;
    yieldsFile << "UL/SM: " << ul/signal_xsec << "    (average BG): " << ul_bgave/signal_xsec << std::endl;
    float ul_bgerr = CLA( lumi_ul, 0., effS, 0., totalBG, 0.05*totalBG );
    float ul_bgerr_bgave = CLA( lumi_ul, 0., effS, 0., totalBG_ave, 0.05*totalBG_ave );
    yieldsFile << std::endl << "5\% error on BG:" << std::endl;
    yieldsFile << "UL: " << ul_bgerr << "    (average BG): " << ul_bgerr_bgave << std::endl;
    yieldsFile << "UL/SM: " << ul_bgerr/signal_xsec << "    (average BG): " << ul_bgerr_bgave/signal_xsec << std::endl;

  }

  std::cout << "-> Saved yields in " << yieldsFileName << std::endl;

}

