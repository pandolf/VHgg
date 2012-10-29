#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"
#include "cl95cms.C"

bool separate_mechanism = false;

void printYields( DrawBase* db, const std::string& suffix, bool doUL=false );

int main(int argc, char* argv[]) {

  if(  argc != 3 && argc != 4 ) {
    std::cout << "USAGE: ./drawTTVHgg [(string)redntpProdVersion][(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string  redntpProdVersion(argv[1]);
  std::string selType(argv[2]);


  std::string bTaggerType = "JP";
  if( argc>=4 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }



  std::string dataDataset = "DATA_Run2012_FULL";

  DrawBase* db_nostack = new DrawBase("ttHgg_nostack");
  db_nostack->set_lumiOnRightSide();
  db_nostack->set_shapeNormalization();


  std::string outputdir_str = "TTVHgg_plots/CS/"+redntpProdVersion+"/TTVHggPlots_MConly_" + selType + "_" + bTaggerType;
  db_nostack->set_outputdir(outputdir_str);

  int signalFillColor = 46;


  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";

  //DATA 
  std::string dataFileName = inputDir+"TTVHgg_" + dataDataset;//check if ok
  dataFileName += "_" + selType;
  dataFileName += "_" + bTaggerType;
  dataFileName += ".root";
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db_nostack->add_dataFile( dataFile, "data");
  
  //Data control sample
  std::string controlSampleFileName ="ttHgg_DATA_Run2012_FULL_cicpfloose_ttHsel3Jets1LeptonBTagVeto_inverted_JP.root";
  TFile* controlSampleFile = TFile::Open(controlSampleFileName.c_str());
  db_nostack->add_mcFile( controlSampleFile, "controlSample","CS",kRed );
  
  
  // inclusive higgs other mechanism
  std::string HToGGFileName = "TTVHgg_HToGG_M-125_8TeV-pythia6";
  HToGGFileName +=  "_cicpfloose_" + selType;
  HToGGFileName += "_" + bTaggerType;
  HToGGFileName += ".root";
  TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());
  db_nostack->add_mcFile( HToGGFile, "HToGG", "H (125) (other)", signalFillColor+2, 0);

  bool log = true;

  
  //  bool doUL = ("_cicpfloose_" + selType != "presel" );
  bool doUL= false;


  db_nostack->drawHisto_fromTree("tree_passedEvents","nvertex_PUW","",50,0,50,"nvertex_PUW","nvertex_PUW");
  db_nostack->drawHisto_fromTree("tree_passedEvents","njets","",10,0,10,"njets", "Number of Jets");
  
  db_nostack->set_flags("zoom");
  db_nostack->drawHisto_fromTree("tree_passedEvents","njets","",8,2,10,"njets", "Number of Jets","","","true");
  db_nostack->reset();
  db_nostack->drawHisto_fromTree("tree_passedEvents","nbjets_loose","",10,0,10,"nbjets_loose", "Number of b-Jets (Loose)");
  db_nostack->drawHisto_fromTree("tree_passedEvents","nbjets_medium","",10,0,10,"nbjets_medium", "Number of b-Jets (Medium)");


  db_nostack->set_rebin(4);
  db_nostack->drawHisto_fromTree("tree_passedEvents","mjj","",50,0,500, "mjj","Dijet Mass", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","mjj_0btag","",50,0,500, "mjj_0btag", "Dijet Mass", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","mjj_1btag","",50,0,500, "mjj_1btag", "Dijet Mass", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","mjj_2btag","",50,0,500, "mjj_2btag", "Dijet Mass", "GeV");

  db_nostack->reset();
  db_nostack->drawHisto_fromTree("tree_passedEvents","mgg","eventWeight*(mgg<120 || mgg>130)",16,100,180,"mgg", "DiPhoton Invariant Mass", "GeV");
  db_nostack->reset();

  db_nostack->set_rebin(4);
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptphot0","",50,0,200,"ptphot0", "Lead Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptphot1","",50,0,200,"ptphot1", "Sublead Photon p_{T}", "GeV");
  db_nostack->reset();

  db_nostack->drawHisto_fromTree("tree_passedEvents","etaphot0","",30,-3.,3.,"etaphot0", "Lead Photon #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etaphot1","",30,-3.,3.,"etaphot1", "Sublead Photon #eta", "");
  db_nostack->reset();


  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet0","",50,0,200,"ptjet0", "Lead Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet1","",50,0,200,"ptjet1", "Sublead Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet2","",50,0,200,"ptjet2", "3rd Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet3","",50,0,200,"ptjet3", "4th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet4","",50,0,200,"ptjet4", "5th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet5","",50,0,200,"ptjet5", "6th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet6","",50,0,200,"ptjet6", "7th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet7","",50,0,200,"ptjet7", "8th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet8","",50,0,200,"ptjet8", "9th Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptjet9","",50,0,200,"ptjet9", "10th Jet p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet0","",30,-3.,3.,"etajet0", "Lead Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet1","",30,-3.,3.,"etajet1", "Sublead Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet2","",30,-3.,3.,"etajet2", "3rd Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet3","",30,-3.,3.,"etajet3", "4th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet4","",30,-3.,3.,"etajet4", "5th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet5","",30,-3.,3.,"etajet5", "6th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet6","",30,-3.,3.,"etajet6", "7th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet7","",30,-3.,3.,"etajet7", "8th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet8","",30,-3.,3.,"etajet8", "9th Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etajet9","",30,-3.,3.,"etajet9", "10th Jet #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents","pfMet","",150,0,150,"pfMet","pfMet","GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","Ht","",100,0,1000, "Ht","Ht", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptele1","",120,0,120,"","ptele1", "Lead electron p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptele2","",120,0,120,"","ptele2", "Sublead electron p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptmu1","",120,0,120, "Lead muon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptmu2","",120,0,120, "Sublead muon p_{T}", "GeV");

  db_nostack->set_rebin(2);

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptDiphot","",50,0,500, "Diphoton p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","deltaPhi","",90,0,3, "#Delta#Phi(dijet-diphoton)", "rad");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptDijet","",50,0,500, "Dijet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptRatio","",30,0,3, "Dijet p_{T} / Diphoton p_{T}","");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptDifference","",100,-200.,200, "Dijet p_{T} - Diphoton p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","deltaEtaJets","",40,-10,10, "Jet-Jet #Delta#eta","");
  db_nostack->drawHisto_fromTree("tree_passedEvents","deltaFabsEtaJets","",40,-5,5, "Jet-Jet #Delta|#eta|","");

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
  std::cout <<  binXmin << std::endl;
  std::cout <<  binXmax << std::endl;

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

