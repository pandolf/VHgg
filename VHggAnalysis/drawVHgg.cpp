#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "CommonTools/DrawBase.h"
#include "CommonTools/fitTools.h"





int main(int argc, char* argv[]) {

  if(  argc != 2 && argc != 3 ) {
    std::cout << "USAGE: ./drawVHgg [(string)selType] [bTaggerType=\"JP\"]" << std::endl;
    exit(23);
  }


  std::string selType(argv[1]);


  std::string bTaggerType = "JP";
  if( argc>=3 ) {
    std::string bTaggerType_str(argv[2]);
    bTaggerType = bTaggerType_str;
  }




  DrawBase* db = new DrawBase("VHgg");

  db->set_lumiOnRightSide();

  std::string outputdir_str = "VHggPlots_MConly_" + selType + "_" + bTaggerType;
  db->set_outputdir(outputdir_str);

  int signalFillColor = 46;

  std::string HToGGFileName = "VHgg_HToGG_M-125_8TeV-pythia6";
  HToGGFileName += "_" + selType;
  HToGGFileName += "_" + bTaggerType;
  HToGGFileName += ".root";
  TFile* HToGGFile = TFile::Open(HToGGFileName.c_str());
  db->add_mcFile( HToGGFile, "HToGG", "H (125)", signalFillColor, 0);

  std::string DiPhotonFileName = "VHgg_DiPhoton_8TeV-pythia6";
  DiPhotonFileName += "_" + selType;
  DiPhotonFileName += "_" + bTaggerType;
  DiPhotonFileName += ".root";
  TFile* DiPhotonFile = TFile::Open(DiPhotonFileName.c_str());
  db->add_mcFile( DiPhotonFile, "DiPhoton", "Diphoton", 29);

  std::string GammaJetFileName = "VHgg_GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  GammaJetFileName += "_" + selType;
  GammaJetFileName += "_" + bTaggerType;
  GammaJetFileName += ".root";
  TFile* GammaJetFile = TFile::Open(GammaJetFileName.c_str());
  db->add_mcFile( GammaJetFile, "GammaJet", "#gamma + Jet", kBlue);


  std::string DiBosonFileName = "VHgg_VV_8TeV";
  DiBosonFileName += "_" + selType;
  DiBosonFileName += "_" + bTaggerType;
  DiBosonFileName += ".root";
  TFile* DiBosonFile = TFile::Open(DiBosonFileName.c_str());
  db->add_mcFile( DiBosonFile, "DiBoson", "Di Boson", 40);

  std::string TriBosonFileName = "VHgg_VGG_8TeV";
  TriBosonFileName += "_" + selType;
  TriBosonFileName += "_" + bTaggerType;
  TriBosonFileName += ".root";
  TFile* TriBosonFile = TFile::Open(TriBosonFileName.c_str());
  db->add_mcFile( TriBosonFile, "TriBoson", "V + #gamma#gamma", 50);

  std::string TTFileName = "VHgg_TT_8TeV";
  TTFileName += "_" + selType;
  TTFileName += "_" + bTaggerType;
  TTFileName += ".root";
  TFile* TTFile = TFile::Open(TTFileName.c_str());
  db->add_mcFile( TTFile, "TT", "Top", kYellow+1);

  std::string QCDFileName = "VHgg_QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6";
  QCDFileName += "_" + selType;
  QCDFileName += "_" + bTaggerType;
  QCDFileName += ".root";
  TFile* QCDFile = TFile::Open(QCDFileName.c_str());
  db->add_mcFile( QCDFile, "QCD", "QCD", kGreen);



  db->set_shapeNormalization();




  bool log = true;


  db->drawHisto("njets", "Number of Jets", "", "Events", true );
  db->drawHisto("nbjets_loose", "Number of b-Jets (Loose)", "", "Events", true );
  db->drawHisto("nbjets_medium", "Number of b-Jets (Medium)", "", "Events", true );

  db->drawHisto("mjj");
  db->drawHisto("mjj_0btag");
  db->drawHisto("mjj_1btag");
  db->drawHisto("mjj_2btag");
  db->drawHisto("qgljet0");
  db->drawHisto("qgljet1");

  db->set_lumiNormalization(30000.);
  db->set_rebin(10);

  db->drawHisto("mgg_0btag");
  db->drawHisto("mgg_1btag");
  db->drawHisto("mgg_2btag");




  return 0;

}
