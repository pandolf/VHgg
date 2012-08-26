#include "RedNtpFinalizer_VHgg.h"
#include "TMath.h"
#include <iostream>




void finalize_oneDataset( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets );
void do_haddCommand( const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType );
std::string getSingleFileName( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType );




int main( int argc, char* argv[] ) {

  if( argc!=3 && argc!=4 ) {
    std::cout << "USAGE: ./finalize_VHgg [dataset] [selectionType] [bTaggerType=\"JP\"]" <<std::endl;
    return 13;
  }


  std::string dataset(argv[1]);
  std::string selectionType(argv[2]);

  std::string bTaggerType="JP";
  if( argc>3 ) {
    std::string bTaggerType_str(argv[3]);
    bTaggerType = bTaggerType_str;
  }


  std::vector<std::string> *datasets = new std::vector<std::string>;


  if( dataset=="DATA_Run2011_FULL" ) {
   
    finalize_oneDataset("DoubleMu_Run2011A_FULL", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset("DoubleMu_Run2011B_v2", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset("DoubleElectron_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset("DoubleElectron_Run2011B_v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset("MuEG_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset("MuEG_Run2011B", selectionType, bTaggerType, datasets);

  } else if( dataset=="DiPhotonBox_8TeV-pythia6" ) {

    finalize_oneDataset("DiPhotonBox_Pt-10To25_8TeV-pythia6", selectionType, bTaggerType, datasets);
    finalize_oneDataset("DiPhotonBox_Pt-250ToInf_8TeV-pythia6", selectionType, bTaggerType, datasets);
    finalize_oneDataset("DiPhotonBox_Pt-25To250_8TeV-pythia6", selectionType, bTaggerType, datasets);

  } else if( dataset=="V_8TeV") {
  //} else if( dataset=="DYJetsToLL_8TeV-madgraph" ) {

    finalize_oneDataset("DYJetsToLL_M-10To50filter_8TeV-madgraph" , selectionType, bTaggerType, datasets);
    finalize_oneDataset("DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball", selectionType, bTaggerType, datasets);
    finalize_oneDataset("WJetsToLNu_TuneZ2Star_8TeV-madgraph", selectionType, bTaggerType, datasets);

  } else if( dataset=="GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6" ) {

    finalize_oneDataset("GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6", selectionType, bTaggerType, datasets);
    finalize_oneDataset("GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6", selectionType, bTaggerType, datasets);

  } else if( dataset=="QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6") {

    finalize_oneDataset("QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6", selectionType, bTaggerType, datasets);
    finalize_oneDataset("QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6", selectionType, bTaggerType, datasets);

  } else if( dataset=="TT_8TeV" ) {

    finalize_oneDataset("TTbarGG_0Jet_S1-cmkuo-TTGG_525_RECO_s46_v1" , selectionType, bTaggerType, datasets);
    finalize_oneDataset("TTGJets_8TeV-madgraph.root redntp_TTGJets_8TeV-madgraph_00", selectionType, bTaggerType, datasets);
    finalize_oneDataset("TTJets_TuneZ2star_8TeV-madgraph-tauola", selectionType, bTaggerType, datasets);

  } else if( dataset=="VG_8TeV" ) {

    //finalize_oneDataset("WGstarToLNu2E_TuneZ2star_8TeV-madgraph-tauola");
    //finalize_oneDataset("WGstarToLNu2Mu_TuneZ2star_7TeV-madgraph-tauola");
    //finalize_oneDataset("WGstarToLNu2Tau_TuneZ2star_7TeV-madgraph-tauola");
    finalize_oneDataset("WGToLNuG_TuneZ2star_8TeV-madgraph-tauola", selectionType, bTaggerType, datasets);
    finalize_oneDataset("ZGToLLG_8TeV-madgraph", selectionType, bTaggerType, datasets);

  } else if( dataset=="VV_8TeV" ) {

    //finalize_oneDataset("WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola");
    finalize_oneDataset("WW_TuneZ2star_8TeV_pythia6_tauola", selectionType, bTaggerType, datasets);
    finalize_oneDataset("ZZ_TuneZ2star_8TeV_pythia6_tauola", selectionType, bTaggerType, datasets);
    finalize_oneDataset("WZ_TuneZ2star_8TeV_pythia6_tauola", selectionType, bTaggerType, datasets);

  } else if( dataset=="VGG_8TeV" ) {

    finalize_oneDataset("WmGG_cmkuo", selectionType, bTaggerType, datasets);
    finalize_oneDataset("WpGG-cmkuo", selectionType, bTaggerType, datasets);
    finalize_oneDataset("ZGG-cmkuo", selectionType, bTaggerType, datasets);

  } else {
  
    finalize_oneDataset( dataset, selectionType, bTaggerType, datasets );

  }


  do_haddCommand( dataset, datasets, selectionType, bTaggerType );

  return 0;

}




void finalize_oneDataset( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets ) {

  std::cout << std::endl << "-> Finalizing " << dataset << std::endl;
  RedNtpFinalizer_VHgg* rf = new RedNtpFinalizer_VHgg( dataset, selectionType, bTaggerType );
  //rf->set_redNtpDir("/xrootdfs/cms/local/pandolf/HiggsGammaGamma/reduced/redntp.52xv5_VH_feasibility_signalOnly.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged");
  rf->set_redNtpDir("/xrootdfs/cms/local/pandolf/HiggsGammaGamma/reduced/redntp.52xv5_VH_feasibility.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged"); 
  rf->addFile(dataset);
  rf->finalize();
  delete rf;

  datasets->push_back(dataset);

}



void do_haddCommand( const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType ) {


  if( datasets->size()<=1 ) return;

  
  std::string hadd_command = "hadd -f " + getSingleFileName( dataset, selectionType, bTaggerType );

  std::string suffix = "";
  for( unsigned i=0; i<datasets->size(); ++i )
    suffix = suffix + " " + getSingleFileName( datasets->at(i), selectionType, bTaggerType );

  hadd_command += suffix;
  std::string rm_command = "rm " + suffix;
  system(hadd_command.c_str());
  system(rm_command.c_str());

}



std::string getSingleFileName( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName = "VHgg_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  return fileName;

}
