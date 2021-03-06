#include "RedNtpFinalizer_VHgg.h"
#include "TMath.h"
#include <iostream>



struct RedntpDirStruct {

  std::string maindir;
  std::string datadir;
  std::string mcdir;

};



void finalize_oneDataset(const std::string& redntpVersion,  const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets );
void do_haddCommand( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType );
std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType );
RedntpDirStruct get_dirs( const std::string& prodVersion );




int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 ) {
    std::cout << "USAGE: ./finalize_VHgg [redntpVersion] [dataset] [selectionType] [bTaggerType=\"JP\"]" <<std::endl;
    return 13;
  }


  std::string redntpVersion(argv[1]);
  std::string dataset(argv[2]);
  std::string selectionType(argv[3]);

  std::string bTaggerType="JP";
  if( argc>4 ) {
    std::string bTaggerType_str(argv[4]);
    bTaggerType = bTaggerType_str;
  }


  std::vector<std::string> *datasets = new std::vector<std::string>;


  if( dataset=="DATA_Run2011_FULL" ) {
   
    finalize_oneDataset(redntpVersion, "DoubleMu_Run2011A_FULL", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset(redntpVersion, "DoubleMu_Run2011B_v2", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset(redntpVersion, "DoubleElectron_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoubleElectron_Run2011B_v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "MuEG_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "MuEG_Run2011B", selectionType, bTaggerType, datasets);

  } else if( dataset=="DiPhotonBox_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="DiPhoton_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DiPhotonJets_8TeV-madgraph-tarball-v2_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="V_8TeV") {
  //} else if( dataset=="DYJetsToLL_8TeV-madgraph" ) {

    finalize_oneDataset(redntpVersion, "DYJetsToLL_M-10To50filter_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1" , selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "GJet_Pt40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6") {

    finalize_oneDataset(redntpVersion, "QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="TT_8TeV" ) {

    finalize_oneDataset(redntpVersion, "TTbarGG_0Jet_Summer12-PU_S7_START52_V9-v1" , selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTGJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  //} else if( dataset=="VG_8TeV" ) {

  } else if( dataset=="VV_8TeV" ) {

    //finalize_oneDataset("WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola");
    finalize_oneDataset(redntpVersion, "WW_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "ZZ_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WZ_TuneZ2star_8TeV_pythia6_tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "ZGToLLG_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="VGG_8TeV" ) {

    finalize_oneDataset(redntpVersion, "WmGG", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "WpGG", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "ZGG", selectionType, bTaggerType, datasets);

  } else if( dataset=="HToGG_M-125_8TeV-pythia6" ) {

    finalize_oneDataset(redntpVersion, "WH_ZH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "GluGluToHToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "VBF_HToGG_M-125_8TeV-powheg-pythia6_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2", selectionType, bTaggerType, datasets);

  } else {
  
    finalize_oneDataset(redntpVersion, dataset, selectionType, bTaggerType, datasets );

  }


  do_haddCommand(redntpVersion, dataset, datasets, selectionType, bTaggerType );

  return 0;

}




void finalize_oneDataset( const std::string& redntpProdVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets ) {

  TString dataset_tstr(dataset);

  std::cout << std::endl << std::endl << std::endl << "####     Finalizing " << dataset << std::endl;
  std::cout << "####     Selection: " << selectionType << std::endl;
  std::cout << "####     b-Tagger: " << bTaggerType << std::endl << std::endl;

  RedntpDirStruct dirs = get_dirs( redntpProdVersion );

  RedNtpFinalizer_VHgg* rf = new RedNtpFinalizer_VHgg( dataset, selectionType, bTaggerType );
  //rf->set_redNtpDir("/xrootdfs/cms/local/pandolf/HiggsGammaGamma/reduced/redntp.52xv5_VH_feasibility_signalOnly.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged");

  bool isData = ( dataset_tstr.Contains("Run2011") || dataset_tstr.Contains("Run2012") );

  std::string redNtpDir = dirs.maindir;
  if( isData ) redNtpDir = redNtpDir + "/" + dirs.datadir;
  else         redNtpDir = redNtpDir + "/" + dirs.mcdir;

  rf->set_redNtpDir(redNtpDir);
  rf->set_outputDir("finalizedTrees_"+redntpProdVersion);
  rf->addFile(dataset);
  rf->finalize();
  delete rf;

  datasets->push_back(dataset);

}



void do_haddCommand( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType ) {


  if( datasets->size()<=1 ) return;

  
  std::string hadd_command = "hadd -f " + getSingleFileName( redntpVersion, dataset, selectionType, bTaggerType );

  std::string suffix = "";
  for( unsigned i=0; i<datasets->size(); ++i )
    suffix = suffix + " " + getSingleFileName( redntpVersion, datasets->at(i), selectionType, bTaggerType );

  hadd_command += suffix;
  std::string rm_command = "rm " + suffix;
  system(hadd_command.c_str());
  if( dataset!="HToGG_M-125_8TeV-pythia6" ) //keep also separate signal processes:
    system(rm_command.c_str());

}



std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName = "finalizedTrees_" + redntpVersion + "/VHgg_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  return fileName;

}


RedntpDirStruct get_dirs( const std::string& prodVersion ) {

  RedntpDirStruct returnStruct;

  if( prodVersion=="Lisbon" ) {

    returnStruct.maindir = "/xrootdfs/cms/local/meridian/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv1_data.cicpfloose.scales-Lisbon-Hgg.v2/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.v1/merged";

  } else if( prodVersion=="micheli_20120930" ) {

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv1_data.cicpfloose.scales-Lisbon-Hgg.ttH_feasibility_v1/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.ttH_feasibility_v1/merged";

  } else {

    std::cout << "-> Unknown prodVersion: '" << prodVersion << "'! Exiting." << std::endl;
    exit(11);
   
  }

  return returnStruct;

}
