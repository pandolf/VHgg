#include "RedNtpFinalizer_TTVHgg.h"
#include "TMath.h"
#include <iostream>



struct RedntpDirStruct {

  std::string maindir;
  std::string datadir;
  std::string mcdir;

};


bool use_inverted_photID = false;


void finalize_oneDataset(const std::string& redntpVersion,  const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType, std::vector<std::string> *datasets );
void do_haddCommand( const std::string& redntpVersion, const std::string& dataset, std::vector<std::string> *datasets, const std::string& selectionType, const std::string& bTaggerType );
std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType );
RedntpDirStruct get_dirs( const std::string& prodVersion );




int main( int argc, char* argv[] ) {

  if( argc!=4 && argc!=5 && argc!=6 ) {
    std::cout << "USAGE: ./finalize_TTVHgg [redntpVersion] [dataset] [selectionType] [use_inverted_photID=\"false\"] [bTaggerType=\"JP\"]" <<std::endl;
    return 13;
  }


  std::string redntpVersion(argv[1]);
  std::string dataset(argv[2]);
  std::string selectionType(argv[3]);
  if( argc>4 ) {
    std::string invPhotID_str(argv[4]);
    use_inverted_photID = (invPhotID_str=="true");
  }

  if( use_inverted_photID )  
    std::cout << std::endl << "-> INVERTING PHOTON ID." << std::endl;


  std::string bTaggerType="JP";
  if( argc>5 ) {
    std::string bTaggerType_str(argv[5]);
    bTaggerType = bTaggerType_str;
  }


  std::vector<std::string> *datasets = new std::vector<std::string>;

  if(dataset=="DATA_Run2012_ALL"){

    finalize_oneDataset(redntpVersion, "Photon_Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "Photon-Run2012A-recover-06Aug2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-24Aug2012-v2", selectionType, bTaggerType, datasets);//-EcalRecover_11Dec2012-v1 is missing
    finalize_oneDataset(redntpVersion, "DoublePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);


  } else if (dataset=="DATA_Run2012_less_than_moriond"){

    finalize_oneDataset(redntpVersion, "Photon_Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-24Aug2012-v2", selectionType, bTaggerType, datasets);//-EcalRecover_11Dec2012-v1 is missing
    finalize_oneDataset(redntpVersion, "DoublePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);


  } else if (dataset=="DATA_Run2012_SinglePhoton_ALL"){
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-EcalRecover_11Dec2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-24Aug2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);

  }else if (dataset=="DATA_Run2012_SinglePhoton_ALL_2"){//run A alone since you do not have double photon
    finalize_oneDataset(redntpVersion, "Photon_Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "Photon-Run2012A-recover-06Aug2012-v1", selectionType, bTaggerType, datasets);
   
  } else if(dataset=="DATA_Run2012_DoublePhoton_ALL"){
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
       finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-24Aug2012-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);


    }else if (dataset=="DATA_Run2012_SinglePhoton_ALL_less_than_moriond"){
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-EcalRecover_11Dec2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-24Aug2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);

  }else if(dataset=="DATA_Run2012_DoublePhoton_ALL_less_than_moriond"){
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-PromptReco-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-24Aug2012-v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton_Run2012D-PromptReco-v1", selectionType, bTaggerType, datasets);


    }else  if( dataset=="DATA_Run2011_FULL" ) {
   
    finalize_oneDataset(redntpVersion, "DoubleMu_Run2011A_FULL", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset(redntpVersion, "DoubleMu_Run2011B_v2", selectionType, bTaggerType, datasets); //first muons! important!
    finalize_oneDataset(redntpVersion, "DoubleElectron_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoubleElectron_Run2011B_v2", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "MuEG_Run2011A_FULL", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "MuEG_Run2011B", selectionType, bTaggerType, datasets);

  } else if( dataset=="DATA_Run2012ABC" ) {
   
    finalize_oneDataset(redntpVersion, "Photon-Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-13Jul2012-v1_b", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012C-13Jul2012-v2", selectionType, bTaggerType, datasets);

  } else if( dataset=="DATA_Run2012AB" ) {
   
    finalize_oneDataset(redntpVersion, "Photon-Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);

  } else if( dataset=="DATA_Run2012SinglePhoton" ) {
   
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "SinglePhoton_Run2012C-PromptReco-v1", selectionType, bTaggerType, datasets);
   

  } else if ( dataset== "DATA_Run2012SinglePhoton_2"){

    finalize_oneDataset(redntpVersion, "Photon-Run2012A-13Jul2012-v1", selectionType, bTaggerType, datasets);

  }   else if(dataset == "DATA_Run2012DoublePhoton"){
    finalize_oneDataset(redntpVersion,"DoublePhoton-Run2012C-13Jul2012-v1_b", selectionType, bTaggerType, datasets); 
    finalize_oneDataset(redntpVersion,"DoublePhoton-Run2012B-13Jul2012-v1", selectionType, bTaggerType, datasets); 
    //problems with this dataset    finalize_oneDataset("DoublePhoton-Run2012C-13Jul2012-v2", selectionType, bTaggerType, datasets); 
    finalize_oneDataset(redntpVersion,"DoublePhoton-Run2012C-13Jul2012-v2", selectionType, bTaggerType, datasets); 

  }else if( dataset=="DiPhotonBox_8TeV-pythia6" ) {

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
    //    finalize_oneDataset(redntpVersion, "TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    //    finalize_oneDataset(redntpVersion, "TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2", selectionType, bTaggerType, datasets);

  //} else if( dataset=="VG_8TeV" ) {

  } else if(dataset == "TTgg_8TeV"){
    // finalize_oneDataset(redntpVersion, "TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2", selectionType, bTaggerType, datasets);
      finalize_oneDataset(redntpVersion, "TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTbarGG_0Jet_Summer12-PU_S7_START52_V9-v1" , selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "TTGJets_8TeV-madgraph_Summer12-PU_S7_START52_V9-v1", selectionType, bTaggerType, datasets);



  }else if( dataset=="VV_8TeV" ) {

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

  } else if(dataset=="onlyTTH"){

    finalize_oneDataset(redntpVersion, "TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2", selectionType, bTaggerType, datasets);

  }else if( dataset=="tHq_mH125_8TeV_testtest" ) {

    finalize_oneDataset(redntpVersion, "tHqHadronic_mH125_8TeV_testtest", selectionType, bTaggerType, datasets);
    finalize_oneDataset(redntpVersion, "tHqLeptonic_mH125_8TeV_testtest", selectionType, bTaggerType, datasets);

  }else {
  
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

  // lousy patch
  //  if( dataset=="tHqLeptonic_mH125_8TeV_testtest" || dataset=="tHqHadronic_mH125_8TeV_testtest" || dataset_tstr.BeginsWith("TT_CT10") || dataset_tstr.BeginsWith("T_") || dataset_tstr.BeginsWith("Tbar_") )
  //dirs.mcdir = "/xrootdfs/cms/local/pandolf/Higgs/reduced/redntp.53xv2.cicpfloose.scales-Lisbon-Hgg.THq_feasibility_v1/merged";

  RedNtpFinalizer_TTVHgg* rf = new RedNtpFinalizer_TTVHgg( dataset, selectionType, bTaggerType );
  //rf->set_redNtpDir("/xrootdfs/cms/local/pandolf/HiggsGammaGamma/reduced/redntp.52xv5_VH_feasibility_signalOnly.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged");

  bool isData = ( dataset_tstr.Contains("Run2011") || dataset_tstr.Contains("Run2012") );

  std::string redNtpDir = dirs.maindir;
  if( isData ) redNtpDir = redNtpDir + "/" + dirs.datadir;
  else         redNtpDir = redNtpDir + "/" + dirs.mcdir;

  rf->setDatasetName(dataset);
  rf->set_redNtpDir(redNtpDir);
  rf->set_outputDir("finalizedTrees_"+redntpProdVersion);
  rf->addFile(dataset);
  if( use_inverted_photID ) 
    rf->invertPhotID();
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
  if( dataset!="HToGG_M-125_8TeV-pythia6"&& dataset!="TT_8TeV" && dataset!="tHq_mH125_8TeV_testtest"  ) //keep also separate signal and top processes:
    system(rm_command.c_str());

}



std::string getSingleFileName( const std::string& redntpVersion, const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) {

  std::string fileName;

  if(!use_inverted_photID){
    fileName = "finalizedTrees_" + redntpVersion + "/TTVHgg_" + dataset + "_" + selectionType + "_" + bTaggerType + ".root";
  }else{
    fileName = "finalizedTrees_" + redntpVersion + "/TTVHgg_" + dataset + "_" + selectionType + "_invertedPhotID_" + bTaggerType + ".root";
  }
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

  } else if(prodVersion=="micheli_noPUID"){
    
    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv1_data.cicpfloose.scales-Lisbon-Hgg.ttH_feasibility_HLT/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.ttH_feasibility_HLT/merged";

  } else if(prodVersion=="micheli_20121112") {

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv1_data.cicpfloose.scales-Lisbon-Hgg.micheli_20121112/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.micheli_20121112/merged";

  } else if(prodVersion=="qg_version_2") {

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv1_data.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";

  }else if(prodVersion=="micheli_SinglePhoton_preselectionCS"){
    
    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv2_data_3.preselectionCS.scales-Lisbon-Hgg.singlePhotonHLT_v2/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.ttH_feasibility_HLT/merged";


  }else if(prodVersion== "micheli_DoublePhoton_preselectionCS"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir="redntp.53xv1_data.preselectionCS.scales-Lisbon-Hgg.ttH_feasibility_v1/merged";
    returnStruct.mcdir = "redntp.52xv5.preselectionCS.scales-Lisbon-Hgg.ttH_feasibility_v1/merged";


  }else if(prodVersion=="micheli_SinglePhoton_preselectionCS_20121112"){
    
    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir="redntp.53xv2_data_all.preselectionCS.scales-Lisbon-Hgg.micheli_20121112/merged";
    returnStruct.mcdir = "redntp.52xv5.preselectionCS.scales-Lisbon-Hgg.micheli_20121112/merged";

  }else if(prodVersion== "micheli_DoublePhoton_preselectionCS_20121112"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir="redntp.53xv1_data.preselectionCS.scales-Lisbon-Hgg.micheli_20121112/merged";
    returnStruct.mcdir = "redntp.52xv5.preselectionCS.scales-Lisbon-Hgg.micheli_20121112/merged";

  }else if( prodVersion=="VH_feasibility" ) {

    returnStruct.maindir = "/xrootdfs/cms/local/pandolf/HiggsGammaGamma/reduced/";
    returnStruct.datadir = "redntp.52xv5_VH_feasibility.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged";
    returnStruct.mcdir = "redntp.52xv5_VH_feasibility.cicpfloose.regrPho_eCorr_20062012.VH_feasibility_v0/merged";

  } else if ( prodVersion=="moriond_dataset"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv2_data.cicpfloose.scales-Moriond-Hgg_2.moriond_dataset/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";

  } else if ( prodVersion=="moriond_dataset_preselectionCS"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv2_data.preselectionCS.scales-Moriond-Hgg_2.moriond_dataset/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";

  } else if ( prodVersion=="less_than_moriond_dataset"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv3_data.cicpfloose.scales-Moriond-Hgg_2.less_than_moriond/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";

  } else if ( prodVersion=="less_than_moriond_dataset_preselectionCS"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv3_data.preselectionCS.scales-Moriond-Hgg_2.less_than_moriond/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.qg_version_2/merged";

  } else if ( prodVersion=="moriond_dataset_elevetoOnOnePho"){

    returnStruct.maindir = "/xrootdfs/cms/local/micheli/Higgs/reduced/";
    returnStruct.datadir = "redntp.53xv2_data.cicpfloosenoelevetoOnOnePho.scales-Lisbon-Hgg.moriond_dataset_eleveto/merged";
    returnStruct.mcdir = "redntp.52xv5.cicpfloosenoelevetoOnOnePho.scales-Lisbon-Hgg.qg_version_2/merged";

  }else if( prodVersion=="THq_feasibility" ) {

    returnStruct.maindir = "/";
    //returnStruct.maindir = "/xrootdfs/cms/local/pandolf/Higgs/reduced/";
    returnStruct.datadir = "/xrootdfs/cms/local/micheli/Higgs/reduced/redntp.53xv2_data.preselectionMVA.scales-Moriond-Hgg_2.moriond_dataset_chargeLeptons/merged";
    //returnStruct.mcdir = "redntp.53xv2.cicpfloose.scales-Lisbon-Hgg.THq_feasibility_v0/merged";
    returnStruct.mcdir   = "/xrootdfs/cms/local/pandolf/Higgs/reduced/redntp.52xv5.cicpfloose.scales-Lisbon-Hgg.THq_feasibility_v1/merged";

  }else{

    std::cout << "-> Unknown prodVersion: '" << prodVersion << "'! Exiting." << std::endl;
    exit(11);
   
  }

  return returnStruct;

}
