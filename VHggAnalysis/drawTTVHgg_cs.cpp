#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
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





  DrawBase* db_nostack = new DrawBase("TTVHgg_nostack");
  db_nostack->set_lumiOnRightSide();
  db_nostack->set_shapeNormalization();

  DrawBase* db_nostack_rw = new DrawBase("TTVHgg_nostack_rw");
  db_nostack_rw->set_lumiOnRightSide();
  db_nostack_rw->set_shapeNormalization();


  std::string outputdir_str = "TTVHgg_plots/CS/"+redntpProdVersion+"/TTVHggPlots_MConly_" + selType + "_" + bTaggerType;
  db_nostack->set_outputdir(outputdir_str);

  db_nostack_rw->set_outputdir(outputdir_str+"_rw/");


  int signalFillColor = 46;


  std::string inputDir="finalizedTrees_"+redntpProdVersion+"/";

  std::string dataDataset = "DATA_Run2012ABC";
  //DATA 
  std::string dataFileName = inputDir+"TTVHgg_" + dataDataset;//check if ok
  dataFileName += "_" + selType;
  dataFileName += "_" + bTaggerType;
  dataFileName += ".root";
  cout<<"using data file"<<dataFileName<<endl;
  TFile* dataFile = TFile::Open(dataFileName.c_str());
  db_nostack->add_dataFile( dataFile, "data");
  db_nostack_rw->add_dataFile( dataFile, "data"); 
 
  //Data control sample
  //  std::string controlSampleFileName ="finalizedTrees_micheli_SinglePhoton_preselectionCS/TTVHgg_DATA_Run2012SinglePhoton_cs_selection_"+bTaggerType+".root";
  //  std::string controlSampleFileName ="tree_forcs/control_sample_had.root";
  //  std::string controlSampleFileName ="tree_forcs/control_sample_lep.root";
  std::string controlSampleFileName ="cs_scaled_DATA2012.root";

  TFile* controlSampleFile = TFile::Open(controlSampleFileName.c_str());
  db_nostack->add_mcFile( controlSampleFile, "controlSample","CS",kRed ); 

  std::string controlSampleFileName_rw =inputDir+"TTVHgg_" + dataDataset;
  controlSampleFileName_rw += "_" + selType+"_inverted";
  controlSampleFileName_rw += "_" + bTaggerType;
  controlSampleFileName_rw += ".root";

  TFile* controlSampleFile_rw = TFile::Open(controlSampleFileName_rw.c_str());
  db_nostack_rw->add_mcFile( controlSampleFile_rw, "controlSample","CS",kRed ); 
  
 
  bool log = true;

  
  //  bool doUL = ("_cicpfloose_" + selType != "presel" );
  bool doUL= false;

  for(int category =0 ; category<5;++category){
    std::stringstream cat_ss;
    cat_ss<<category;
    string cat_str("(category=="+cat_ss.str()+")*(mgg>100 && mgg<180)");
    //    cout<<cat_str<<endl;
    db_nostack->set_flags("category_"+cat_ss.str());
    //  db_nostack->drawHisto_fromTree("tree_passedEvents","nvertex_PUW",cat_str,50,0,50,"nvertex_PUW","nvertex_PUW");
  db_nostack->drawHisto_fromTree("tree_passedEvents","njets",cat_str,10,0,10,"njets", "Number of Jets");
  
  db_nostack->set_flags("zoom");
  db_nostack->drawHisto_fromTree("tree_passedEvents","njets",cat_str,8,2,10,"njets", "Number of Jets","","","true");
  db_nostack->reset();
  db_nostack->drawHisto_fromTree("tree_passedEvents","nbjets_xloose",cat_str,10,0,10,"nbjets_loose", "Number of b-Jets (Loose)");
  db_nostack->drawHisto_fromTree("tree_passedEvents","nbjets_medium",cat_str,10,0,10,"nbjets_medium", "Number of b-Jets (Medium)");


  db_nostack->set_rebin(4);
 db_nostack->set_flags("category_"+cat_ss.str());
  db_nostack->drawHisto_fromTree("tree_passedEvents","mjj",cat_str,50,0,500, "mjj","Dijet Mass", "GeV");
  db_nostack->reset();
 db_nostack->set_flags("category_"+cat_ss.str());
  db_nostack->drawHisto_fromTree("tree_passedEvents","mgg",cat_str+"*(eventWeight*(mgg<120 || mgg>130))",16,100,180,"mgg", "DiPhoton Invariant Mass", "GeV");
  db_nostack->reset();

  db_nostack->set_flags("category_"+cat_ss.str());
  //  db_nostack->set_rebin(4);
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptPhot1",cat_str,25,60,260,"ptphot0", "Lead Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptPhot2",cat_str,12,20,108,"ptphot1", "Sublead Photon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etaPhot1",cat_str,30,-3.,3.,"etaphot0", "Lead Photon #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etaPhot2",cat_str,30,-3.,3.,"etaphot1", "Sublead Photon #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptJet1",cat_str,50,0,200,"ptjet0", "Lead Jet p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptJet2",cat_str,50,0,200,"ptjet1", "Sublead Jet p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","etaJet1",cat_str,30,-3.,3.,"etajet0", "Lead Jet #eta", "");
  db_nostack->drawHisto_fromTree("tree_passedEvents","etaJet2",cat_str,30,-3.,3.,"etajet1", "Sublead Jet #eta", "");

  db_nostack->drawHisto_fromTree("tree_passedEvents","epfMet",cat_str,150,0,150,"pfMet","pfMet","GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","Ht",cat_str,50,0,1000, "Ht","Ht", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptele1",cat_str,120,0,120,"ptele1", "Lead electron p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptele2",cat_str,120,0,120,"ptele2", "Sublead electron p_{T}", "GeV");

  db_nostack->drawHisto_fromTree("tree_passedEvents","ptmu1",cat_str,120,0,120,"ptmu1", "Lead muon p_{T}", "GeV");
  db_nostack->drawHisto_fromTree("tree_passedEvents","ptmu2",cat_str,120,0,120,"ptmu2", "Sublead muon p_{T}", "GeV");

  db_nostack->reset();
 
  }

  //plots not divided by category
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)",50,60,260,"ptphot0_inclusive","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)",25,20,120,"ptphot1_inclusive","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","(mgg<120 || mgg>130)",16,100,180,"mgg_inclusive", "DiPhoton Invariant Mass", "GeV");

  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)",30,-3,3,"etaphot0_inclusive","Lead Photon #eta", "");
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)",30,-3,3,"etaphot1_inclusive","Sublead Photon #eta", "");



  //eta weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight",50,60,260,"ptphot0_scaled_2D_etaBoth","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight",25,20,120,"ptphot1_scaled_2D_etaBoth","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","eta_scaled_2D_weight*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_etaBoth", "DiPhoton Invariant Mass", "GeV");

  //pt weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight",50,60,260,"ptphot0_scaled_2D_ptBoth","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight",25,20,120,"ptphot1_scaled_2D_ptBoth","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","pt_scaled_2D_weight*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_ptBoth", "DiPhoton Invariant Mass", "GeV");                       
  //double weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight*pt_scaled_2D_weight",50,60,260,"ptphot0_scaled_2D_ptEtaBoth","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight*pt_scaled_2D_weight",25,20,120,"ptphot1_scaled_2D_ptEtaBoth","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","eta_scaled_2D_weight*pt_scaled_2D_weight*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_ptEtaBoth", "DiPhoton Invariant Mass", "GeV");

  //data weights
  //eta weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data",50,60,260,"ptphot0_scaled_2D_etaBoth_data","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data",25,20,120,"ptphot1_scaled_2D_etaBoth_data","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","eta_scaled_2D_weight_data*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_etaBoth_data", "DiPhoton Invariant Mass", "GeV");

  //pt weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight_data",50,60,260,"ptphot0_scaled_2D_ptBoth_data","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight_data",25,20,120,"ptphot1_scaled_2D_ptBoth_data","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","pt_scaled_2D_weight_data*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_ptBoth_data", "DiPhoton Invariant Mass", "GeV");

  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight_data",30,-3.,3.,"etaphot0_scaled_2D_ptBoth_data","Lead Photon #eta", "");
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*pt_scaled_2D_weight_data",30,-3.,3.,"etaphot1_scaled_2D_ptBoth_data","Sublead Photon #eta", "");


  //eta weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data",30,-3.,3.,"etaphot0_scaled_2D_etaBoth_data","Lead Photon #eta", "");
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data",30,-3.,3.,"etaphot1_scaled_2D_etaBoth_data","Sublead Photon #eta", "");
  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","eta_scaled_2D_weight_data*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_etaBoth_data", "DiPhoton Invariant Mass", "GeV");
                       
  //double weight
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data*pt_scaled_2D_weight_data",50,60,260,"ptphot0_scaled_2D_ptEtaBoth_data","Lead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","ptPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data*pt_scaled_2D_weight_data",25,20,120,"ptphot1_scaled_2D_ptEtaBoth_data","Sublead Photon p_{T}", "GeV");
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot1","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data*pt_scaled_2D_weight_data",30,-3.,3.,"etaphot0_scaled_2D_ptEtaBoth_data","Lead Photon #eta", "");
  db_nostack_rw->drawHisto_fromTree("tree_weights","etaPhot2","(mgg>100 && mgg<180)*(mgg<120||mgg>130)*eta_scaled_2D_weight_data*pt_scaled_2D_weight_data",30,-3.,3.,"etaphot1_scaled_2D_ptEtaBoth_data","Sublead Photon #eta", "");



  db_nostack_rw->drawHisto_fromTree("tree_weights","mgg","eta_scaled_2D_weight_data*pt_scaled_2D_weight_data*(mgg<120 || mgg>130)",16,100,180,"mgg_scaled_2D_ptEtaBoth_data", "DiPhoton Invariant Mass", "GeV");



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

