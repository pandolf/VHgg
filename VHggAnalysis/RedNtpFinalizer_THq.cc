#include "RedNtpFinalizer_THq.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"
#include "THqStudies/THqLeptonicLikelihood/interface/THqLeptonicLikelihoodCalculator.h"

#include "BTagSFUtil/interface/BTagSFUtil.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"



#define DEBUG_EVENT_NUMBER_ 40124



RedNtpFinalizer_THq::RedNtpFinalizer_THq( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "THq", dataset ) 
{

  bTaggerType_ = bTaggerType;

  setSelectionType(selectionType);

  BLIND_ = true;

  use_thqPartonInfo_ = false;

}

RedNtpFinalizer_THq::~RedNtpFinalizer_THq()
{

   outFile_->Close();


   if (!tree_) return;
   delete tree_->GetCurrentFile();
}





void RedNtpFinalizer_THq::finalize()
{

   this->Init();

   TString dataset_tstr(dataset_);
   bool isSignalMC =  dataset_tstr.BeginsWith("tHq");


   std::string fullFlags = selectionType_ + "_";
   if(!invert_photonCuts_){
     fullFlags+=bTaggerType_;
   }else{
     fullFlags+="invertedPhotID_"+bTaggerType_;

   }
   this->set_flags(fullFlags); //this is for the outfile name
   this->createOutputFile();

   outFile_->cd();


   TH1D*  h1_nGenEvents = new TH1D("nGenEvents", "", 1, 0., 1.);
   h1_nGenEvents->Sumw2();

   TH1D*  h1_cutFlow_lept = new TH1D("cutFlow_lept", "", 6, 0., 6.);
   h1_cutFlow_lept->Sumw2();
   h1_cutFlow_lept->GetXaxis()->SetLabelSize(0.1);
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(1, "#gamma presel");
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(2, "#gamma running cut");
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(3, "b-tag");
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(4, "q-Jet");
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(5, "CJV");
   h1_cutFlow_lept->GetXaxis()->SetBinLabel(6, "LD");

   TH1D*  h1_cutFlow_hadr = new TH1D("cutFlow_hadr", "", 6, 0., 6.);
   h1_cutFlow_hadr->Sumw2();


   TH1D*  h1_nvertex = new TH1D("nvertex", "", 51, -0.5, 50.5);
   h1_nvertex->Sumw2();
   TH1D*  h1_nvertex_PUW = new TH1D("nvertex_PUW", "", 51, -0.5, 50.5);
   h1_nvertex_PUW->Sumw2();

   TH1D*  h1_njets_presel = new TH1D("njets_presel", "", 11, -0.5, 10.5);
   h1_njets_presel->Sumw2();
   TH1D*  h1_njets = new TH1D("njets", "", 11, -0.5, 10.5);
   h1_njets->Sumw2();
   TH1D*  h1_nbjets_loose = new TH1D("nbjets_loose", "", 11, -0.5, 10.5);
   h1_nbjets_loose->Sumw2();
   TH1D*  h1_nbjets_medium = new TH1D("nbjets_medium", "", 11, -0.5, 10.5);
   h1_nbjets_medium->Sumw2();


   TH1D* h1_ptphot0 = new TH1D("ptphot0", "",200, 0., 200.);
   h1_ptphot0->Sumw2();
   TH1D* h1_ptphot1 = new TH1D("ptphot1", "", 200, 0., 200.);
   h1_ptphot1->Sumw2();
   TH1D* h1_ptrunphot0 = new TH1D("ptrunphot0", "",200, 0., 200.);
   h1_ptrunphot0->Sumw2();
   TH1D* h1_ptrunphot1 = new TH1D("ptrunphot1", "", 200, 0., 200.);
   h1_ptrunphot1->Sumw2();
   TH1D* h1_etaphot0 = new TH1D("etaphot0", "", 100, -5., 5.);
   h1_etaphot0->Sumw2();
   TH1D* h1_etaphot1 = new TH1D("etaphot1", "", 100, -5., 5.);
   h1_etaphot1->Sumw2();


   //MET
   TH1D* h1_pfMet = new TH1D("pfMet", "", 150, 0., 150.);
   h1_pfMet->Sumw2();


   //HT
   TH1D* h1_Ht = new TH1D("Ht", "", 200, 0., 2000.);
   h1_Ht->Sumw2();

   //leptons
   TH1D* h1_pt_lept = new TH1D("pt_lept", "", 100, 0., 120.);
   h1_pt_lept->Sumw2();
   TH1D* h1_eta_lept = new TH1D("eta_lept", "", 100, -2.5, 2.5);
   h1_eta_lept->Sumw2();
   TH1D* h1_m_ele_phot = new TH1D("m_ele_phot", "", 100, 0., 200.);
   h1_m_ele_phot->Sumw2();
   TH1D* h1_deltaR_lept_phot1 = new TH1D("deltaR_lept_phot1", "", 100, 0., 10.);
   h1_deltaR_lept_phot1->Sumw2();
   TH1D* h1_deltaR_lept_phot2 = new TH1D("deltaR_lept_phot2", "", 100, 0., 10.);
   h1_deltaR_lept_phot2->Sumw2();


   // jets:
   TH1D* h1_pt_qJet = new TH1D("pt_qJet", "", 100, 0., 500.);
   h1_pt_qJet->Sumw2();
   TH1D* h1_eta_qJet = new TH1D("eta_qJet", "", 100, -5., 5.);
   h1_eta_qJet->Sumw2();
   TH1D* h1_qgl_qJet = new TH1D("qgl_qJet", "", 100, 0., 1.0001);
   h1_qgl_qJet->Sumw2();

   TH1D* h1_pt_bJet = new TH1D("pt_bJet", "", 100, 0., 500.);
   h1_pt_bJet->Sumw2();
   TH1D* h1_eta_bJet = new TH1D("eta_bJet", "", 100, -5., 5.);
   h1_eta_bJet->Sumw2();

   TH1D* h1_nCentralJets = new TH1D("nCentralJets", "", 7, -0.5, 6.5 );
   h1_nCentralJets->Sumw2();
   TH1D* h1_hardestCentralJetPt = new TH1D("hardestCentralJetPt", "", 100, 0., 200.);
   h1_hardestCentralJetPt->Sumw2();



   // top
   TH1D* h1_pt_top = new TH1D("pt_top", "", 100, 0., 120.);
   h1_pt_top->Sumw2();
   TH1D* h1_eta_top = new TH1D("eta_top", "", 100, -2.5, 2.5);
   h1_eta_top->Sumw2();
   TH1D* h1_m_top = new TH1D("m_top", "", 100, 0., 500.);
   h1_m_top->Sumw2();
   TH1D*  h1_deltaPhi_met_bJet = new TH1D("deltaPhi_met_bJet", "", 100, 0., 10.);
   h1_deltaPhi_met_bJet->Sumw2();



   TH1D* h1_mgg_prepresel = new TH1D("mgg_prepresel", "", 80, 100., 180.);
   h1_mgg_prepresel->Sumw2();
   TH1D* h1_mgg_presel = new TH1D("mgg_presel", "", 80, 100., 180.);
   h1_mgg_presel->Sumw2();

   TH1D* h1_mgg= new TH1D("mgg", "", 160, 100., 180.);
   h1_mgg->Sumw2();
   TH1D* h1_mgg_lept= new TH1D("mgg_lept", "", 160, 100., 180.);
   h1_mgg_lept->Sumw2();
   TH1D* h1_mgg_lept_BDT= new TH1D("mgg_lept_BDT", "", 160, 100., 180.);
   h1_mgg_lept_BDT->Sumw2();
   TH1D* h1_mgg_lept_LD= new TH1D("mgg_lept_LD", "", 160, 100., 180.);
   h1_mgg_lept_LD->Sumw2();
   TH1D* h1_mgg_hadr= new TH1D("mgg_hadr", "", 160, 100., 180.);
   h1_mgg_hadr->Sumw2();
   TH1D* h1_mgg_hadr_centralJetVeto= new TH1D("mgg_hadr_centralJetVeto", "", 160, 100., 180.);
   h1_mgg_hadr_centralJetVeto->Sumw2();




   // higgs:
   TH1D*  h1_ptDiphot = new TH1D("ptDiphot", "", 100, 0., 500.);
   h1_ptDiphot->Sumw2();
   TH1D*  h1_ptRunDiphot = new TH1D("ptRunDiphot", "", 100, 0., 500.);
   h1_ptRunDiphot->Sumw2();

   TH1D*  h1_deltaPhi_top_higgs = new TH1D("deltaPhi_top_higgs", "", 100, 0., 3.15);
   h1_deltaPhi_top_higgs->Sumw2();


   TH1D* h1_BDTlept_lept = new TH1D("BDTlept_lept", "", 100, -1., 1.0001);
   h1_BDTlept_lept->Sumw2();

   TH1D* h1_LDlept_lept = new TH1D("LDlept_lept", "", 100, 0., 1.0001);
   h1_LDlept_lept->Sumw2();


   // MC:

   TH1D* h1_deltaR_q = new TH1D("deltaR_q", "", 100, 0., 10.);
   h1_deltaR_q->Sumw2();
   TH1D* h1_deltaR_t = new TH1D("deltaR_t", "", 100, 0., 10.);
   h1_deltaR_t->Sumw2();
   TH1D* h1_deltaR_h = new TH1D("deltaR_h", "", 100, 0., 10.);
   h1_deltaR_h->Sumw2();
   TH1D* h1_deltaR_b = new TH1D("deltaR_b", "", 100, 0., 10.);
   h1_deltaR_b->Sumw2();
   TH1D* h1_deltaR_Wq = new TH1D("deltaR_Wq", "", 100, 0., 10.);
   h1_deltaR_Wq->Sumw2();
   TH1D* h1_deltaR_Wqbar = new TH1D("deltaR_Wqbar", "", 100, 0., 10.);
   h1_deltaR_Wqbar->Sumw2();
   TH1D* h1_deltaPhi_t = new TH1D("deltaPhi_t", "", 100, 0., 3.2);
   h1_deltaPhi_t->Sumw2();




   int njets_t;
   int nbjets_loose_t;
   int nbjets_medium_t;
   float ptPhot1_t;
   float ptPhot2_t;
   float ptRunPhot1_t;
   float ptRunPhot2_t;
   float etaPhot1_t;
   float etaPhot2_t;

   float mgg_t;
   float ptgg_t;
   float etagg_t;
   float ptRungg_t;

   float Ht_t, mHt_t, mSt_t;

   int hasPassedSinglePhot_t;
   int hasPassedDoublePhot_t;

   float mt_top_t;
   float m_top_t;
   float mt_W_t;
   float m_W_t;
   float BDT_lept_t;
   float BDT2_lept_t;
   float BDTold_lept_t;
   float LD_lept_t;

   bool isLeptonic_t;
   bool isMu_t;
   int   charge_lept_t;
   float pt_lept_t;
   float eta_lept_t;
   float pt_qJet_t;
   float eta_qJet_t;
   float qgl_qJet_t;
   float pt_jetW1_t;
   float eta_jetW1_t;
   float qgl_jetW1_t;
   float pt_jetW2_t;
   float eta_jetW2_t;
   float qgl_jetW2_t;
   float pt_bJet_t;
   float eta_bJet_t;
   float pt_top_t;
   float eta_top_t;
   int   nCentralJets_t;
   int   nFwdJets_t;
   float hardestCentralJetPt_t;
   float deltaPhi_top_higgs_t;
   float deltaEta_lept_qJet_t;
   float deltaEta_bJet_qJet_t;
   float deltaR_lept_phot_t;
   float m_ele_phot_t;
   float zeppen_t;
   float cosThetaStar_t;
   float cosThetaStar_tH_t;
   float deltaR_top_phot_min_hadr_t;
   float deltaR_qJet_phot_min_hadr_t;


   float eventWeight = 1.;

   TTree* tree_passedEvents = new TTree();
   //tree_passedEvents->SetName("thq_m125_8TeV");
   tree_passedEvents->SetName("tree_passedEvents");
   tree_passedEvents->Branch( "run", &run, "run/I" );
   tree_passedEvents->Branch( "lumi", &lumi, "lumi/I" );
   tree_passedEvents->Branch( "event", &event, "event/I" );
   tree_passedEvents->Branch( "eventWeight", &eventWeight, "eventWeight/F" );
   tree_passedEvents->Branch( "njets", &njets_t, "njets_t/I" );
   tree_passedEvents->Branch( "nbjets_loose", &nbjets_loose_t, "nbjets_loose_t/I" );
   tree_passedEvents->Branch( "nbjets_medium", &nbjets_medium_t, "nbjets_medium_t/I" );
   tree_passedEvents->Branch( "ptPhot1", &ptPhot1_t, "ptPhot1_t/F" );
   tree_passedEvents->Branch( "ptPhot2", &ptPhot2_t, "ptPhot2_t/F" );
   tree_passedEvents->Branch( "ptRunPhot1", &ptRunPhot1_t, "ptRunPhot1_t/F" );
   tree_passedEvents->Branch( "ptRunPhot2", &ptRunPhot2_t, "ptRunPhot2_t/F" );
   tree_passedEvents->Branch( "etaPhot1", &etaPhot1_t, "etaPhot1_t/F" );
   tree_passedEvents->Branch( "etaPhot2", &etaPhot2_t, "etaPhot2_t/F" );
   tree_passedEvents->Branch( "mgg", &mgg_t, "mgg_t/F" );
   tree_passedEvents->Branch( "ptgg", &ptgg_t, "ptgg_t/F" );
   tree_passedEvents->Branch( "etagg", &etagg_t, "etagg_t/F" );
   tree_passedEvents->Branch( "ptRungg", &ptRungg_t, "ptRungg_t/F" );
   tree_passedEvents->Branch("epfMet", &epfMet, "epfMet");
   tree_passedEvents->Branch( "hasPassedSinglePhot", &hasPassedSinglePhot_t, "hasPassedSinglePhot_t/I" );
   tree_passedEvents->Branch( "hasPassedDoublePhot", &hasPassedDoublePhot_t, "hasPassedDoublePhot_t/I" );
   tree_passedEvents->Branch("isLeptonic",          &isLeptonic_t,         "isLeptonic_t/O");
   tree_passedEvents->Branch("isMu",          &isMu_t,         "isMu_t/O");
   tree_passedEvents->Branch("charge_lept",         &charge_lept_t,        "charge_lept_t/I");
   tree_passedEvents->Branch("pt_lept",             &pt_lept_t,            "pt_lept_t/F");
   tree_passedEvents->Branch("eta_lept",            &eta_lept_t,           "eta_lept_t/F");
   tree_passedEvents->Branch("pt_qJet",             &pt_qJet_t,            "pt_qJet_t/F");
   tree_passedEvents->Branch("eta_qJet",            &eta_qJet_t,           "eta_qJet_t/F");
   tree_passedEvents->Branch("qgl_qJet",            &qgl_qJet_t,           "qgl_qJet_t/F");
   tree_passedEvents->Branch("pt_jetW1",             &pt_jetW1_t,            "pt_jetW1_t/F");
   tree_passedEvents->Branch("eta_jetW1",            &eta_jetW1_t,           "eta_jetW1_t/F");
   tree_passedEvents->Branch("qgl_jetW1",            &qgl_jetW1_t,           "qgl_jetW1_t/F");
   tree_passedEvents->Branch("pt_jetW2",             &pt_jetW2_t,            "pt_jetW2_t/F");
   tree_passedEvents->Branch("eta_jetW2",            &eta_jetW2_t,           "eta_jetW2_t/F");
   tree_passedEvents->Branch("qgl_jetW2",            &qgl_jetW2_t,           "qgl_jetW2_t/F");
   tree_passedEvents->Branch("pt_bJet",             &pt_bJet_t,            "pt_bJet_t/F");
   tree_passedEvents->Branch("eta_bJet",            &eta_bJet_t,           "eta_bJet_t/F");
   tree_passedEvents->Branch("pt_top",              &pt_top_t,             "pt_top_t/F");
   tree_passedEvents->Branch("eta_top",             &eta_top_t,            "eta_top_t/F");
   tree_passedEvents->Branch("nCentralJets",        &nCentralJets_t,       "nCentralJets_t/I");
   tree_passedEvents->Branch("nFwdJets",        &nFwdJets_t,       "nFwdJets_t/I");
   tree_passedEvents->Branch("hardestCentralJetPt", &hardestCentralJetPt_t,"hardestCentralJetPt_t/F");
   tree_passedEvents->Branch("deltaPhi_top_higgs",  &deltaPhi_top_higgs_t, "deltaPhi_top_higgs_t/F");
   tree_passedEvents->Branch("deltaEta_lept_qJet",  &deltaEta_lept_qJet_t, "deltaEta_lept_qJet_t/F");
   tree_passedEvents->Branch("deltaEta_bJet_qJet",  &deltaEta_bJet_qJet_t, "deltaEta_bJet_qJet_t/F");
   tree_passedEvents->Branch("mt_top",  &mt_top_t, "mt_top_t/F");
   tree_passedEvents->Branch("m_top",  &m_top_t, "m_top_t/F");
   tree_passedEvents->Branch("mt_W",  &mt_W_t, "mt_W_t/F");
   tree_passedEvents->Branch("m_W",  &m_W_t, "m_W_t/F");
   tree_passedEvents->Branch("BDT_lept",  &BDT_lept_t, "BDT_lept_t/F");
   tree_passedEvents->Branch("BDT2_lept",  &BDT2_lept_t, "BDT2_lept_t/F");
   tree_passedEvents->Branch("BDTold_lept",  &BDTold_lept_t, "BDTold_lept_t/F");
   tree_passedEvents->Branch("LD_lept",  &LD_lept_t, "LD_lept_t/F");
   tree_passedEvents->Branch("deltaR_lept_phot",  &deltaR_lept_phot_t, "deltaR_lept_phot_t/F");
   tree_passedEvents->Branch("m_ele_phot",  &m_ele_phot_t, "m_ele_phot_t/F");
   tree_passedEvents->Branch("zeppen",  &zeppen_t, "zeppen_t/F");
   tree_passedEvents->Branch("deltaR_top_phot_min_hadr",  &deltaR_top_phot_min_hadr_t, "deltaR_top_phot_min_hadr_t/F");
   tree_passedEvents->Branch("deltaR_qJet_phot_min_hadr",  &deltaR_qJet_phot_min_hadr_t, "deltaR_qJet_phot_min_hadr_t/F");



   tree_passedEvents->Branch("cosThetaStar",  &cosThetaStar_t, "cosThetaStar/F");

   tree_passedEvents->Branch("cosThetaStar",  &cosThetaStar_t, "cosThetaStar/F");
   tree_passedEvents->Branch("cosThetaStar_tH",  &cosThetaStar_tH_t, "cosThetaStar_tH/F");

   tree_passedEvents->Branch("pt_q",   &pt_q,  "pt_q/F");
   tree_passedEvents->Branch("eta_q",  &eta_q, "eta_q/F");
   tree_passedEvents->Branch("phi_q",  &phi_q, "phi_q/F");
   tree_passedEvents->Branch("e_q",    &e_q,   "e_q/F");
 
   tree_passedEvents->Branch("pt_h",   &pt_h,  "pt_h/F");
   tree_passedEvents->Branch("eta_h",  &eta_h, "eta_h/F");
   tree_passedEvents->Branch("phi_h",  &phi_h, "phi_h/F");
   tree_passedEvents->Branch("e_h",    &e_h,   "e_h/F");

   tree_passedEvents->Branch("pt_t",   &pt_t,  "pt_t/F");
   tree_passedEvents->Branch("eta_t",  &eta_t, "eta_t/F");
   tree_passedEvents->Branch("phi_t",  &phi_t, "phi_t/F");
   tree_passedEvents->Branch("e_t",    &e_t,   "e_t/F");

   tree_passedEvents->Branch("pt_b",  &pt_b, "pt_b/F");
   tree_passedEvents->Branch("eta_b",  &eta_b, "eta_b/F");
   tree_passedEvents->Branch("phi_b",  &phi_b, "phi_b/F");
   tree_passedEvents->Branch("e_b",  &e_b, "e_b/F");

   tree_passedEvents->Branch("pt_Wq",  &pt_Wq, "pt_Wq/F");
   tree_passedEvents->Branch("eta_Wq",  &eta_Wq, "eta_Wq/F");
   tree_passedEvents->Branch("phi_Wq",  &phi_Wq, "phi_Wq/F");
   tree_passedEvents->Branch("e_Wq",  &e_Wq, "e_Wq/F");

   tree_passedEvents->Branch("pt_Wqbar",  &pt_Wqbar, "pt_Wqbar/F");
   tree_passedEvents->Branch("eta_Wqbar",  &eta_Wqbar, "eta_Wqbar/F");
   tree_passedEvents->Branch("phi_Wqbar",  &phi_Wqbar, "phi_Wqbar/F");
   tree_passedEvents->Branch("e_Wqbar",  &e_Wqbar, "e_Wqbar/F");

   tree_passedEvents->Branch("Ht",  &Ht_t, "Ht_t/F");
   tree_passedEvents->Branch("mHt",  &mHt_t, "mHt_t/F");
   tree_passedEvents->Branch("mSt",  &mSt_t, "mSt_t/F");



   TFile* file_thq_hadronic_weights = TFile::Open("thq_hadronic_puweights.root");
   TH1F* h1_thq_hadronic_weights = (TH1F*)file_thq_hadronic_weights->Get("puweights");
 
//   std::string qglFileName = "/afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1.root";
   //std::string qglFileName="/afs/cern.ch/user/m/micheli/public/ttH/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1.root";
   //std::string qglFileName="/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root";
   std::string qglFileName="/afs/cern.ch/work/p/pandolf/public/ReducedHisto_2012.root";

   QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator( qglFileName );

   //std::string thqFile_LD = "/afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/UserCode/pandolf/VHgg/VHggAnalysis/finalizedTrees_THq_feasibility/THq_tHqLeptonic_mH125_8TeV_testtest_presel_isLeptonic_CSV.root";
   //std::string tthFile_LD = "/afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/UserCode/pandolf/VHgg/VHggAnalysis/finalizedTrees_THq_feasibility/THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_isLeptonic_CSV.root";
   //THqLeptonicLikelihoodCalculator * thqlikeli = new THqLeptonicLikelihoodCalculator( thqFile_LD, tthFile_LD );

   std::string thqFile_LD = "TMVA/THq_tHqLeptonic_mH125_8TeV_testtest_presel_isLeptonic_CSV.root";
   std::string tthFile_LD = "TMVA/THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_isLeptonic_CSV.root";
   THqLeptonicLikelihoodCalculator * thqlikeli = new THqLeptonicLikelihoodCalculator( thqFile_LD, tthFile_LD );


   BTagSFUtil* btsfutil = new BTagSFUtil(bTaggerType_, 13);


     TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 
     TMVA::Reader *reader2 = new TMVA::Reader( "!Color:!Silent" ); 
     TMVA::Reader *readerold = new TMVA::Reader( "!Color:!Silent" ); 
  
     float njets_bdt;
     float nbjets_loose_bdt;
     float deltaPhi_top_higgs_bdt;
     float mt_top_bdt;
     float eta_qJet_bdt;
     float charge_lept_bdt;
     float pt_bJet_bdt;
     float deltaEta_lept_qJet_bdt;
  
  
     reader->AddVariable( "njets", &njets_bdt );
     //reader->AddVariable( "deltaPhi_top_higgs", &deltaPhi_top_higgs_bdt );
     reader->AddVariable( "mt_top", &mt_top_bdt );
     reader->AddVariable( "eta_qJet", &eta_qJet_bdt );
     reader->AddVariable( "charge_lept", &charge_lept_bdt );
     reader->AddVariable( "deltaEta_lept_qJet", &deltaEta_lept_qJet_bdt );
  
  
     reader2->AddVariable( "njets", &njets_bdt );
     reader2->AddVariable( "nbjets_loose", &nbjets_loose_bdt );
     reader2->AddVariable( "deltaPhi_top_higgs", &deltaPhi_top_higgs_bdt );
     reader2->AddVariable( "mt_top", &mt_top_bdt );
     reader2->AddVariable( "eta_qJet", &eta_qJet_bdt );
     reader2->AddVariable( "charge_lept", &charge_lept_bdt );
     reader2->AddVariable( "deltaEta_lept_qJet", &deltaEta_lept_qJet_bdt );
     reader2->AddVariable( "pt_bJet", &pt_bJet_bdt );
  
     readerold->AddVariable( "njets", &njets_bdt );
     readerold->AddVariable( "nbjets_loose", &nbjets_loose_bdt );
     readerold->AddVariable( "deltaPhi_top_higgs", &deltaPhi_top_higgs_bdt );
     readerold->AddVariable( "mt_top", &mt_top_bdt );
     readerold->AddVariable( "eta_qJet", &eta_qJet_bdt );
     readerold->AddVariable( "charge_lept", &charge_lept_bdt );
     readerold->AddVariable( "pt_bJet", &pt_bJet_bdt );

  
     TString methodName = "BDTG method";
     //TString weightfile = "TMVA/weights/prova_nonjets_BDTG.weights.xml";
     TString weightfileNEW = "TMVA/weights/provaNEW_BDTG.weights.xml";
     TString weightfile2 = "TMVA/weights/provaNEW3_BDTG.weights.xml";
     TString weightfileold = "TMVA/weights/prova_BDTG.weights.xml";
     
  
     std::cout << "-> Booking BDT" << std::endl;
     reader->BookMVA( methodName, weightfileNEW ); 
     reader2->BookMVA( methodName, weightfile2 ); 
     readerold->BookMVA( methodName, weightfileold ); 



   if (tree_ == 0) return;

   Long64_t nentries = tree_->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   if( xSection_!=1. && xSection_!=0. ) {
     std::cout << std::endl;
     std::cout << "-> Cross-Section: " << xSection_ << " pb" << std::endl;
     std::cout << "-> # Generated Events: " << nGenEvents_ << std::endl;
     std::cout << "-> Event Weight: " << xSection_/nGenEvents_ << std::endl << std::endl;
   }


   h1_nGenEvents->SetBinContent(1, nGenEvents_);









   float N_all_lept = 0.;
   float N_all_hadr = 0.;

   float N_qMatched_lept = 0.;
   float N_tMatched_lept = 0.;
   float N_hMatched_lept = 0.;
   float N_bMatched_lept = 0.;
   float N_WqMatched_lept = 0.;
   float N_WqbarMatched_lept = 0.;

   float N_qMatched_hadr = 0.;
   float N_tMatched_hadr = 0.;
   float N_hMatched_hadr = 0.;
   float N_bMatched_hadr = 0.;
   float N_WqMatched_hadr = 0.;
   float N_WqbarMatched_hadr = 0.;






   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if( ientry % 20000 == 0 ) std::cout << "-> Entry " << ientry << " / " << nentries << std::endl;

      nb = tree_->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      bool debug = (event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999);

      bool isMC = ( run<5 );


      if( isMC ) {

        if( nGenEvents_<=0 ) {
          std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
          exit(99);
        }

        eventWeight = xSection_ / nGenEvents_ ;
 

        h1_nvertex->Fill( nvtx, eventWeight );

        // pu/pt reeventWeighting
        if( dataset_ == "tHqHadronic_mH125_8TeV_testtest" )
          eventWeight *= h1_thq_hadronic_weights->GetBinContent( nvtx );

        if(dopureeventWeight_) eventWeight *= pu_weight; 
        if(doptreeventWeight_) eventWeight *= pt_weight; 


      }


      h1_nvertex_PUW->Fill( nvtx, eventWeight );



      if( useGenPhotons_ ) {
        ptphot1 = gen_pt_gamma1;
        ptphot2 = gen_pt_gamma2;
      }


      TLorentzVector diphot;
      if( useGenPhotons_ ) {
        TLorentzVector genphot1, genphot2;
        genphot1.SetPtEtaPhiM( gen_pt_gamma1, gen_eta_gamma1, gen_phi_gamma1, 0. );
        genphot2.SetPtEtaPhiM( gen_pt_gamma2, gen_eta_gamma2, gen_phi_gamma2, 0. );
        diphot = genphot1 + genphot2;
      } else {
        diphot.SetPtEtaPhiM( ptggnewvtx, etagg, phigg, massggnewvtx);
      }



      if( debug ) {
        std::cout << std::endl << std::endl << "#### DEBUG LOG FOR EVENT " << DEBUG_EVENT_NUMBER_ << std::endl << std::endl;
        std::cout << std::endl << "Here are the photons: " << std::endl;
        std::cout << "   Phot1:  pt: " << ptphot1 << " eta: " << etascphot1 << std::endl;
        std::cout << "   Phot2:  pt: " << ptphot2 << " eta: " << etascphot2 << std::endl;
        std::cout << std::endl << "Diphoton pt: " << diphot.Pt() << " mass: " << diphot.M() << std::endl;
      }

      TLorentzVector phot1, phot2;
      phot1.SetPtEtaPhiM( ptphot1, etaphot1, phiphot1, 0. );
      phot2.SetPtEtaPhiM( ptphot2, etaphot2, phiphot2, 0. );


      if( debug ) {
        std::cout << "Checking deltaRToTrack: " << std::endl;
        std::cout << "deltaRToTrackphot1: " << deltaRToTrackphot1 << std::endl;
        std::cout << "deltaRToTrackphot2: " << deltaRToTrackphot2 << std::endl;
      }


      if( deltaRToTrackphot1<1. || deltaRToTrackphot2<1. ) continue;

      if( debug ) {
        std::cout << "ptele1: " << ptele1 << " etaele1: " << etaele1 << std::endl;
        std::cout << "ptele2: " << ptele2 << " etaele2: " << etaele2 << std::endl;
        std::cout << "ptmu1: " << ptmu1 << " etamu1: " << etamu1 << std::endl;
        std::cout << "ptmu2: " << ptmu2 << " etamu2: " << etamu2 << std::endl;
      }
      
      // tag lepton right away
      TLorentzVector lept;
      std::vector<TLorentzVector> electrons, muons;
      //if( ptelenontr901>ptLept_thresh_ ) {
      //  TLorentzVector l;
      //  l.SetPtEtaPhiE( ptelenontr901, etaelenontr901, phielenontr901, eneelenontr901 );
      //  electrons.push_back(l);
      //}
      //if( ptelenontr902>ptLept_thresh_ ) {
      //  TLorentzVector l;
      //  l.SetPtEtaPhiE( ptelenontr902, etaelenontr902, phielenontr902, eneelenontr902 );
      //  electrons.push_back(l);
      //}
      if( ptele1>ptLept_thresh_ ) {
        TLorentzVector l;
        l.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
        if( l.DeltaR(phot1)>0.5 && l.DeltaR(phot2)>0.5 )
          electrons.push_back(l);
      }
      if( ptele2>ptLept_thresh_ ) {
        TLorentzVector l;
        l.SetPtEtaPhiE( ptele2, etaele2, phiele2, eneele2 );
        if( l.DeltaR(phot1)>0.5 && l.DeltaR(phot2)>0.5 )
          electrons.push_back(l);
      }
      if( ptmu1>ptLept_thresh_ ) {
        TLorentzVector l;
        l.SetPtEtaPhiE( ptmu1, etamu1, phimu1, enemu1 );
        muons.push_back(l);
      }
      if( ptmu2>ptLept_thresh_ ) {
        TLorentzVector l;
        l.SetPtEtaPhiE( ptmu2, etamu2, phimu2, enemu2 );
        muons.push_back(l);
      }
      //if( ptmu1>ptLept_thresh_ ) {
      //  TLorentzVector l;
      //  l.SetPtEtaPhiE( ptmu1, etamu1, phimu1, enemu1 );
      //  muons.push_back(l);
      //}
      //if( ptmu2>ptLept_thresh_ ) {
      //  TLorentzVector l;
      //  l.SetPtEtaPhiE( ptmu2, etamu2, phimu2, enemu2 );
      //  muons.push_back(l);
      //}



      if( electrons.size() + muons.size() > 1 ) {

        continue;

      } else if( electrons.size() > 0 ) {

        //if( electrons[0].DeltaR( phot1 ) < 0.5 || electrons[0].DeltaR( phot2 ) < 0.5 ) continue;

        isLeptonic_t = true;
        isMu_t = false;
        lept = electrons[0];
        charge_lept_t = chargeele1;
        //charge_lept_t = chargeelenontr901;

      } else if( muons.size() > 0 ) {

        isLeptonic_t = true;
        isMu_t = true;
        lept = muons[0];
        charge_lept_t = chargemu1;

      } else {

        isLeptonic_t = false;
        isMu_t = false;
        charge_lept_t = 0;

      }

      
      
//      float leptVetoPt_thresh=5.;
//      if( (ptele1>leptVetoPt_thresh && ptmu1>leptVetoPt_thresh) || ptele2>leptVetoPt_thresh || ptmu2>leptVetoPt_thresh ) continue;
//
//        float deltaR1 = lept.DeltaR(phot1);
//        //if( deltaR1<1. ) continue;
//        float deltaR2 = lept.DeltaR(phot2);
//        //if( deltaR2<1. ) continue;
//
//        h1_deltaR_lept_phot1->Fill( deltaR1, eventWeight );
//        h1_deltaR_lept_phot2->Fill( deltaR2, eventWeight );
//        deltaR_lept_phot_t = TMath::Min( deltaR1, deltaR2 );
//        if( deltaR_lept_phot_t < 0.5 ) continue;
//
//
//      isLeptonic_t = (ptele1>ptLept_thresh_ || ptmu1>ptLept_thresh_);//cut on pt lepton
//      isMu_t=ptmu1>ptLept_thresh_;
//      if( isLeptonic_t ) {
//        if( isMu_t )
//          lept.SetPtEtaPhiE( ptmu1, etamu1, phimu1, enemu1 );
//        else
//          lept.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
//      }


      if( selectionType_=="presel_isLeptonic" && !isLeptonic_t ) continue;


      // analysis cuts:

      //if(npu>=30) continue;

      // if(massggnewvtx<90 || massggnewvtx>190) continue;
      //      if(diphot.M()<100 || diphot.M()>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
         || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance


      if( debug ) {
        std::cout << "-> Passed photon eta acceptance." << std::endl;
      }
      //     if(ptphot1<ptphot1cut) continue; //pt first photon

      bool ebeb = TMath::Abs(etascphot1)<1.4442 && TMath::Abs(etascphot2)<1.4442;

      float triggerThreshPhot1=33.;
      float triggerThreshPhot2=25.;

      if(ptphot1 < triggerThreshPhot1) continue; //pt first photon
      if(ptphot2 < triggerThreshPhot2) continue; //pt second photon


      if( debug ) {
        std::cout << "-> Passed photon trigger pt thresholds." << std::endl;
        std::cout << std::endl << "Checking analysis thresholds: " << std::endl;
        std::cout << ptphot1 << " > " << ptphot1cut_ << "*" << massggnewvtx << "/120. = " << ptphot1cut_*massggnewvtx/120. << " ? " << std::endl;
        std::cout << ptphot2 << " > " << ptphot2cut_ << "*" << massggnewvtx << "/120. = " << ptphot2cut_*massggnewvtx/120. << " ? " << std::endl;

      }

      if( debug ) {
        std::cout << "-> Passed photon analysis pt thresholds." << std::endl;
      }

      //if(ptgg_thresh_>0 && diphot.Pt()< pthiggsmincut_) continue; //pt higgs min



       // photon identification

       // r9 categories
       bool isr9phot1(0), isr9phot2(0);
 
       if(TMath::Abs(etascphot1)<1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)<1.4442 && r9phot2>.94) isr9phot2 = 1;
       if(TMath::Abs(etascphot1)>1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)>1.4442 && r9phot2>.94) isr9phot2 = 1;
 
 
       // photon id
       bool idphot1(0), idphot2(0), looseidphot1(0), looseidphot2(0), pxlphot1(1), pxlphot2(1);
 

      if(!invert_photonCuts_){       
       idphot1 = (idcicpfphot1 >= photonID_thresh_);
       idphot2 = (idcicpfphot2 >= photonID_thresh_);
      }else{
	idphot1 = (idcicpfphot1 >= photonID_thresh_);
	idphot2 = (idcicpfphot2 < photonID_thresh_);
      }


      if( !useGenPhotons_ ) {

        if(!cs_){ // photon id no control sample
   
          if(photonID_thresh_>0) {
            if(!(idphot1)) continue;
            if(!(idphot2)) continue;
          }else{
            if(!(idphot1 && pxlphot1)) continue;
            if(!(idphot2 && pxlphot2)) continue;
          }
        
        }else{ // photon id for control sample
         
          looseidphot1 = (idcicpfphot1 > 0 );
          looseidphot2 = (idcicpfphot2 > 0 );
          //	  if( !( (idphot1 && looseidphot2 && !idphot2) || (idphot2 && looseidphot1 && !idphot1) ) ) continue; 
          // Not perfect should be using the same electronVeto wrt CiC selection (now using matchedPromptEle veto)
          if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 
   
        }

      } // if !useGenPhotons_


      if( debug ) {
        std::cout << "-> Passed photon ID." << std::endl;
      }

      if(diphot.M()<100. || diphot.M()>180.) continue;
      bool blindCondition = !isMC && BLIND_ && massggnewvtx>115. && massggnewvtx<135.;
      if( blindCondition ) continue;




      if( isLeptonic_t )
        h1_cutFlow_lept->Fill(0., eventWeight);
      else
        h1_cutFlow_hadr->Fill(0., eventWeight);

      if( debug ) {
        std::cout << "-> Passed mgg 100-180 cut." << std::endl;
      }


      if(ptphot1 < ptphot1cut_*massggnewvtx/120.) continue; //pt first photon (running)
      if(ptphot2 < ptphot2cut_ ) continue; //pt second photon (not running)

      if( isLeptonic_t )
        h1_cutFlow_lept->Fill(1., eventWeight);
      else
        h1_cutFlow_hadr->Fill(1., eventWeight);


      if( debug ) {
        if( isLeptonic_t ) {
          std::cout << "-> Lepton: pt: " << lept.Pt() << " eta: " << lept.Eta() << " charge: " << charge_lept_t << std::endl;
          if( isMu_t ) std::cout << "It's a muon" << std::endl;
          else std::cout << "It's an electron" << std::endl;
        } else {
          std::cout << "-> No leptons: this is a hadronic event." << std::endl;
        }
      }
      
      
      
      h1_mgg_prepresel->Fill( massggnewvtx, eventWeight );
      

      
      
      // jets
      int njets_selected = 0;
      int njets_selected_btagloose = 0;
      int njets_selected_btagmedium = 0;
      std::vector<int> index_MatchedJet;
      std::vector<int> index_selected;
      std::vector<int> index_selected_btagloose;
      std::vector<int> index_selected_btagmedium;
      Ht_t = 0.;
      
      if( debug ) std::cout << "jets: " << std::endl;


      for( unsigned ijet=0; ijet<njets; ++ijet ) {

        if( ptcorrjet[ijet] < ptjetthresh_count_ ) continue;
        if( fabs(etajet[ijet]) > etajetthresh_count_ ) continue;
        

        //PU ID
        bool passedPUID = true;
        if(use_PUID_){

          if(TMath::Abs(etajet[ijet]) > 4.7) passedPUID = 0;

          if(TMath::Abs(etajet[ijet]) < 2.5) {
            if(betastarjet[ijet] > 0.2 * log( nvtx - PUID_betastar_thresh_ ) ) passedPUID = 0;
            if(rmsjet[ijet] > 0.06) passedPUID = 0;
          } else if(TMath::Abs(etajet[ijet]) < 2.75){
            if(betastarjet[ijet] > 0.3 * log( nvtx - PUID_betastar_thresh_ ) ) passedPUID = 0;
            if(rmsjet[ijet] > 0.05) passedPUID = 0;
          } else if(TMath::Abs(etajet[ijet]) < 3.){
            if(rmsjet[ijet] > 0.05) passedPUID = 0;
          } else {
            if(rmsjet[ijet] > 0.055) passedPUID = 0;
          }
        }

      
        if( isLeptonic_t ) {
          TLorentzVector thisJet;
          thisJet.SetPtEtaPhiE(ptcorrjet[ijet],etajet[ijet],phijet[ijet],ecorrjet[ijet]);
          if( thisJet.DeltaR(lept)<0.5 ) continue;
        }
      

        if( debug ) {
           std::cout << ijet << "/" << njets << " pt: " << ptcorrjet[ijet] << " eta: " << etajet[ijet] << " btagcsvjet: " << btagcsvjet[ijet] << std::endl;
        }
      
        bool btagged_loose = false;
        bool btagged_medium = false;
      
        if(bTaggerType_=="JP") {
          btagged_loose = btagjprobjet[ijet]>0.275;
          btagged_medium = btagjprobjet[ijet]>0.545;
        } else if(bTaggerType_=="CSV") {
          btagged_loose = btagcsvjet[ijet]>0.244;
          btagged_medium = btagcsvjet[ijet]>0.679;
        } else {
          std::cout<<"WARNING: btag type "<<bTaggerType_<<"not found"<<std::endl;
        }
      
        //// then modify btags with Scale Factors:
        //if( isMC )
        //  btsfutil->modifyBTagsWithSF_fast(btagged_loose, btagged_medium, ptcorrjet[ijet], etajet[ijet], partPdgIDjet[ijet], "mean" );
      
      
        if( btagged_loose ) {
          njets_selected_btagloose++;
          index_selected_btagloose.push_back(ijet);
        }

        if( btagged_medium ) {
          njets_selected_btagmedium++;
          index_selected_btagmedium.push_back(ijet);
        }
      
      
        index_selected.push_back(ijet);
      
      
        njets_selected++;
      
      
        //AnalysisJet thisJet;
        //thisJet.set
        //selectedJets.push_back( thisJet );
      
      } //for jets
      


      if( index_selected_btagmedium.size()==0 ) { 
        if( debug ) std::cout << "didn't find bJet" << std::endl;
        continue;
      }

      if( njets_selected_btagmedium>nbtagmedium_upper_thresh_ ) continue;


      if( isLeptonic_t )
        h1_cutFlow_lept->Fill(2., eventWeight);
      else
        h1_cutFlow_hadr->Fill(2., eventWeight);



      int bJet_index = index_selected_btagmedium[0];
      TLorentzVector bJet;
      bJet.SetPtEtaPhiE(ptcorrjet[bJet_index],etajet[bJet_index],phijet[bJet_index],ecorrjet[bJet_index]);

      if( debug ) std::cout << "bJet: pt: " << bJet.Pt() << " eta: " << bJet.Eta() << std::endl;


      //// veto events with 2 leptons:
      //float leptVetoPt_thresh=5.;
      //TLorentzVector lept;

      //// *** TWO MUONS
      //// (this is bad)
      //if( ptmu2>leptVetoPt_thresh ) { 

      //  continue; 
      //

      //// *** ONE ELECTRON AND ONE MUON
      //// (this is bad only if electron is not a photon)
      //} else if( ptmu1>leptVetoPt_thresh && ptele1>leptVetoPt_thresh ) { // electron and muon is bad only if electron is not a photon

      //  TLorentzVector electron;
      //  electron.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
      //  if( electron.DeltaR(phot1)>0.1 && electron.DeltaR(phot2)>0.1 ) continue;
      //

      //// *** TWO ELECTRONS
      //// (this is bad only if both or neither electrons are photons)
      //} else if( ptele1>leptVetoPt_thresh && ptele2>leptVetoPt_thresh ) { // two electrons is bad only if both different from photons

      //  TLorentzVector ele1, ele2;
      //  ele1.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
      //  ele2.SetPtEtaPhiE( ptele2, etaele2, phiele2, eneele2 );
      //  bool ele1isphot = (ele1.DeltaR(phot1)<0.1 || ele1.DeltaR(phot2)<0.1);
      //  bool ele2isphot = (ele2.DeltaR(phot2)<0.1 || ele2.DeltaR(phot2)<0.1);
      //  if( !ele1isphot && !ele2isphot ) { // two leptons
      //    continue; 
      //  } else if( ele1isphot && ele2isphot ) { //zero leptons
      //    continue; 
      //  } else if( !ele1isphot && ele2isphot ) {
      //    if( ptele1>ptLept_thresh_ ) {
      //      lept.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
      //      isLeptonic_t = true;
      //      isMu_t = false;
      //    } else { //there was one which was not a photon but didnt pass the pt cut
      //      continue;
      //    }
      //  } else if( ele1isphot && !ele2isphot ) {
      //    if( ptele1>ptLept_thresh_ ) {
      //      lept.SetPtEtaPhiE( ptele2, etaele2, phiele2, eneele2 );
      //      isLeptonic_t = true;
      //      isMu_t = false;
      //    } else { //there was one which was not a photon but didnt pass the pt cut
      //      continue;
      //    }
      //  } else { //this shouldnt be possible
      //    continue;
      //  }


      //// *** AND NOW THE ONE LEPTON CASE (MUON OR ELECTRON)
      //} else if( ptele1>ptLept_thresh_ && !(ptemu1>leptVetoPt_thresh) ) {

      //  lept.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );
      //  isLeptonic_t = true;
      //  isMu_t = false;

      //} else if( !(ptele1>leptVetoPt_thresh) && ptemu1>ptLept_thresh_ ) {

      //  lept.SetPtEtaPhiE( ptmu1, etamu1, phimu1, enemu1 );
      //  isLeptonic_t = true;
      //  isMu_t = true;

      //} else {

      //  continue; //what else is there?

      //}



      int index_jetW1=-1;
      int index_jetW2=-1;
      

      TLorentzVector top, top2;
      TLorentzVector neutrino;
      TLorentzVector jetW1,jetW2;
  

      // now look for q-jet:
      int index_qJet=-1;
      TLorentzVector qJet;
      nCentralJets_t = 0;
      nFwdJets_t = 0;
      float hardestCentralJetPt=0.;
      float eta_thresh_qJet = (isLeptonic_t) ? 1. : eta_qJet_thresh_hadr_;
      for( unsigned ii=0; ii<index_selected.size(); ++ii ) {

        int i = index_selected[ii];

        if( i==index_selected_btagmedium[0] ) continue;

        if( ptcorrjet[i]<ptjetthresh_count_ ) continue;
        if( fabs(etajet[i])<eta_thresh_qJet ) {
          if( nCentralJets_t==0 ) hardestCentralJetPt = ptcorrjet[i]; //hardest central jet
          nCentralJets_t++;
        } else {
          if( index_qJet<0 ) { //hardest forward jet
            index_qJet=i;
            qJet.SetPtEtaPhiE(ptcorrjet[i],etajet[i],phijet[i],ecorrjet[i]);
          } else {
            nFwdJets_t++;
          }
        }

      }
      if( index_qJet<0 ) continue;

      if( debug ) std::cout << "qJet: pt: " << qJet.Pt() << " eta: " << qJet.Eta() << std::endl;



      deltaEta_bJet_qJet_t = fabs( bJet.Eta()- qJet.Eta());


      if( isLeptonic_t )
        h1_cutFlow_lept->Fill(3., eventWeight);
      else
        h1_cutFlow_hadr->Fill(3., eventWeight);





      // ------------------------------------
      // ---                              ---
      // ---   ** LEPTONIC CHANNEL **     ---
      // ---                              ---
      // ------------------------------------

      if( isLeptonic_t ) { 

        deltaEta_lept_qJet_t = fabs( lept.Eta() - qJet.Eta() );

        if( fabs(qJet.Eta())<1. ) continue;
        if( njets_selected > njets_upper_thresh_lept_ ) continue;
        if( njets_selected < njets_thresh_lept_ ) continue;
        if( nCentralJets_t > nCentralJets_upper_thresh_lept_ ) continue;

        h1_cutFlow_lept->Fill(4., eventWeight);

        
        if( deltaEta_lept_qJet_t < deltaEta_lept_qJet_thresh_ ) continue;


        zeppen_t = -999.;



        float px_neutr = epfMet*cos(phipfMet);
        float py_neutr = epfMet*sin(phipfMet);
        neutrino.SetPtEtaPhiE( epfMet,0,phipfMet,epfMet );

        top = bJet + lept + neutrino;

        float pz_neutr = getNeutrinoPz( lept, px_neutr, py_neutr );
        neutrino.SetPxPyPzE( px_neutr, py_neutr, pz_neutr, epfMet );
        top2 = bJet + lept + neutrino;
        if( debug ) {
          std::cout << "bJet: pt: " << bJet.Pt() << " eta: " << bJet.Eta() << std::endl;
          std::cout << "qJet: pt: " << qJet.Pt() << " eta: " << qJet.Eta() << std::endl;
          std::cout << "neutrino: pt: " << neutrino.Pt() << " phi: " << neutrino.Phi() << std::endl;
          std::cout << "top: pt: " << top.Pt() << " eta: " << top.Eta() << std::endl;
        }

        


        if( !isMu_t ) {
 
          TLorentzVector ele_phot1 = lept + phot1;
          TLorentzVector ele_phot2 = lept + phot2;

          if( fabs(ele_phot1.M()-91.19) < fabs(ele_phot2.M()-91.19) ) {
            h1_m_ele_phot->Fill( ele_phot1.M(), eventWeight );
            m_ele_phot_t = ele_phot1.M();
          } else {
            h1_m_ele_phot->Fill( ele_phot2.M(), eventWeight );
            m_ele_phot_t = ele_phot2.M();
          }

        }


        // compute costhetastar:
        TLorentzVector lept_topF(lept);
        lept_topF.Boost(-top.BoostVector());
        TLorentzVector qJet_topF(qJet);
        qJet_topF.Boost(-top.BoostVector());
        cosThetaStar_t = cos(lept_topF.Angle(qJet_topF.Vect()));


        TLorentzVector th = top + diphot;

        TLorentzVector lept_thF(lept);
        lept_thF.Boost(-th.BoostVector());
        TLorentzVector qJet_thF(qJet);
        qJet_thF.Boost(-th.BoostVector());
        cosThetaStar_tH_t = cos(lept_thF.Angle(qJet_thF.Vect()));



      // ------------------------------------
      // ---                              ---
      // ---   ** HADRONIC CHANNEL **     ---
      // ---                              ---
      // ------------------------------------

      } else { 


        if( njets_selected < njets_thresh_hadr_ ) continue;


        // choose top selection method:
        bool top_max_deltaphi = true;
        bool top_best_mass = false;
        bool top_hardest = false;
 
        float topmass = 172.5;

        if( top_max_deltaphi ) {

          float maxDeltaPhi=0;
          
          // define top candidate as the one maximising deltaPhi(top, higgs):
          for( unsigned ii=0; ii<index_selected.size(); ++ii ) {
            for( unsigned jj=ii+1; jj<index_selected.size(); ++jj ) {

              int i = index_selected[ii];
              int j = index_selected[jj];

              if( i==index_selected_btagmedium[0] || j==index_selected_btagmedium[0] ) continue;
              if( i==index_qJet || j==index_qJet ) continue;


              TLorentzVector jetW1_tmp,jetW2_tmp;
              jetW1_tmp.SetPtEtaPhiE(ptcorrjet[i],etajet[i],phijet[i],ecorrjet[i]);
              jetW2_tmp.SetPtEtaPhiE(ptcorrjet[j],etajet[j],phijet[j],ecorrjet[j]);

              TLorentzVector top_tmp = jetW1_tmp+jetW2_tmp+bJet;

              float thisDeltaPhi = fabs(top_tmp.DeltaPhi(diphot));
              if( thisDeltaPhi>maxDeltaPhi ) {
                maxDeltaPhi=thisDeltaPhi;
                index_jetW1 = i;
                index_jetW2 = j;
              }

            } //j
          } //i

        } else if( top_best_mass ) {

          float best_m_top = -9999.;
          
          // define top candidate as the one with mass closest to nominal
          for( unsigned ii=0; ii<index_selected.size(); ++ii ) {
            for( unsigned jj=ii+1; jj<index_selected.size(); ++jj ) {

              int i = index_selected[ii];
              int j = index_selected[jj];

              if( i==index_selected_btagmedium[0] || j==index_selected_btagmedium[0] ) continue;
              if( i==index_qJet || j==index_qJet ) continue;


              TLorentzVector jetW1_tmp,jetW2_tmp;
              jetW1_tmp.SetPtEtaPhiE(ptcorrjet[i],etajet[i],phijet[i],ecorrjet[i]);
              jetW2_tmp.SetPtEtaPhiE(ptcorrjet[j],etajet[j],phijet[j],ecorrjet[j]);

              TLorentzVector top_tmp = jetW1_tmp+jetW2_tmp+bJet;

              float this_m_top = top.M();
              if( fabs(this_m_top-topmass)<fabs(best_m_top-topmass) ) {
                best_m_top=this_m_top;
                index_jetW1 = i;
                index_jetW2 = j;
              }

            } //j
          } //i

        } else if( top_hardest ) {

          
          // define top candidate as the one with mass closest to nominal
          for( unsigned ii=0; ii<index_selected.size(); ++ii ) {
            for( unsigned jj=ii+1; jj<index_selected.size(); ++jj ) {

              int i = index_selected[ii];
              int j = index_selected[jj];

              if( i==index_selected_btagmedium[0] || j==index_selected_btagmedium[0] ) continue;
              if( i==index_qJet || j==index_qJet ) continue;

              index_jetW1 = i;
              index_jetW2 = j;
              break;

            } //j
          } //i

        }


        if( index_jetW1<0 || index_jetW2<0 ) continue;

        jetW1.SetPtEtaPhiE(ptcorrjet[index_jetW1],etajet[index_jetW1],phijet[index_jetW1],ecorrjet[index_jetW1]);
        jetW2.SetPtEtaPhiE(ptcorrjet[index_jetW2],etajet[index_jetW2],phijet[index_jetW2],ecorrjet[index_jetW2]);

        if( fabs(jetW1.Eta())<2. ) nCentralJets_t--;
        if( fabs(jetW2.Eta())<2. ) nCentralJets_t--;
        
        if( fabs(jetW1.Eta())>2. ) nFwdJets_t--;
        if( fabs(jetW2.Eta())>2. ) nFwdJets_t--;
        

        //for( unsigned ii=0; ii<index_selected.size(); ++ii ) {
      
        //  int i = index_selected[ii];
      
        //  if( i==index_selected_btagmedium[0] ) continue;
        //  if( i==index_jetW1 ) continue;
        //  if( i==index_jetW2 ) continue;
      
        //  if( ptcorrjet[i]<20. ) continue;
        //  if( fabs(etajet[i])<eta_thresh_qJet ) {
        //    if( nCentralJets==0 ) hardestCentralJetPt = ptcorrjet[i]; //hardest central jet
        //    nCentralJets++;
        //  } else if( index_qJet<0 ) { //hardest forward jet
        //    index_qJet=i;
        //    qJet.SetPtEtaPhiE(ptcorrjet[i],etajet[i],phijet[i],ecorrjet[i]);
        //  }
      
        //}
        //if( index_qJet<0 ) continue;


        if( qJet.Pt() < pt_qJet_thresh_hadr_ ) continue;

        TLorentzVector W = jetW1 + jetW2;

        top = bJet + W;

        float Wmass = 80.4;

        if( fabs(top.M()-topmass)>m_top_thresh_hadr_ ) continue;
        if( fabs(W.M()-Wmass)>m_W_thresh_hadr_ ) continue;


        zeppen_t = diphot.Eta() - 0.5*( qJet.Eta() + top.Eta() );


        //if( nCentralJetsHadr_t > nCentralJets_upper_thresh_hadr_ ) continue;

        if( deltaEta_bJet_qJet_t<deltaEta_bJet_qJet_thresh_hadr_ ) continue;


        float deltaR_top_phot1 = top.DeltaR(phot1);
        float deltaR_top_phot2 = top.DeltaR(phot2);
        deltaR_top_phot_min_hadr_t = TMath::Min( deltaR_top_phot1, deltaR_top_phot2 );

        float deltaR_qJet_phot1 = qJet.DeltaR(phot1);
        float deltaR_qJet_phot2 = qJet.DeltaR(phot2);
        deltaR_qJet_phot_min_hadr_t = TMath::Min( deltaR_qJet_phot1, deltaR_qJet_phot2 );

        // compute costhetastar:
        // first boost bjet and qjet in top mass frame
        TLorentzVector bJet_topF(bJet);
        bJet_topF.Boost(-top.BoostVector());
        TLorentzVector qJet_topF(qJet);
        qJet_topF.Boost(-top.BoostVector());
        cosThetaStar_t = cos(bJet_topF.Angle(qJet_topF.Vect()));

        TLorentzVector th = top + diphot;

        TLorentzVector bJet_thF(bJet);
        bJet_thF.Boost(-th.BoostVector());
        TLorentzVector qJet_thF(qJet);
        qJet_thF.Boost(-th.BoostVector());
        cosThetaStar_tH_t = cos(bJet_thF.Angle(qJet_thF.Vect()));

        Ht_t = bJet.Pt() + qJet.Pt() + jetW1.Pt() + jetW2.Pt();

        TLorentzVector totalHardEvent = top + qJet + diphot;
        mHt_t = totalHardEvent.Pt();

        TLorentzVector higgs_top = top + diphot;
        mSt_t = higgs_top.Pt();
        

      } //if is hadronic
      
      


            



      
      if( isSignalMC && use_thqPartonInfo_ ) {

        TLorentzVector h, t, q, b, Wq, Wqbar;
        h.SetPtEtaPhiE( pt_h, eta_h, phi_h, e_h );
        t.SetPtEtaPhiE( pt_t, eta_t, phi_t, e_t );
        b.SetPtEtaPhiE( pt_b, eta_b, phi_b, e_b );
        q.SetPtEtaPhiE( pt_q, eta_q, phi_q, e_q );
        Wq.SetPtEtaPhiE( pt_Wq, eta_Wq, phi_Wq, e_Wq );
        Wqbar.SetPtEtaPhiE( pt_Wqbar, eta_Wqbar, phi_Wqbar, e_Wqbar );


        if( isLeptonic_t )
          N_all_lept+=eventWeight;
        else
          N_all_hadr+=eventWeight;

        float deltaR_q = q.DeltaR(qJet);
        float deltaR_t = t.DeltaR(top);
        float deltaPhi_t = fabs(t.DeltaPhi(top));
        float deltaR_b = b.DeltaR(bJet);
        float deltaR_h = h.DeltaR(diphot);

        if( deltaR_q<0.5 && isLeptonic_t )
          N_qMatched_lept+=eventWeight;
        if( deltaPhi_t<0.5 && isLeptonic_t )
          N_tMatched_lept+=eventWeight;
        if( deltaR_b<0.5 && isLeptonic_t )
          N_bMatched_lept+=eventWeight;
        if( deltaR_h<0.5 && isLeptonic_t )
          N_hMatched_lept+=eventWeight;

        if( deltaR_q<0.5 && !isLeptonic_t )
          N_qMatched_hadr+=eventWeight;
        if( deltaPhi_t<0.5 && !isLeptonic_t )
          N_tMatched_hadr+=eventWeight;
        if( deltaR_b<0.5 && !isLeptonic_t )
          N_bMatched_hadr+=eventWeight;
        if( deltaR_h<0.5 && !isLeptonic_t )
          N_hMatched_hadr+=eventWeight;


        //if( Wq.DeltaR(qJet)<0.5 )
        //  N_WqMatched+=eventWeight
        //if( Wq.DeltaR(qJet)<0.5 )
        //  N_qMatched+=eventWeight

        h1_deltaR_q->  Fill( deltaR_q, eventWeight );
        h1_deltaR_t->  Fill( deltaR_t, eventWeight );
        h1_deltaPhi_t->Fill( deltaPhi_t, eventWeight );
        h1_deltaR_b->  Fill( deltaR_b, eventWeight );
        h1_deltaR_h->  Fill( deltaR_h, eventWeight );

      } // if isMC



      
      if(njets_selected<njets_thresh_) continue;


      if( debug ) {
        std::cout << "-> Passed njets cut." << std::endl;
      }

      if(njets_selected>njets_upper_thresh_)continue;
      if( debug ) {
        std::cout << "-> Passed njets upper cut." << std::endl;
      }

      
      
      // set tree vars:
      float deltaPhi_top_higgs = fabs(top.DeltaPhi(diphot));

      njets_t  = njets_selected;
      nbjets_loose_t = njets_selected_btagloose;
      nbjets_medium_t = njets_selected_btagmedium;
      ptPhot1_t = ptphot1;
      ptPhot2_t = ptphot2;
      ptRunPhot1_t = ptphot1*120./massggnewvtx;
      ptRunPhot2_t = ptphot2*120./massggnewvtx;
      etaPhot1_t = etaphot1;
      etaPhot2_t = etaphot2;
      mgg_t = diphot.M();

      ptgg_t = diphot.Pt();
      etagg_t = diphot.Eta();
      ptRungg_t = diphot.Pt()*120./massggnewvtx;
      
      pt_lept_t = (isLeptonic_t) ? lept.Pt() : 0.;
      eta_lept_t = (isLeptonic_t) ? lept.Eta() : 0.;

      pt_qJet_t = qJet.Pt();
      eta_qJet_t = qJet.Eta();
      float qgl_qJet = qglikeli->computeQGLikelihood2012( qJet.Pt(), qJet.Eta(), rhoPF, nChg_QCjet[index_qJet]+nNeutral_ptCutjet[index_qJet], ptD_QCjet[index_qJet], axis2_QCjet[index_qJet] );
      qgl_qJet_t = qgl_qJet;
      

      if( !isLeptonic_t ) {

        pt_jetW1_t = jetW1.Pt();
        eta_jetW1_t = jetW1.Eta();
        float qgl_jetW1 = qglikeli->computeQGLikelihood2012( jetW1.Pt(), jetW1.Eta(), rhoPF, nChg_QCjet[index_jetW1]+nNeutral_ptCutjet[index_jetW1], ptD_QCjet[index_jetW1], axis2_QCjet[index_jetW1] );
        qgl_jetW1_t = qgl_jetW1;
        
        pt_jetW2_t = jetW2.Pt();
        eta_jetW2_t = jetW2.Eta();
        float qgl_jetW2 = qglikeli->computeQGLikelihood2012( jetW2.Pt(), jetW2.Eta(), rhoPF, nChg_QCjet[index_jetW2]+nNeutral_ptCutjet[index_jetW2], ptD_QCjet[index_jetW2], axis2_QCjet[index_jetW2] );
        qgl_jetW2_t = qgl_jetW2;

      } else {
      
        pt_jetW1_t = 0.;
        eta_jetW1_t = -10.;
        qgl_jetW1_t = -1.;
        
        pt_jetW2_t = 0.;
        eta_jetW2_t = -10.;
        qgl_jetW2_t = -1.;

      }

      pt_bJet_t = bJet.Pt();
      eta_bJet_t = bJet.Eta();
      
      pt_top_t = top.Pt();
      eta_top_t = top.Eta();
      
      hardestCentralJetPt_t = hardestCentralJetPt;

      deltaPhi_top_higgs_t = deltaPhi_top_higgs;
      
      mt_top_t = top.Mt();
      m_top_t = top2.M();

      TLorentzVector W = (isLeptonic_t) ? (lept+neutrino) : (jetW1+jetW2);

      mt_W_t = W.Mt();
      m_W_t = W.M();
 
      njets_bdt              = njets_t;
      nbjets_loose_bdt       = nbjets_loose_t;
      deltaPhi_top_higgs_bdt = deltaPhi_top_higgs_t;
      mt_top_bdt             = mt_top_t;
      eta_qJet_bdt           = eta_qJet_t;
      charge_lept_bdt        = charge_lept_t;
      pt_bJet_bdt            = pt_bJet_t;
      deltaEta_lept_qJet_bdt = deltaEta_lept_qJet_t;
  
      BDT_lept_t = reader->EvaluateMVA( "BDTG method" );
      BDT2_lept_t = reader2->EvaluateMVA( "BDTG method" );
      BDTold_lept_t = readerold->EvaluateMVA( "BDTG method" );

      // choose old BDT:
      //BDT_lept_t = readerold->EvaluateMVA( "BDTG method" );

      LD_lept_t = thqlikeli->computeLikelihood( njets_t, eta_qJet_t, mt_top_t, charge_lept_t, deltaEta_lept_qJet_t );
      if( debug ) {
        std::cout << "njets_t: " <<  njets_t << std::endl;
        std::cout << "eta_qJet_t: " <<  eta_qJet_t << std::endl;
        std::cout << "mt_top_t: " <<  mt_top_t << std::endl;
        std::cout << "charge_lept_t: " <<  charge_lept_t << std::endl;
        std::cout << "deltaEta_lept_qJet_t: " <<  deltaEta_lept_qJet_t << std::endl;
      }

      //if( BDT_lept_t < bdt_lept_thresh_ ) continue; 



      
      // and now fill histos:      
      
      h1_ptphot0->Fill(ptphot1, eventWeight);
      h1_ptphot1->Fill(ptphot2, eventWeight);
      
      h1_ptrunphot0->Fill(ptphot1*120/massggnewvtx, eventWeight);
      h1_ptrunphot1->Fill(ptphot2*120/massggnewvtx, eventWeight);
      
      h1_etaphot0->Fill(etaphot1, eventWeight);
      h1_etaphot1->Fill(etaphot2, eventWeight);
      
      
      h1_njets->Fill( njets_selected, eventWeight );
      h1_nbjets_loose->Fill( njets_selected_btagloose, eventWeight );
      h1_nbjets_medium->Fill( njets_selected_btagmedium, eventWeight );
      
      
      
      h1_mgg_presel->Fill( massggnewvtx, eventWeight );
      h1_njets_presel->Fill( njets_selected, eventWeight ); 
      
      

      
      h1_pt_qJet->Fill( qJet.Pt(), eventWeight ); 
      h1_eta_qJet->Fill( qJet.Eta(), eventWeight ); 
      h1_qgl_qJet->Fill( qgl_qJet, eventWeight ); 
      
      h1_pt_bJet->Fill( bJet.Pt(), eventWeight ); 
      h1_eta_bJet->Fill( bJet.Eta(), eventWeight ); 
      
      h1_nCentralJets->Fill( nCentralJets_t, eventWeight );
      h1_hardestCentralJetPt->Fill( hardestCentralJetPt, eventWeight );
      
      
      h1_Ht->Fill(Ht_t,eventWeight);
      
      
      
      h1_ptDiphot->Fill( diphot.Pt(), eventWeight );
      h1_ptRunDiphot->Fill( diphot.Pt()*120./massggnewvtx, eventWeight );
      
      
      
      
      h1_mgg->Fill( massggnewvtx, eventWeight );
      
      
      
      if( isLeptonic_t ) {

        h1_BDTlept_lept->Fill( BDT_lept_t, eventWeight );
        h1_LDlept_lept->Fill( LD_lept_t, eventWeight );

        //leptons
        h1_pt_lept->Fill(lept.Pt(),eventWeight);
        h1_eta_lept->Fill(lept.Eta(),eventWeight);
      
        //MET
        h1_pfMet->Fill(epfMet,eventWeight);
        float deltaPhi_met_bJet = fabs(bJet.DeltaPhi(neutrino));
        h1_deltaPhi_met_bJet->Fill(deltaPhi_met_bJet,eventWeight);

        h1_mgg_lept->Fill( massggnewvtx, eventWeight );
        if( BDT_lept_t>bdt_lept_thresh_ ) {
          h1_mgg_lept_BDT->Fill( massggnewvtx, eventWeight );
        }
//std::cout << "*****************************************************************************           run: " << run << " event: " << event << " mgg: " << massggnewvtx << " LD: " << LD_lept_t << " charge: " << charge_lept_t << std::endl;
        if( LD_lept_t>ld_lept_thresh_ ) {
          h1_cutFlow_lept->Fill(5., eventWeight);
          h1_mgg_lept_LD->Fill( massggnewvtx, eventWeight );
        }
        

      } else { // hadronic channel

        h1_mgg_hadr->Fill( massggnewvtx, eventWeight );
        if( nCentralJets_t < nCentralJets_upper_thresh_hadr_ ) {
          h1_mgg_hadr_centralJetVeto->Fill( massggnewvtx, eventWeight );
        }

      }
      
      h1_deltaPhi_top_higgs->Fill(deltaPhi_top_higgs,eventWeight);
     

      


      tree_passedEvents->Fill();


   } //for entries


   if( isSignalMC ) {
     std::cout << std::endl;
     std::cout << "LEPTONIC CHANNEL:" << std::endl;
     std::cout << "-> Matched Higgs in " << (float)N_hMatched_lept*100./N_all_lept << "% of the cases." << std::endl;
     std::cout << "-> Matched bJet  in " << (float)N_bMatched_lept*100./N_all_lept << "% of the cases." << std::endl;
     std::cout << "-> Matched top   in " << (float)N_tMatched_lept*100./N_all_lept << "% of the cases." << std::endl;
     std::cout << "-> Matched qJet  in " << (float)N_qMatched_lept*100./N_all_lept << "% of the cases." << std::endl;

     std::cout << std::endl;
     std::cout << "HADRONIC CHANNEL:" << std::endl;
     std::cout << "-> Matched Higgs in " << (float)N_hMatched_hadr*100./N_all_hadr << "% of the cases." << std::endl;
     std::cout << "-> Matched bJet  in " << (float)N_bMatched_hadr*100./N_all_hadr << "% of the cases." << std::endl;
     std::cout << "-> Matched top   in " << (float)N_tMatched_hadr*100./N_all_hadr << "% of the cases." << std::endl;
     std::cout << "-> Matched qJet  in " << (float)N_qMatched_hadr*100./N_all_hadr << "% of the cases." << std::endl;
   }


   outFile_->cd();


   tree_passedEvents->Write();





   h1_nGenEvents->Write();

   h1_cutFlow_lept->Write();
   h1_cutFlow_hadr->Write();

   h1_nvertex->Write();
   h1_nvertex_PUW->Write();

   h1_njets_presel->Write();
   h1_njets->Write();
   h1_nbjets_loose->Write();
   h1_nbjets_medium->Write();

   h1_ptphot0->Write();
   h1_ptphot1->Write();
   h1_ptrunphot0->Write();
   h1_ptrunphot1->Write();
   h1_etaphot0->Write();
   h1_etaphot1->Write();


   h1_pfMet->Write();
   h1_Ht->Write();

   h1_pt_qJet->Write();
   h1_eta_qJet->Write();
   h1_qgl_qJet->Write();
   
   h1_pt_bJet->Write();
   h1_eta_bJet->Write();

   h1_pt_top->Write();
   h1_eta_top->Write();

   h1_nCentralJets->Write();
   h1_hardestCentralJetPt->Write();

   h1_ptDiphot->Write();
   h1_ptRunDiphot->Write();


   h1_mgg_prepresel->Write();
   h1_mgg_presel->Write();
   h1_mgg->Write();
   h1_mgg_lept->Write();
   h1_mgg_lept_BDT->Write();
   h1_mgg_lept_LD->Write();
   h1_mgg_hadr->Write();
   h1_mgg_hadr_centralJetVeto->Write();


   h1_pt_lept->Write();
   h1_eta_lept->Write();

   h1_deltaR_lept_phot1->Write();
   h1_deltaR_lept_phot2->Write();
   h1_m_ele_phot->Write();

   h1_deltaPhi_met_bJet->Write();
   h1_deltaPhi_top_higgs->Write();

   h1_deltaR_q->Write();
   h1_deltaR_t->Write();
   h1_deltaR_h->Write();
   h1_deltaR_b->Write();
   h1_deltaR_Wq->Write();
   h1_deltaR_Wqbar->Write();
   h1_deltaPhi_t->Write();

   h1_BDTlept_lept->Write();
   h1_LDlept_lept->Write();

} //finalize


















Int_t RedNtpFinalizer_THq::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!tree_) return 0;
   return tree_->GetEntry(entry);
}
Long64_t RedNtpFinalizer_THq::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!tree_) return -5;
   Long64_t centry = tree_->LoadTree(entry);
   if (centry < 0) return centry;
   if (tree_->GetTreeNumber() != fCurrent) {
      fCurrent = tree_->GetTreeNumber();
   }
   return centry;
}

void RedNtpFinalizer_THq::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   fCurrent = -1;
   tree_->SetMakeClass(1);

   tree_->SetBranchAddress("run", &run, &b_run);
   tree_->SetBranchAddress("event", &event, &b_event);
   tree_->SetBranchAddress("lumi", &lumi, &b_lumi);
   tree_->SetBranchAddress("H_event", &H_event, &b_H_event);
   tree_->SetBranchAddress("WH_event", &WH_event, &b_WH_event);
   tree_->SetBranchAddress("ZH_event", &ZH_event, &b_ZH_event);
   tree_->SetBranchAddress("V_event", &V_event, &b_V_event);
   tree_->SetBranchAddress("Vqq_event", &Vqq_event, &b_Vqq_event);
   tree_->SetBranchAddress("Zbb_event", &Zbb_event, &b_Zbb_event);
   tree_->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   tree_->SetBranchAddress("massgg", &massgg, &b_massgg);
   tree_->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   tree_->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
   tree_->SetBranchAddress("phigg", &phigg, &b_phigg);
   tree_->SetBranchAddress("etagg", &etagg, &b_etagg);
   tree_->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
   tree_->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
   tree_->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
   tree_->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
   tree_->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
   tree_->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
   tree_->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
   tree_->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
   tree_->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
   tree_->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
   tree_->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
   tree_->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
   tree_->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
   tree_->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
   tree_->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
   tree_->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
   tree_->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
   tree_->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
   tree_->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
   tree_->SetBranchAddress("energyErrphot1", &energyErrphot1, &b_energyErrphot1);
   tree_->SetBranchAddress("energyErrphot2", &energyErrphot2, &b_energyErrphot2);
   tree_->SetBranchAddress("energySmearingphot1", &energySmearingphot1, &b_energySmearingphot1);
   tree_->SetBranchAddress("energySmearingphot2", &energySmearingphot2, &b_energySmearingphot2);
   tree_->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
   tree_->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
   tree_->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
   tree_->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
   tree_->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
   tree_->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
   tree_->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
   tree_->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
   tree_->SetBranchAddress("idmvaphot1", &idmvaphot1, &b_idmvaphot1);
   tree_->SetBranchAddress("idmvaphot2", &idmvaphot2, &b_idmvaphot2);
   tree_->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
   tree_->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
   tree_->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
   tree_->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
   tree_->SetBranchAddress("idcicpfphot1", &idcicpfphot1, &b_idcicpfphot1);
   tree_->SetBranchAddress("idcicpfphot2", &idcicpfphot2, &b_idcicpfphot2);
   tree_->SetBranchAddress("idcicpfnoelvetophot1", &idcicpfnoelvetophot1, &b_idcicpfnoelvetophot1);
   tree_->SetBranchAddress("idcicpfnoelvetophot2", &idcicpfnoelvetophot2, &b_idcicpfnoelvetophot2);
   tree_->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
   tree_->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
   tree_->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
   tree_->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
   tree_->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
   tree_->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
   tree_->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
   tree_->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
   tree_->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
   tree_->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
   tree_->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
   tree_->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
   tree_->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
   tree_->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
   tree_->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
   tree_->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
   tree_->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
   tree_->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
   tree_->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
   tree_->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
   tree_->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
   tree_->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
   tree_->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
   tree_->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
   tree_->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
   tree_->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
   tree_->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
   tree_->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
   tree_->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
   tree_->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
   tree_->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
   tree_->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
   tree_->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
   tree_->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
   tree_->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
   tree_->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
   tree_->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
   tree_->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
   tree_->SetBranchAddress("njets", &njets, &b_njets);
   tree_->SetBranchAddress("ptjet", ptjet, &b_ptjet);
   tree_->SetBranchAddress("ptcorrjet", ptcorrjet, &b_ptcorrjet);
   tree_->SetBranchAddress("ecorrjet", ecorrjet, &b_ecorrjet);
   tree_->SetBranchAddress("etajet", etajet, &b_etajet);
   tree_->SetBranchAddress("phijet", phijet, &b_phijet);
   tree_->SetBranchAddress("betajet", betajet, &b_betajet);
   tree_->SetBranchAddress("betastarjet", betastarjet, &b_betastarjet);
   tree_->SetBranchAddress("btagvtxjet", btagvtxjet, &b_btagvtxjet);
   tree_->SetBranchAddress("btagcsvjet", btagcsvjet, &b_btagcsvjet);
   tree_->SetBranchAddress("btagtrkjet", btagtrkjet, &b_btagtrkjet);
   tree_->SetBranchAddress("btagjprobjet", btagjprobjet, &b_btagjprobjet);
   tree_->SetBranchAddress("ptDjet", ptDjet, &b_ptDjet);
   tree_->SetBranchAddress("ptD_QCjet", ptD_QCjet, &b_ptD_QCjet);
   tree_->SetBranchAddress("axis2_QCjet", axis2_QCjet, &b_axis2_QCjet);
   tree_->SetBranchAddress("rmsjet", rmsjet, &b_rmsjet);
   tree_->SetBranchAddress("ntrkjet", ntrkjet, &b_ntrkjet);
   tree_->SetBranchAddress("nneutjet", nneutjet, &b_nneutjet);
   tree_->SetBranchAddress("nChg_QCjet", nChg_QCjet, &b_nChg_QCjet);
   tree_->SetBranchAddress("nNeutral_ptCutjet", nNeutral_ptCutjet, &b_nNeutral_ptCutjet);
   tree_->SetBranchAddress("jetIdSimple_mvajet", jetIdSimple_mvajet, &b_jetIdSimple_mvajet);
   tree_->SetBranchAddress("jetIdFull_mvajet", jetIdFull_mvajet, &b_jetIdFull_mvajet);
   tree_->SetBranchAddress("jetId_dR2Meanjet", jetId_dR2Meanjet, &b_jetId_dR2Meanjet);
   tree_->SetBranchAddress("jetId_betaStarClassicjet", jetId_betaStarClassicjet, &b_jetId_betaStarClassicjet);
   tree_->SetBranchAddress("jetId_frac01jet", jetId_frac01jet, &b_jetId_frac01jet);
   tree_->SetBranchAddress("jetId_frac02jet", jetId_frac02jet, &b_jetId_frac02jet);
   tree_->SetBranchAddress("jetId_frac03jet", jetId_frac03jet, &b_jetId_frac03jet);
   tree_->SetBranchAddress("jetId_frac04jet", jetId_frac04jet, &b_jetId_frac04jet);
   tree_->SetBranchAddress("jetId_frac05jet", jetId_frac05jet, &b_jetId_frac05jet);
   tree_->SetBranchAddress("jetId_betajet", jetId_betajet, &b_jetId_betajet);
   tree_->SetBranchAddress("jetId_betaStarjet", jetId_betaStarjet, &b_jetId_betaStarjet);
   tree_->SetBranchAddress("jetIdCutBased_wpjet", jetIdCutBased_wpjet, &b_jetIdCutBased_wpjet);
   tree_->SetBranchAddress("jetIdSimple_wpjet", jetIdSimple_wpjet, &b_jetIdSimple_wpjet);
   tree_->SetBranchAddress("jetIdFull_wpjet", jetIdFull_wpjet, &b_jetIdFull_wpjet);
   tree_->SetBranchAddress("assjet", assjet, &b_assjet);
   tree_->SetBranchAddress("partPdgIDjet", partPdgIDjet, &b_partPdgIDjet);
   tree_->SetBranchAddress("partMomPdgIDjet", partMomPdgIDjet, &b_partMomPdgIDjet);
   tree_->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   tree_->SetBranchAddress("zeppenjet", &zeppenjet, &b_zeppenjet);
   tree_->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   tree_->SetBranchAddress("deltaphinewvtx", &deltaphinewvtx, &b_deltaphinewvtx);
   tree_->SetBranchAddress("deltaphigg", &deltaphigg, &b_deltaphigg);
   tree_->SetBranchAddress("invmassjet", &invmassjet, &b_invmassjet);
   tree_->SetBranchAddress("invmass2g1j", &invmass2g1j, &b_invmass2g1j);
   tree_->SetBranchAddress("invmass2g2j", &invmass2g2j, &b_invmass2g2j);
   tree_->SetBranchAddress("pt2g2j", &pt2g2j, &b_pt2g2j);
   tree_->SetBranchAddress("eta2j", &eta2j, &b_eta2j);
   tree_->SetBranchAddress("phi2j", &phi2j, &b_phi2j);
   tree_->SetBranchAddress("pt2j", &pt2j, &b_pt2j);
   tree_->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   tree_->SetBranchAddress("vtxId", &vtxId, &b_vtxId);
   tree_->SetBranchAddress("vtxPos_x", &vtxPos_x, &b_vtxPos_x);
   tree_->SetBranchAddress("vtxPos_y", &vtxPos_y, &b_vtxPos_y);
   tree_->SetBranchAddress("vtxPos_z", &vtxPos_z, &b_vtxPos_z);
   tree_->SetBranchAddress("vtxIdMVA", &vtxIdMVA, &b_vtxIdMVA);
   tree_->SetBranchAddress("vtxIdEvtProb", &vtxIdEvtProb, &b_vtxIdEvtProb);
   tree_->SetBranchAddress("sMet", &sMet, &b_sMet);
   tree_->SetBranchAddress("eMet", &eMet, &b_eMet);
   tree_->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   tree_->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   tree_->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
   tree_->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
   tree_->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
   tree_->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
   tree_->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
   tree_->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
   tree_->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
   tree_->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
   tree_->SetBranchAddress("eShiftedScaledMetPUcorr", &eShiftedScaledMetPUcorr, &b_eShiftedScaledMetPUcorr);
   tree_->SetBranchAddress("phiShiftedScaledMetPUcorr", &phiShiftedScaledMetPUcorr, &b_phiShiftedScaledMetPUcorr);
   tree_->SetBranchAddress("eSmearedShiftedMePUcorrt", &eSmearedShiftedMePUcorrt, &b_eSmearedShiftedMetPUcorr);
   tree_->SetBranchAddress("phiSmearedShiftedMetPUcorr", &phiSmearedShiftedMetPUcorr, &b_phiSmearedShiftedMetPUcorr);
   tree_->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
   tree_->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
   tree_->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   tree_->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   tree_->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
   tree_->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
   tree_->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   tree_->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   tree_->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
   tree_->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
   tree_->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   tree_->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   tree_->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
   tree_->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
   tree_->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   tree_->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   tree_->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
   tree_->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
   tree_->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
   tree_->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
   tree_->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
   tree_->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
   tree_->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
   tree_->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
   tree_->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
   tree_->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
   tree_->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
   tree_->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
   tree_->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
   tree_->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
   tree_->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
   tree_->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
   tree_->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
   tree_->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
   tree_->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
   tree_->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
   tree_->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
   tree_->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
   tree_->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
   tree_->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
   tree_->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
   tree_->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
   tree_->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
   tree_->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   tree_->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   tree_->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
   tree_->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
   tree_->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
   tree_->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
   tree_->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
   tree_->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
   tree_->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   tree_->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   tree_->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
   tree_->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
   tree_->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   tree_->SetBranchAddress("npu", &npu, &b_npu);
   tree_->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
   tree_->SetBranchAddress("xsection", &xsection, &b_xsection);
   tree_->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
   tree_->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
   tree_->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   tree_->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
   tree_->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
   tree_->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
   tree_->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
   tree_->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
   tree_->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
   tree_->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
   tree_->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
   tree_->SetBranchAddress("gen_pt_genjet1", &gen_pt_genjet1, &b_gen_pt_genjet1);
   tree_->SetBranchAddress("gen_pt_genjet2", &gen_pt_genjet2, &b_gen_pt_genjet2);
   tree_->SetBranchAddress("gen_eta_genjet1", &gen_eta_genjet1, &b_gen_eta_genjet1);
   tree_->SetBranchAddress("gen_eta_genjet2", &gen_eta_genjet2, &b_gen_eta_genjet2);
   tree_->SetBranchAddress("gen_phi_genjet1", &gen_phi_genjet1, &b_gen_phi_genjet1);
   tree_->SetBranchAddress("gen_phi_genjet2", &gen_phi_genjet2, &b_gen_phi_genjet2);
   tree_->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
   tree_->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
   tree_->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
   tree_->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
   tree_->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
   tree_->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
   tree_->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
   tree_->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
   tree_->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
   tree_->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
   tree_->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
   tree_->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
   tree_->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
   tree_->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
   tree_->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
   tree_->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
   tree_->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
   tree_->SetBranchAddress("chargeele1", &chargeele1, &b_chargeele1);
   tree_->SetBranchAddress("chargeele2", &chargeele2, &b_chargeele2);
   tree_->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
   tree_->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
   tree_->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
   tree_->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
   tree_->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
   tree_->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
   tree_->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
   tree_->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
   tree_->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
   tree_->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
   tree_->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
   tree_->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
   tree_->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
   tree_->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
   tree_->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
   tree_->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
   tree_->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
   tree_->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
   tree_->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
   tree_->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
   tree_->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
   tree_->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
   tree_->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
   tree_->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
   tree_->SetBranchAddress("oEmoPele1", &oEmoPele1, &b_oEmoPele1);
   tree_->SetBranchAddress("oEmoPele2", &oEmoPele2, &b_oEmoPele2);
   tree_->SetBranchAddress("mvanotrigele1", &mvanotrigele1, &b_mvanotrigele1);
   tree_->SetBranchAddress("mvanotrigele2", &mvanotrigele2, &b_mvanotrigele2);
   tree_->SetBranchAddress("mvatrigele1", &mvatrigele1, &b_mvatrigele1);
   tree_->SetBranchAddress("mvatrigele2", &mvatrigele2, &b_mvatrigele2);
   tree_->SetBranchAddress("matchconvele1", &matchconvele1, &b_matchconvele1);
   tree_->SetBranchAddress("matchconvele2", &matchconvele2, &b_matchconvele2);
   tree_->SetBranchAddress("chHadIso03ele1", &chHadIso03ele1, &b_chHadIso03ele1);
   tree_->SetBranchAddress("chHadIso03ele2", &chHadIso03ele2, &b_chHadIso03ele2);
   tree_->SetBranchAddress("nHadIso03ele1", &nHadIso03ele1, &b_nHadIso03ele1);
   tree_->SetBranchAddress("nHadIso03ele2", &nHadIso03ele2, &b_nHadIso03ele2);
   tree_->SetBranchAddress("photIso03ele1", &photIso03ele1, &b_photIso03ele1);
   tree_->SetBranchAddress("photIso03ele2", &photIso03ele2, &b_photIso03ele2);
   tree_->SetBranchAddress("isorelele1", &isorelele1, &b_isorelele1);
   tree_->SetBranchAddress("pteleloose1", &pteleloose1, &b_pteleloose1);
   tree_->SetBranchAddress("pteleloose2", &pteleloose2, &b_pteleloose2);
   tree_->SetBranchAddress("etaeleloose1", &etaeleloose1, &b_etaeleloose1);
   tree_->SetBranchAddress("etaeleloose2", &etaeleloose2, &b_etaeleloose2);
   tree_->SetBranchAddress("phieleloose1", &phieleloose1, &b_phieleloose1);
   tree_->SetBranchAddress("phieleloose2", &phieleloose2, &b_phieleloose2);
   tree_->SetBranchAddress("eneeleloose1", &eneeleloose1, &b_eneeleloose1);
   tree_->SetBranchAddress("eneeleloose2", &eneeleloose2, &b_eneeleloose2);
   tree_->SetBranchAddress("sIeIeeleloose1", &sIeIeeleloose1, &b_sIeIeeleloose1);
   tree_->SetBranchAddress("sIeIeeleloose2", &sIeIeeleloose2, &b_sIeIeeleloose2);
   tree_->SetBranchAddress("dphieleloose1", &dphieleloose1, &b_dphieleloose1);
   tree_->SetBranchAddress("dphieleloose2", &dphieleloose2, &b_dphieleloose2);
   tree_->SetBranchAddress("detaeleloose1", &detaeleloose1, &b_detaeleloose1);
   tree_->SetBranchAddress("detaeleloose2", &detaeleloose2, &b_detaeleloose2);
   tree_->SetBranchAddress("mhitseleloose1", &mhitseleloose1, &b_mhitseleloose1);
   tree_->SetBranchAddress("mhitseleloose2", &mhitseleloose2, &b_mhitseleloose2);
   tree_->SetBranchAddress("d0eleloose1", &d0eleloose1, &b_d0eleloose1);
   tree_->SetBranchAddress("d0eleloose2", &d0eleloose2, &b_d0eleloose2);
   tree_->SetBranchAddress("dzeleloose1", &dzeleloose1, &b_dzeleloose1);
   tree_->SetBranchAddress("dzeleloose2", &dzeleloose2, &b_dzeleloose2);
   tree_->SetBranchAddress("invMasseleloose1g1", &invMasseleloose1g1, &b_invMasseleloose1g1);
   tree_->SetBranchAddress("invMasseleloose1g2", &invMasseleloose1g2, &b_invMasseleloose1g2);
   tree_->SetBranchAddress("invMasseleloose2g1", &invMasseleloose2g1, &b_invMasseleloose2g1);
   tree_->SetBranchAddress("invMasseleloose2g2", &invMasseleloose2g2, &b_invMasseleloose2g2);
   tree_->SetBranchAddress("oEmoPeleloose1", &oEmoPeleloose1, &b_oEmoPeleloose1);
   tree_->SetBranchAddress("oEmoPeleloose2", &oEmoPeleloose2, &b_oEmoPeleloose2);
   tree_->SetBranchAddress("mvanotrigeleloose1", &mvanotrigeleloose1, &b_mvanotrigeleloose1);
   tree_->SetBranchAddress("mvanotrigeleloose2", &mvanotrigeleloose2, &b_mvanotrigeleloose2);
   tree_->SetBranchAddress("mvatrigeleloose1", &mvatrigeleloose1, &b_mvatrigeleloose1);
   tree_->SetBranchAddress("mvatrigeleloose2", &mvatrigeleloose2, &b_mvatrigeleloose2);
   tree_->SetBranchAddress("matchconveleloose1", &matchconveleloose1, &b_matchconveleloose1);
   tree_->SetBranchAddress("matchconveleloose2", &matchconveleloose2, &b_matchconveleloose2);
   tree_->SetBranchAddress("chHadIso03eleloose1", &chHadIso03eleloose1, &b_chHadIso03eleloose1);
   tree_->SetBranchAddress("chHadIso03eleloose2", &chHadIso03eleloose2, &b_chHadIso03eleloose2);
   tree_->SetBranchAddress("nHadIso03eleloose1", &nHadIso03eleloose1, &b_nHadIso03eleloose1);
   tree_->SetBranchAddress("nHadIso03eleloose2", &nHadIso03eleloose2, &b_nHadIso03eleloose2);
   tree_->SetBranchAddress("photIso03eleloose1", &photIso03eleloose1, &b_photIso03eleloose1);
   tree_->SetBranchAddress("photIso03eleloose2", &photIso03eleloose2, &b_photIso03eleloose2);
   tree_->SetBranchAddress("chargeelenontr901", &chargeelenontr901, &b_chargeelenontr901);
   tree_->SetBranchAddress("chargeelenontr902", &chargeelenontr902, &b_chargeelenontr902);
   tree_->SetBranchAddress("ptelenontr901", &ptelenontr901, &b_ptelenontr901);
   tree_->SetBranchAddress("ptelenontr902", &ptelenontr902, &b_ptelenontr902);
   tree_->SetBranchAddress("etaelenontr901", &etaelenontr901, &b_etaelenontr901);
   tree_->SetBranchAddress("etaelenontr902", &etaelenontr902, &b_etaelenontr902);
   tree_->SetBranchAddress("phielenontr901", &phielenontr901, &b_phielenontr901);
   tree_->SetBranchAddress("phielenontr902", &phielenontr902, &b_phielenontr902);
   tree_->SetBranchAddress("eneelenontr901", &eneelenontr901, &b_eneelenontr901);
   tree_->SetBranchAddress("eneelenontr902", &eneelenontr902, &b_eneelenontr902);
   tree_->SetBranchAddress("chargemu1", &chargemu1, &b_chargemu1);
   tree_->SetBranchAddress("chargemu2", &chargemu2, &b_chargemu2);
   tree_->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   tree_->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   tree_->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   tree_->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   tree_->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   tree_->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   tree_->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
   tree_->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
   tree_->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
   tree_->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
   tree_->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
   tree_->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
   tree_->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
   tree_->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
   tree_->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
   tree_->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
   tree_->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
   tree_->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
   tree_->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
   tree_->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
   tree_->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
   tree_->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
   tree_->SetBranchAddress("chHadmu1", &chHadmu1, &b_chHadmu1);
   tree_->SetBranchAddress("chHadmu2", &chHadmu2, &b_chHadmu2);
   tree_->SetBranchAddress("nHadmu1", &nHadmu1, &b_nHadmu1);
   tree_->SetBranchAddress("nHadmu2", &nHadmu2, &b_nHadmu2);
   tree_->SetBranchAddress("photmu1", &photmu1, &b_photmu1);
   tree_->SetBranchAddress("photmu2", &photmu2, &b_photmu2);
   tree_->SetBranchAddress("puptmu1", &puptmu1, &b_puptmu1);
   tree_->SetBranchAddress("puptmu2", &puptmu2, &b_puptmu2);
   tree_->SetBranchAddress("ptmuloose1", &ptmuloose1, &b_ptmuloose1);
   tree_->SetBranchAddress("ptmuloose2", &ptmuloose2, &b_ptmuloose2);
   tree_->SetBranchAddress("etamuloose1", &etamuloose1, &b_etamuloose1);
   tree_->SetBranchAddress("etamuloose2", &etamuloose2, &b_etamuloose2);
   tree_->SetBranchAddress("phimuloose1", &phimuloose1, &b_phimuloose1);
   tree_->SetBranchAddress("phimuloose2", &phimuloose2, &b_phimuloose2);
   tree_->SetBranchAddress("enemuloose1", &enemuloose1, &b_enemuloose1);
   tree_->SetBranchAddress("enemuloose2", &enemuloose2, &b_enemuloose2);
   tree_->SetBranchAddress("pixhitsmuloose1", &pixhitsmuloose1, &b_pixhitsmuloose1);
   tree_->SetBranchAddress("pixhitsmuloose2", &pixhitsmuloose2, &b_pixhitsmuloose2);
   tree_->SetBranchAddress("trkhitsmuloose1", &trkhitsmuloose1, &b_trkhitsmuloose1);
   tree_->SetBranchAddress("trkhitsmuloose2", &trkhitsmuloose2, &b_trkhitsmuloose2);
   tree_->SetBranchAddress("hitsmuloose1", &hitsmuloose1, &b_hitsmuloose1);
   tree_->SetBranchAddress("hitsmuloose2", &hitsmuloose2, &b_hitsmuloose2);
   tree_->SetBranchAddress("chi2muloose1", &chi2muloose1, &b_chi2muloose1);
   tree_->SetBranchAddress("chi2muloose2", &chi2muloose2, &b_chi2muloose2);
   tree_->SetBranchAddress("matchmuloose1", &matchmuloose1, &b_matchmuloose1);
   tree_->SetBranchAddress("matchmuloose2", &matchmuloose2, &b_matchmuloose2);
   tree_->SetBranchAddress("d0muloose1", &d0muloose1, &b_d0muloose1);
   tree_->SetBranchAddress("d0muloose2", &d0muloose2, &b_d0muloose2);
   tree_->SetBranchAddress("dzmuloose1", &dzmuloose1, &b_dzmuloose1);
   tree_->SetBranchAddress("dzmuloose2", &dzmuloose2, &b_dzmuloose2);
   tree_->SetBranchAddress("chHadmuloose1", &chHadmuloose1, &b_chHadmuloose1);
   tree_->SetBranchAddress("chHadmuloose2", &chHadmuloose2, &b_chHadmuloose2);
   tree_->SetBranchAddress("nHadmuloose1", &nHadmuloose1, &b_nHadmuloose1);
   tree_->SetBranchAddress("nHadmuloose2", &nHadmuloose2, &b_nHadmuloose2);
   tree_->SetBranchAddress("photmuloose1", &photmuloose1, &b_photmuloose1);
   tree_->SetBranchAddress("photmuloose2", &photmuloose2, &b_photmuloose2);
   tree_->SetBranchAddress("puptmuloose1", &puptmuloose1, &b_puptmuloose1);
   tree_->SetBranchAddress("puptmuloose2", &puptmuloose2, &b_puptmuloose2);
   tree_->SetBranchAddress("hasPassedSinglePhot", &hasPassedSinglePhot, &b_hasPassedSinglePhot);
   tree_->SetBranchAddress("hasPassedDoublePhot", &hasPassedDoublePhot, &b_hasPassedDoublePhot);

   tree_->SetBranchAddress("pt_h", &pt_h, &b_pt_h);
   tree_->SetBranchAddress("eta_h", &eta_h, &b_eta_h);
   tree_->SetBranchAddress("phi_h", &phi_h, &b_phi_h);
   tree_->SetBranchAddress("e_h", &e_h, &b_e_h);

   tree_->SetBranchAddress("pt_t", &pt_t, &b_pt_t);
   tree_->SetBranchAddress("eta_t", &eta_t, &b_eta_t);
   tree_->SetBranchAddress("phi_t", &phi_t, &b_phi_t);
   tree_->SetBranchAddress("e_t", &e_t, &b_e_t);

   tree_->SetBranchAddress("pt_b", &pt_b, &b_pt_b);
   tree_->SetBranchAddress("eta_b", &eta_b, &b_eta_b);
   tree_->SetBranchAddress("phi_b", &phi_b, &b_phi_b);
   tree_->SetBranchAddress("e_b", &e_b, &b_e_b);

   tree_->SetBranchAddress("pt_q", &pt_q, &b_pt_q);
   tree_->SetBranchAddress("eta_q", &eta_q, &b_eta_q);
   tree_->SetBranchAddress("phi_q", &phi_q, &b_phi_q);
   tree_->SetBranchAddress("e_q", &e_q, &b_e_q);

   tree_->SetBranchAddress("pt_Wq", &pt_Wq, &b_pt_Wq);
   tree_->SetBranchAddress("eta_Wq", &eta_Wq, &b_eta_Wq);
   tree_->SetBranchAddress("phi_Wq", &phi_Wq, &b_phi_Wq);
   tree_->SetBranchAddress("e_Wq", &e_Wq, &b_e_Wq);

   tree_->SetBranchAddress("pt_Wqbar", &pt_Wqbar, &b_pt_Wqbar);
   tree_->SetBranchAddress("eta_Wqbar", &eta_Wqbar, &b_eta_Wqbar);
   tree_->SetBranchAddress("phi_Wqbar", &phi_Wqbar, &b_phi_Wqbar);
   tree_->SetBranchAddress("e_Wqbar", &e_Wqbar, &b_e_Wqbar);

}


void RedNtpFinalizer_THq::setDatasetName( const std::string& dataset ) {
  dataset_=dataset;
}



void RedNtpFinalizer_THq::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;


  // default values 
  
  dopureeventWeight_ = false;
  doptreeventWeight_ = true;
  r9cat_ = 1;
  photonID_thresh_ = 4;
  cs_ = false;
  useGenPhotons_ = false;


  ptphot1cut_ = 30.;
  ptphot2cut_ = 25.;

  ptjetthresh_count_ = 20.;
  etajetthresh_count_ = 5.;


  njets_thresh_=0;
  njets_upper_thresh_=1000;
  njets_thresh_lept_=0;
  njets_thresh_hadr_=0;
  njets_upper_thresh_lept_ = 1000;
  nbtagloose_thresh_=0;
  nbtagmedium_thresh_=0;
  nbtagmedium_upper_thresh_=1000;
  nCentralJets_upper_thresh_lept_=1000;

  bdt_lept_thresh_ = -10;
  ld_lept_thresh_ = -10;

  invert_photonCuts_=false;
  use_PUID_=true;

  njets_PUID_thresh_=0;
  ptLept_thresh_=10;

  PUID_betastar_thresh_=0.67;
  chooseBtaggedJets_=true;//be careful with this

  m_top_thresh_hadr_ = 10000.;
  m_W_thresh_hadr_ = 10000.;

  pt_qJet_thresh_hadr_ = 20.;
  eta_qJet_thresh_hadr_ = 2.;
  nCentralJets_upper_thresh_hadr_ = 1000;
  deltaEta_bJet_qJet_thresh_hadr_ = 0.;
  deltaEta_lept_qJet_thresh_ = 0.;


  if( selectionType=="presel" ) {

    // leave all cuts to default
    
  } else if( selectionType=="presel_isLeptonic" ) {

    // leave all cuts to default
    
  } else if( selectionType=="presel_plusPhotonCuts" ) {

    ptphot1cut_ = 50.;
    
  } else if( selectionType=="preselGEN" ) {

    useGenPhotons_ = true;
    
  } else if( selectionType=="presel_1CSVM" ) {

    nbtagmedium_upper_thresh_ = 1;
    
  } else if( selectionType=="sel0" ) {

    njets_thresh_lept_ = 2;
    njets_upper_thresh_lept_ = 3;
    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;
    
  } else if( selectionType=="sel1" ) {

    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    //njets_thresh_lept_ = 3;
    njets_thresh_hadr_ = 4;

    nbtagmedium_upper_thresh_ = 1;
    //nCentralJets_upper_thresh_lept_ = 1;
    bdt_lept_thresh_ = 0.2;

    m_top_thresh_hadr_ = 50.;
    //m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 40.;
    nCentralJets_upper_thresh_hadr_ = 1; 

  } else if( selectionType=="sel2" ) {

    njets_thresh_lept_ = 3;
    njets_thresh_hadr_ = 4;
    nbtagmedium_upper_thresh_ = 1;
    nCentralJets_upper_thresh_lept_ = 1;
    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    m_top_thresh_hadr_ = 50.;
    m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 40.;
    nCentralJets_upper_thresh_hadr_ = 1; 

  } else if( selectionType=="sel3" ) {

    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    //nbtagmedium_upper_thresh_ = 1;

    bdt_lept_thresh_ = 0.2;

    njets_thresh_hadr_ = 4;
    m_top_thresh_hadr_ = 50.;
    //m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 40.;
    nCentralJets_upper_thresh_hadr_ = 1; 


  } else if( selectionType=="sel4" ) { //this is the flagship selection

    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    //deltaEta_lept_qJet_thresh_ = 1.;

    bdt_lept_thresh_ = -0.4;
    ld_lept_thresh_ = 0.25;

    njets_thresh_hadr_ = 4;
    m_top_thresh_hadr_ = 40.;
    //m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 45.;
    nCentralJets_upper_thresh_hadr_ = 1; 


  } else if( selectionType=="sel4_pt25" ) { //like sel4 but with pt(jets)>25

    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    ptjetthresh_count_ = 25.;
    //deltaEta_lept_qJet_thresh_ = 1.;

    bdt_lept_thresh_ = -0.4;
    ld_lept_thresh_ = 0.25;

    njets_thresh_hadr_ = 4;
    m_top_thresh_hadr_ = 40.;
    //m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 45.;
    nCentralJets_upper_thresh_hadr_ = 1; 


  } else if( selectionType=="sel4_2" ) { //like sel4 but with lept CJV back in

    ptphot1cut_ = 50.;
    ptphot2cut_ = 25.;

    //deltaEta_lept_qJet_thresh_ = 1.;

    nCentralJets_upper_thresh_lept_ = 1;
    bdt_lept_thresh_ = 0.2;

    njets_thresh_hadr_ = 4;
    m_top_thresh_hadr_ = 50.;
    //m_W_thresh_hadr_ = 30.;
    pt_qJet_thresh_hadr_ = 40.;
    nCentralJets_upper_thresh_hadr_ = 1; 



  } else {
    std::cout << std::endl << std::endl << "Selection '" << selectionType << "' currently not implemented. Exiting." << std::endl;
    exit(12345);

  }


}




float RedNtpFinalizer_THq::getNeutrinoPz( TLorentzVector lepton, float pxPFMet, float pyPFMet) {

  
  float mW = 80.399;
  float app = lepton.E()*lepton.E() + pxPFMet*pxPFMet + pyPFMet*pyPFMet - lepton.Px()*lepton.Px() + pxPFMet*pxPFMet - lepton.Py()*lepton.Py() + pyPFMet*pyPFMet - lepton.Pz()*lepton.Pz() - mW*mW;
  float a= lepton.E()*lepton.E() - lepton.Pz()*lepton.Pz();
  float b= lepton.Pz()*app;
  float c= ( pxPFMet*pxPFMet + pyPFMet*pyPFMet )*lepton.E()*lepton.E() - app*app/4.;

  float pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  float pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
  if ( b*b-4.*a*c < 0. ){ 
    pznp=-b/(2.*a); 
    pznm =-b/(2.*a); 
  }
  

  float pn;
  if( fabs(pznp) < fabs(pznm) ){
    pn=pznp;
  } else{ 
    pn=pznm;
  }

  return pn;

}

