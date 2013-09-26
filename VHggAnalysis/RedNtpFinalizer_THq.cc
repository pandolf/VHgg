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
   tree_passedEvents->SetName("thq_m125_8TeV");
   //tree_passedEvents->SetName("tree_passedEvents");
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

   int category =11;
   bool isLep_ele;
   tree_passedEvents->Branch("category",        &category, "category/I");
   tree_passedEvents->Branch("isLep_ele",        &isLep_ele, "isLep_ele/O");
   tree_passedEvents->Branch("ptele1",        &ptele1, "ptele1/F");
   tree_passedEvents->Branch("etaele1",       &etaele1, "etaele1/F");
   tree_passedEvents->Branch("phiele1",       &phiele1, "phiele1/F");
   tree_passedEvents->Branch("sIeIeele1",     &sIeIeele1, "sIeIeele1/F");
   tree_passedEvents->Branch("dphiele1",      &dphiele1, "dphiele1/F");
   tree_passedEvents->Branch("detaele1",      &detaele1, "detaele1/F");
   tree_passedEvents->Branch("mhitsele1",     &mhitsele1, "mhitsele1/I");
   tree_passedEvents->Branch("d0ele1",        &d0ele1, "d0ele1/F");
   tree_passedEvents->Branch("dzele1",        &dzele1, "dzele1/F");
   tree_passedEvents->Branch("oEmoPele1",     &oEmoPele1, "oEmoPele1/F");
   tree_passedEvents->Branch("matchconvele1", &matchconvele1, "matchconvele1/I");
   tree_passedEvents->Branch("isorelele1",    &isorelele1, "isorelele1/F");


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

if(
!( event!=     1510 
&& event!=     1566 
&& event!=     1578 
&& event!=     1599 
&& event!=     1638 
&& event!=     1688 
&& event!=     1700 
&& event!=     1719 
&& event!=     1739 
&& event!=     1751 
&& event!=     1780 
&& event!=     1783 
&& event!=     1787 
&& event!=     1798 
&& event!=     1810 
&& event!=     1851 
&& event!=     1853 
&& event!=     1855 
&& event!=     1858 
&& event!=     1871 
&& event!=     1902 
&& event!=     1907 
&& event!=     1923 
&& event!=     1931 
&& event!=     1949 
&& event!=     1961 
&& event!=     1969 
&& event!=     1975 
&& event!=     1979 
&& event!=     1985 
&& event!=     1986 
&& event!=     1992 
&& event!=       11 
&& event!=       33 
&& event!=       45 
&& event!=       58 
&& event!=       63 
&& event!=       77 
&& event!=      118 
&& event!=      127 
&& event!=      135 
&& event!=      149 
&& event!=      155 
&& event!=      178 
&& event!=      200 
&& event!=      209 
&& event!=      221 
&& event!=      265 
&& event!=      284 
&& event!=      339 
&& event!=      344 
&& event!=      349 
&& event!=      350 
&& event!=      357 
&& event!=      366 
&& event!=      382 
&& event!=      409 
&& event!=      414 
&& event!=      425 
&& event!=      434 
&& event!=      451 
&& event!=      473 
&& event!=      480 
&& event!=      488 
&& event!=      506 
&& event!=      511 
&& event!=      518 
&& event!=      525 
&& event!=      531 
&& event!=      543 
&& event!=      577 
&& event!=      698 
&& event!=      709 
&& event!=      764 
&& event!=      771 
&& event!=      779 
&& event!=      810 
&& event!=      847 
&& event!=      891 
&& event!=      901 
&& event!=      928 
&& event!=      972 
&& event!=      983 
&& event!=      990 
&& event!=      998 
&& event!=     1036 
&& event!=     1041 
&& event!=     1108 
&& event!=     1114 
&& event!=     1140 
&& event!=     1167 
&& event!=     1171 
&& event!=     1173 
&& event!=     1210 
&& event!=     1211 
&& event!=     1217 
&& event!=     1219 
&& event!=     1252 
&& event!=     1255 
&& event!=     1275 
&& event!=     1294 
&& event!=     1316 
&& event!=     1341 
&& event!=     1349 
&& event!=     1351 
&& event!=     1357 
&& event!=     1361 
&& event!=     1366 
&& event!=     1388 
&& event!=     1417 
&& event!=     1441 
&& event!=     1462 
&& event!=     1466 
&& event!=     1477 
&& event!=     2046 
&& event!=     2060 
&& event!=     2067 
&& event!=     2127 
&& event!=     2160 
&& event!=     2181 
&& event!=     2198 
&& event!=     2214 
&& event!=     2231 
&& event!=     2241 
&& event!=     2244 
&& event!=     2269 
&& event!=     2291 
&& event!=     2340 
&& event!=     2365 
&& event!=     2375 
&& event!=     2377 
&& event!=     2381 
&& event!=     2387 
&& event!=     2397 
&& event!=     2431 
&& event!=     2433 
&& event!=     2437 
&& event!=    36512 
&& event!=    36540 
&& event!=    36583 
&& event!=    36585 
&& event!=    36631 
&& event!=    36643 
&& event!=    36666 
&& event!=    36694 
&& event!=    36698 
&& event!=    36724 
&& event!=    36734 
&& event!=    36759 
&& event!=    36771 
&& event!=    36776 
&& event!=    36786 
&& event!=    36790 
&& event!=    36798 
&& event!=    36799 
&& event!=    36805 
&& event!=    36816 
&& event!=    36829 
&& event!=    36843 
&& event!=    36885 
&& event!=    36890 
&& event!=    36898 
&& event!=    36905 
&& event!=    36917 
&& event!=    36969 
&& event!=    36974 
&& event!=    36984 
&& event!=    37040 
&& event!=    37041 
&& event!=    37115 
&& event!=    37141 
&& event!=    37162 
&& event!=    37164 
&& event!=    37165 
&& event!=    37171 
&& event!=    37189 
&& event!=    37203 
&& event!=    37212 
&& event!=    37238 
&& event!=    37294 
&& event!=    37302 
&& event!=    37327 
&& event!=    37330 
&& event!=    37339 
&& event!=    37356 
&& event!=    37372 
&& event!=    37388 
&& event!=    37423 
&& event!=    37454 
&& event!=    37455 
&& event!=    37458 
&& event!=    37461 
&& event!=    37464 
&& event!=    37466 
&& event!=    37480 
&& event!=    37496 
&& event!=    37515 
&& event!=    37564 
&& event!=    37572 
&& event!=    37573 
&& event!=    37605 
&& event!=    37607 
&& event!=    37623 
&& event!=    37644 
&& event!=    37656 
&& event!=    37697 
&& event!=    37707 
&& event!=    37713 
&& event!=    37752 
&& event!=    37758 
&& event!=    37765 
&& event!=    37808 
&& event!=    37835 
&& event!=    37864 
&& event!=    37869 
&& event!=    37884 
&& event!=    37905 
&& event!=    37923 
&& event!=    37937 
&& event!=    37990 
&& event!=    38001 
&& event!=    38015 
&& event!=    38027 
&& event!=    38068 
&& event!=    38071 
&& event!=    38076 
&& event!=    38082 
&& event!=    38086 
&& event!=    38110 
&& event!=    38142 
&& event!=    38146 
&& event!=    38194 
&& event!=    38201 
&& event!=    38225 
&& event!=    38235 
&& event!=    38237 
&& event!=    38260 
&& event!=    38262 
&& event!=    38278 
&& event!=    38287 
&& event!=    38299 
&& event!=    38303 
&& event!=    38340 
&& event!=    38346 
&& event!=    38352 
&& event!=    38353 
&& event!=    38357 
&& event!=    38384 
&& event!=    38400 
&& event!=    38417 
&& event!=    38448 
&& event!=    38461 
&& event!=    38487 
&& event!=    38492 
&& event!=    38498 
&& event!=    38519 
&& event!=    38524 
&& event!=    38608 
&& event!=    38621 
&& event!=    38634 
&& event!=    38660 
&& event!=    38732 
&& event!=    38739 
&& event!=    38754 
&& event!=    38760 
&& event!=    38769 
&& event!=    38787 
&& event!=    38804 
&& event!=    38810 
&& event!=    38859 
&& event!=    38862 
&& event!=    38873 
&& event!=    38889 
&& event!=    38918 
&& event!=    38922 
&& event!=    38932 
&& event!=    38941 
&& event!=    38946 
&& event!=    38951 
&& event!=    39001 
&& event!=    39083 
&& event!=    39096 
&& event!=    39097 
&& event!=    39168 
&& event!=    39172 
&& event!=    39190 
&& event!=    39191 
&& event!=    39196 
&& event!=    39240 
&& event!=    39257 
&& event!=    39300 
&& event!=    39311 
&& event!=    39345 
&& event!=    39385 
&& event!=    39402 
&& event!=    39420 
&& event!=    39449 
&& event!=    39464 
&& event!=    39476 
&& event!=    39498 
&& event!=    39500 
&& event!=    39507 
&& event!=    39539 
&& event!=    39548 
&& event!=    39578 
&& event!=    39583 
&& event!=    39656 
&& event!=    39665 
&& event!=    39695 
&& event!=    39727 
&& event!=    39729 
&& event!=    39756 
&& event!=    39811 
&& event!=    39913 
&& event!=    39942 
&& event!=    39973 
&& event!=    39988 
&& event!=    39992 
&& event!=    40007 
&& event!=    40011 
&& event!=    40064 
&& event!=    40160 
&& event!=    40174 
&& event!=    40196 
&& event!=    40228 
&& event!=    40247 
&& event!=    40301 
&& event!=    40319 
&& event!=    40325 
&& event!=    40341 
&& event!=    40357 
&& event!=    40359 
&& event!=    40373 
&& event!=    40393 
&& event!=    40401 
&& event!=    40410 
&& event!=    40442 
&& event!=    40459 
&& event!=    40468 
&& event!=    40517 
&& event!=    40531 
&& event!=    40535 
&& event!=    40538 
&& event!=    40556 
&& event!=    40618 
&& event!=    40645 
&& event!=    40700 
&& event!=    40772 
&& event!=    40775 
&& event!=    40777 
&& event!=    40835 
&& event!=    40861 
&& event!=    40867 
&& event!=    40878 
&& event!=    40890 
&& event!=    40905 
&& event!=    40929 
&& event!=    40933 
&& event!=    40973 
&& event!=    40974 
&& event!=    41066 
&& event!=    41108 
&& event!=    41115 
&& event!=    41139 
&& event!=    41214 
&& event!=    41215 
&& event!=    41218 
&& event!=    41225 
&& event!=    41229 
&& event!=    41255 
&& event!=    41283 
&& event!=    41285 
&& event!=    41318 
&& event!=    41320 
&& event!=    41332 
&& event!=    41335 
&& event!=    41352 
&& event!=    41444 
&& event!=    41455 
&& event!=    41480 
&& event!=    41509 
&& event!=    41528 
&& event!=    41540 
&& event!=    41545 
&& event!=    41563 
&& event!=    41574 
&& event!=    41584 
&& event!=    41632 
&& event!=    41685 
&& event!=    41698 
&& event!=    41704 
&& event!=    41730 
&& event!=    41762 
&& event!=    41763 
&& event!=    41767 
&& event!=    41778 
&& event!=    41816 
&& event!=    41823 
&& event!=    41837 
&& event!=    41843 
&& event!=    41885 
&& event!=    41913 
&& event!=    41941 
&& event!=    41944 
&& event!=    41950 
&& event!=    41963 
&& event!=    41967 
&& event!=    42013 
&& event!=    42019 
&& event!=    42043 
&& event!=    42044 
&& event!=    42055 
&& event!=    42080 
&& event!=    42103 
&& event!=    42106 
&& event!=    42109 
&& event!=    42148 
&& event!=    42164 
&& event!=    42171 
&& event!=    42235 
&& event!=    42253 
&& event!=    42262 
&& event!=    42264 
&& event!=    42332 
&& event!=    42358 
&& event!=    42376 
&& event!=    42422 
&& event!=    42430 
&& event!=    42451 
&& event!=    42504 
&& event!=    42522 
&& event!=    42546 
&& event!=    42592 
&& event!=    42604 
&& event!=    42631 
&& event!=    42663 
&& event!=    42669 
&& event!=    42670 
&& event!=    42671 
&& event!=    42676 
&& event!=    42711 
&& event!=    42723 
&& event!=    42736 
&& event!=    42748 
&& event!=    42762 
&& event!=    42790 
&& event!=    42801 
&& event!=    42812 
&& event!=    42826 
&& event!=    42840 
&& event!=    42873 
&& event!=    42879 
&& event!=    42882 
&& event!=    42884 
&& event!=    42888 
&& event!=    42895 
&& event!=    42920 
&& event!=    42949 
&& event!=    42952 
&& event!=    42999 
&& event!=    43021 
&& event!=    43044 
&& event!=    43067 
&& event!=    43097 
&& event!=    43112 
&& event!=    43120 
&& event!=    43135 
&& event!=    43195 
&& event!=    43204 
&& event!=    43249 
&& event!=    43318 
&& event!=    43330 
&& event!=    43339 
&& event!=    43367 
&& event!=    43390 
&& event!=    43391 
&& event!=    43404 
&& event!=    43412 
&& event!=    43439 
&& event!=    43441 
&& event!=    43453 
&& event!=    43505 
&& event!=    43522 
&& event!=    43565 
&& event!=    43625 
&& event!=    43691 
&& event!=    43741 
&& event!=    43832 
&& event!=    43835 
&& event!=    43879 
&& event!=    43885 
&& event!=    44005 
&& event!=    44023 
&& event!=    44041 
&& event!=    44126 
&& event!=    44128 
&& event!=    44144 
&& event!=    44155 
&& event!=    44221 
&& event!=    44222 
&& event!=    44231 
&& event!=    44250 
&& event!=    44254 
&& event!=    44262 
&& event!=    44340 
&& event!=    44353 
&& event!=    44355 
&& event!=    44367 
&& event!=    44401 
&& event!=    44406 
&& event!=    44423 
&& event!=    44425 
&& event!=    44435 
&& event!=    44444 
&& event!=    44482 
&& event!=    44494 
&& event!=    44532 
&& event!=    44544 
&& event!=    44546 
&& event!=    44547 
&& event!=    44551 
&& event!=    44557 
&& event!=    44571 
&& event!=    44591 
&& event!=    44601 
&& event!=    44646 
&& event!=    44647 
&& event!=    44693 
&& event!=    44697 
&& event!=    44709 
&& event!=    44751 
&& event!=    44777 
&& event!=    44804 
&& event!=    44840 
&& event!=    44844 
&& event!=    44849 
&& event!=    44932 
&& event!=    44975 
&& event!=    44978 
&& event!=    45518 
&& event!=    45519 
&& event!=    45528 
&& event!=    45556 
&& event!=    45570 
&& event!=    45574 
&& event!=    45601 
&& event!=    45665 
&& event!=    45666 
&& event!=    45668 
&& event!=    45683 
&& event!=    45687 
&& event!=    45701 
&& event!=    45709 
&& event!=    45732 
&& event!=    45763 
&& event!=    45770 
&& event!=    45781 
&& event!=    45782 
&& event!=    45790 
&& event!=    45797 
&& event!=    45815 
&& event!=    45822 
&& event!=    45829 
&& event!=    45879 
&& event!=    45885 
&& event!=    45900 
&& event!=    45920 
&& event!=    45944 
&& event!=    45953 
&& event!=    45975 
&& event!=    45980 
&& event!=    45985 
&& event!=    45999 
&& event!=    46006 
&& event!=    46064 
&& event!=    46081 
&& event!=    46082 
&& event!=    46086 
&& event!=    46094 
&& event!=    46096 
&& event!=    46103 
&& event!=    46142 
&& event!=    46153 
&& event!=    46161 
&& event!=    46173 
&& event!=    46186 
&& event!=    46190 
&& event!=    46196 
&& event!=    46230 
&& event!=    46232 
&& event!=    46233 
&& event!=    46285 
&& event!=    46308 
&& event!=    46326 
&& event!=    46342 
&& event!=    46352 
&& event!=    46354 
&& event!=    46359 
&& event!=    46367 
&& event!=    46372 
&& event!=    46374 
&& event!=    46406 
&& event!=    46442 
&& event!=    46453 
&& event!=    46517 
&& event!=    46518 
&& event!=    46526 
&& event!=    46570 
&& event!=    46572 
&& event!=    46577 
&& event!=    46581 
&& event!=    46582 
&& event!=    46608 
&& event!=    46615 
&& event!=    46634 
&& event!=    46653 
&& event!=    46655 
&& event!=    46658 
&& event!=    46662 
&& event!=    46719 
&& event!=    46761 
&& event!=    46818 
&& event!=    46830 
&& event!=    46863 
&& event!=    46879 
&& event!=    46891 
&& event!=    46917 
&& event!=    46928 
&& event!=    46993 
&& event!=    46996 
&& event!=    47026 
&& event!=    47031 
&& event!=    47040 
&& event!=    47043 
&& event!=    47052 
&& event!=    47077 
&& event!=    47091 
&& event!=    47097 
&& event!=    47105 
&& event!=    47151 
&& event!=    47156 
&& event!=    47230 
&& event!=    47278 
&& event!=    47283 
&& event!=    47291 
&& event!=    47324 
&& event!=    47348 
&& event!=    47349 
&& event!=    47361 
&& event!=    47381 
&& event!=    47389 
&& event!=    47403 
&& event!=    47411 
&& event!=    47451 
&& event!=    47457 
&& event!=    47476 
&& event!=    47518 
&& event!=    47548 
&& event!=    47553 
&& event!=    47565 
&& event!=    47578 
&& event!=    47591 
&& event!=    47593 
&& event!=    47598 
&& event!=    47602 
&& event!=    47611 
&& event!=    47633 
&& event!=    47646 
&& event!=    47649 
&& event!=    47693 
&& event!=    47711 
&& event!=    47738 
&& event!=    47741 
&& event!=    47752 
&& event!=    47767 
&& event!=    47770 
&& event!=    47780 
&& event!=    47782 
&& event!=    47792 
&& event!=    47797 
&& event!=    47824 
&& event!=    47826 
&& event!=    47840 
&& event!=    47862 
&& event!=    47877 
&& event!=    47899 
&& event!=    47917 
&& event!=    47990 
&& event!=    48005 
&& event!=    48008 
&& event!=    48009 
&& event!=    48013 
&& event!=    48099 
&& event!=    48104 
&& event!=    48137 
&& event!=    48165 
&& event!=    48207 
&& event!=    48209 
&& event!=    48217 
&& event!=    48279 
&& event!=    48290 
&& event!=    48293 
&& event!=    48306 
&& event!=    48348 
&& event!=    48458 
&& event!=    48465 
&& event!=    48473 
&& event!=    48474 
&& event!=    48475 
&& event!=    48487 
&& event!=    48508 
&& event!=    48525 
&& event!=    48527 
&& event!=    48575 
&& event!=    48582 
&& event!=    48601 
&& event!=    48701 
&& event!=    48702 
&& event!=    48704 
&& event!=    48728 
&& event!=    48782 
&& event!=    48789 
&& event!=    48797 
&& event!=    48811 
&& event!=    48825 
&& event!=    48875 
&& event!=    48907 
&& event!=    48969 
&& event!=    48976 
&& event!=    48981 
&& event!=    48985 
&& event!=    48987 
&& event!=    48992 
&& event!=    49003 
&& event!=    49022 
&& event!=    49064 
&& event!=    49073 
&& event!=    49156 
&& event!=    49223 
&& event!=    49226 
&& event!=    49227 
&& event!=    49247 
&& event!=    49259 
&& event!=    49263 
&& event!=    49297 
&& event!=    49321 
&& event!=    49332 
&& event!=    49390 
&& event!=    49402 
&& event!=    49411 
&& event!=    49417 
&& event!=    49439 
&& event!=    49451 
&& event!=    49463 
&& event!=    49507 
&& event!=    49510 
&& event!=    49556 
&& event!=    49583 
&& event!=    49645 
&& event!=    49676 
&& event!=    49697 
&& event!=    49702 
&& event!=    49705 
&& event!=    49735 
&& event!=    49736 
&& event!=    49746 
&& event!=    49771 
&& event!=    49787 
&& event!=    49792 
&& event!=    49818 
&& event!=    49834 
&& event!=    49859 
&& event!=    49931 
&& event!=    49932 
&& event!=    49934 
&& event!=    49945 
&& event!=    49948 
&& event!=    49960 
&& event!=    49995 
&& event!=    49998 
&& event!=    11006 
&& event!=    11025 
&& event!=    11027 
&& event!=    11034 
&& event!=    11048 
&& event!=    11054 
&& event!=    11060 
&& event!=    11096 
&& event!=    11109 
&& event!=    11110 
&& event!=    11130 
&& event!=    11153 
&& event!=    11154 
&& event!=    11163 
&& event!=    11189 
&& event!=    11239 
&& event!=    11266 
&& event!=    11272 
&& event!=    11286 
&& event!=    11297 
&& event!=    11305 
&& event!=    11311 
&& event!=    11319 
&& event!=    11321 
&& event!=    11362 
&& event!=    11375 
&& event!=    11376 
&& event!=    11384 
&& event!=    11424 
&& event!=    11442 
&& event!=    11453 
&& event!=    11474 
&& event!=    11500 
&& event!=    11514 
&& event!=    11530 
&& event!=    11551 
&& event!=    11595 
&& event!=    11598 
&& event!=    11599 
&& event!=    11629 
&& event!=    11641 
&& event!=    11664 
&& event!=    11685 
&& event!=    11725 
&& event!=    11754 
&& event!=    11768 
&& event!=    11775 
&& event!=    11777 
&& event!=    11798 
&& event!=    11805 
&& event!=    11816 
&& event!=    11826 
&& event!=    11829 
&& event!=    11834 
&& event!=    11866 
&& event!=    11875 
&& event!=    11884 
&& event!=    11885 
&& event!=    11909 
&& event!=    11962 
&& event!=    11982 
&& event!=    12020 
&& event!=    12060 
&& event!=    12074 
&& event!=    12108 
&& event!=    12125 
&& event!=    12171 
&& event!=    12175 
&& event!=    12182 
&& event!=    12200 
&& event!=    12223 
&& event!=    12234 
&& event!=    12239 
&& event!=    12245 
&& event!=    12257 
&& event!=    12259 
&& event!=    12260 
&& event!=    12287 
&& event!=    12291 
&& event!=    12323 
&& event!=    12333 
&& event!=    12335 
&& event!=    12336 
&& event!=    12348 
&& event!=    12349 
&& event!=    12397 
&& event!=    12399 
&& event!=    12403 
&& event!=    12414 
&& event!=    12464 
&& event!=    12492 
&& event!=    12493 
&& event!=    12538 
&& event!=    12549 
&& event!=    12553 
&& event!=    12560 
&& event!=    12576 
&& event!=    12588 
&& event!=    12613 
&& event!=    12626 
&& event!=    12644 
&& event!=    12649 
&& event!=    12669 
&& event!=    12673 
&& event!=    12688 
&& event!=    12697 
&& event!=    12714 
&& event!=    12719 
&& event!=    12724 
&& event!=    12738 
&& event!=    12742 
&& event!=    12755 
&& event!=    12764 
&& event!=    12770 
&& event!=    12775 
&& event!=    12794 
&& event!=    12802 
&& event!=    12809 
&& event!=    12818 
&& event!=    12824 
&& event!=    12842 
&& event!=    12859 
&& event!=    12933 
&& event!=    12941 
&& event!=    12943 
&& event!=    12947 
&& event!=    12975 
&& event!=    13003 
&& event!=    13015 
&& event!=    13021 
&& event!=    13059 
&& event!=    13067 
&& event!=    13084 
&& event!=    13094 
&& event!=    13096 
&& event!=    13100 
&& event!=    13102 
&& event!=    13110 
&& event!=    13114 
&& event!=    13122 
&& event!=    13123 
&& event!=    13150 
&& event!=    13162 
&& event!=    13169 
&& event!=    13226 
&& event!=    13244 
&& event!=    13271 
&& event!=    13284 
&& event!=    13296 
&& event!=    13317 
&& event!=    13329 
&& event!=    13331 
&& event!=    13347 
&& event!=    13378 
&& event!=    13387 
&& event!=    13419 
&& event!=    13438 
&& event!=    13457 
&& event!=    13527 
&& event!=    13548 
&& event!=    13606 
&& event!=    13607 
&& event!=    13615 
&& event!=    13655 
&& event!=    13656 
&& event!=    13680 
&& event!=    13694 
&& event!=    13698 
&& event!=    13761 
&& event!=    13816 
&& event!=    13822 
&& event!=    13853 
&& event!=    13860 
&& event!=    13916 
&& event!=    13925 
&& event!=    13927 
&& event!=    13959 
&& event!=    13996 
&& event!=    17524 
&& event!=    17531 
&& event!=    17592 
&& event!=    17610 
&& event!=    17636 
&& event!=    17647 
&& event!=    17653 
&& event!=    17656 
&& event!=    17673 
&& event!=    17679 
&& event!=    17680 
&& event!=    17688 
&& event!=    17706 
&& event!=    17721 
&& event!=    17732 
&& event!=    17763 
&& event!=    17775 
&& event!=    17779 
&& event!=    17797 
&& event!=    17824 
&& event!=    17834 
&& event!=    17856 
&& event!=    17879 
&& event!=    17893 
&& event!=    17905 
&& event!=    17926 
&& event!=    17945 
&& event!=    17984 
&& event!=    14028 
&& event!=    14050 
&& event!=    14065 
&& event!=    14066 
&& event!=    14075 
&& event!=    14081 
&& event!=    14103 
&& event!=    14119 
&& event!=    14123 
&& event!=    14190 
&& event!=    14205 
&& event!=    14225 
&& event!=    14272 
&& event!=    14281 
&& event!=    14293 
&& event!=    14333 
&& event!=    14363 
&& event!=    14366 
&& event!=    14390 
&& event!=    14437 
&& event!=    14449 
&& event!=    14477 
&& event!=    16527 
&& event!=    16569 
&& event!=    16598 
&& event!=    16611 
&& event!=    16635 
&& event!=    16703 
&& event!=    16734 
&& event!=    16738 
&& event!=    16779 
&& event!=    16805 
&& event!=    16823 
&& event!=    16862 
&& event!=    16926 
&& event!=    16967 
&& event!=    17037 
&& event!=    17093 
&& event!=    17102 
&& event!=    17109 
&& event!=    17136 
&& event!=    17157 
&& event!=    17227 
&& event!=    17241 
&& event!=    17257 
&& event!=    17263 
&& event!=    17271 
&& event!=    17274 
&& event!=    17301 
&& event!=    17305 
&& event!=    17421 
&& event!=    17439 
&& event!=    17450 
&& event!=    17466 
&& event!=    17473 
&& event!=    17477 
&& event!=    17484 
&& event!=    17488 
&& event!=    17498 
&& event!=    18003 
&& event!=    18024 
&& event!=    18025 
&& event!=    18046 
&& event!=    18071 
&& event!=    18075 
&& event!=    18083 
&& event!=    18098 
&& event!=    18109 
&& event!=    18112 
&& event!=    18128 
&& event!=    18144 
&& event!=    18159 
&& event!=    18179 
&& event!=    18197 
&& event!=    18290 
&& event!=    18326 
&& event!=    18350 
&& event!=    18377 
&& event!=    18396 
&& event!=    18400 
&& event!=    18413 
&& event!=    18417 
&& event!=    18430 
&& event!=    18443 
&& event!=    18453 
&& event!=    18532 
&& event!=    18550 
&& event!=    18556 
&& event!=    18557 
&& event!=    18566 
&& event!=    18568 
&& event!=    18571 
&& event!=    18588 
&& event!=    18589 
&& event!=    18618 
&& event!=    18623 
&& event!=    18624 
&& event!=    18634 
&& event!=    18665 
&& event!=    18672 
&& event!=    18700 
&& event!=    18725 
&& event!=    18825 
&& event!=    18832 
&& event!=    18833 
&& event!=    18841 
&& event!=    18855 
&& event!=    18856 
&& event!=    18863 
&& event!=    18887 
&& event!=    18898 
&& event!=    18900 
&& event!=    18916 
&& event!=    18921 
&& event!=    18922 
&& event!=    18926 
&& event!=    18940 
&& event!=    18965 
&& event!=    19004 
&& event!=    19016 
&& event!=    19022 
&& event!=    19028 
&& event!=    19032 
&& event!=    19044 
&& event!=    19119 
&& event!=    19128 
&& event!=    19136 
&& event!=    19152 
&& event!=    19201 
&& event!=    19228 
&& event!=    19232 
&& event!=    19256 
&& event!=    19265 
&& event!=    19285 
&& event!=    19292 
&& event!=    19309 
&& event!=    19312 
&& event!=    19362 
&& event!=    19410 
&& event!=    19434 
&& event!=    19443 
&& event!=    19456 
&& event!=    19469 
&& event!=    19480 
&& event!=    19493 
&& event!=    19508 
&& event!=    19516 
&& event!=    19521 
&& event!=    19526 
&& event!=    19534 
&& event!=    19558 
&& event!=    19607 
&& event!=    19613 
&& event!=    19618 
&& event!=    19621 
&& event!=    19673 
&& event!=    19682 
&& event!=    19705 
&& event!=    19712 
&& event!=    19714 
&& event!=    19730 
&& event!=    19773 
&& event!=    19780 
&& event!=    19781 
&& event!=    19818 
&& event!=    19858 
&& event!=    19869 
&& event!=    19885 
&& event!=    19916 
&& event!=    19982 
&& event!=    19994 
&& event!=    20009 
&& event!=    20040 
&& event!=    20054 
&& event!=    20057 
&& event!=    20079 
&& event!=    20086 
&& event!=    20088 
&& event!=    20091 
&& event!=    20097 
&& event!=    20108 
&& event!=    20144 
&& event!=    20147 
&& event!=    20148 
&& event!=    20162 
&& event!=    20164 
&& event!=    20177 
&& event!=    20192 
&& event!=    20213 
&& event!=    20222 
&& event!=    20231 
&& event!=    20244 
&& event!=    20262 
&& event!=    20273 
&& event!=    20294 
&& event!=    20323 
&& event!=    20370 
&& event!=    20395 
&& event!=    20450 
&& event!=    20476 
&& event!=    20519 
&& event!=    20533 
&& event!=    20585 
&& event!=    20591 
&& event!=    20604 
&& event!=    20608 
&& event!=    20629 
&& event!=    20640 
&& event!=    20644 
&& event!=    20665 
&& event!=    20674 
&& event!=    20684 
&& event!=    20741 
&& event!=    20750 
&& event!=    20752 
&& event!=    20809 
&& event!=    20820 
&& event!=    20852 
&& event!=    20854 
&& event!=    20866 
&& event!=    20872 
&& event!=    20891 
&& event!=    20910 
&& event!=    20941 
&& event!=    20942 
&& event!=    20946 
&& event!=    20977 
&& event!=    20978 
&& event!=    20982 
&& event!=    21027 
&& event!=    21038 
&& event!=    21069 
&& event!=    21077 
&& event!=    21098 
&& event!=    21104 
&& event!=    21117 
&& event!=    21175 
&& event!=    21177 
&& event!=    21234 
&& event!=    21259 
&& event!=    21304 
&& event!=    21326 
&& event!=    21353 
&& event!=    21361 
&& event!=    21419 
&& event!=    21424 
&& event!=    21494 
&& event!=    21497 
&& event!=    21532 
&& event!=    21608 
&& event!=    21612 
&& event!=    21660 
&& event!=    21678 
&& event!=    21682 
&& event!=    21748 
&& event!=    21755 
&& event!=    21778 
&& event!=    21786 
&& event!=    21797 
&& event!=    21812 
&& event!=    21814 
&& event!=    21850 
&& event!=    21853 
&& event!=    21928 
&& event!=    21929 
&& event!=    21984 
&& event!=    22051 
&& event!=    22064 
&& event!=    22073 
&& event!=    22076 
&& event!=    22078 
&& event!=    22086 
&& event!=    22106 
&& event!=    22138 
&& event!=    22151 
&& event!=    22158 
&& event!=    22170 
&& event!=    22187 
&& event!=    22194 
&& event!=    22277 
&& event!=    22284 
&& event!=    22310 
&& event!=    22321 
&& event!=    22338 
&& event!=    22346 
&& event!=    22360 
&& event!=    22367 
&& event!=    22371 
&& event!=    22414 
&& event!=    22457 
&& event!=    22495 
&& event!=    22511 
&& event!=    22556 
&& event!=    22629 
&& event!=    22715 
&& event!=    22719 
&& event!=    22720 
&& event!=    22757 
&& event!=    22786 
&& event!=    22826 
&& event!=    22862 
&& event!=    22867 
&& event!=    22897 
&& event!=    22912 
&& event!=    22925 
&& event!=    22994 
&& event!=    23009 
&& event!=    23061 
&& event!=    23066 
&& event!=    23075 
&& event!=    23085 
&& event!=    23099 
&& event!=    23154 
&& event!=    23199 
&& event!=    23223 
&& event!=    23229 
&& event!=    23231 
&& event!=    23235 
&& event!=    23245 
&& event!=    23268 
&& event!=    23288 
&& event!=    23300 
&& event!=    23318 
&& event!=    23324 
&& event!=    23336 
&& event!=    23337 
&& event!=    23378 
&& event!=    23386 
&& event!=    23391 
&& event!=    23437 
&& event!=    23451 
&& event!=    23473 
&& event!=    23504 
&& event!=    23530 
&& event!=    23563 
&& event!=    23594 
&& event!=    23609 
&& event!=    23656 
&& event!=    23668 
&& event!=    23679 
&& event!=    23700 
&& event!=    23706 
&& event!=    23708 
&& event!=    23721 
&& event!=    23727 
&& event!=    23789 
&& event!=    23792 
&& event!=    23794 
&& event!=    23805 
&& event!=    23817 
&& event!=    23836 
&& event!=    23881 
&& event!=    23950 
&& event!=    23974 
&& event!=    23984 
&& event!=    23992 
&& event!=    24004 
&& event!=    24025 
&& event!=    24033 
&& event!=    24061 
&& event!=    24065 
&& event!=    24071 
&& event!=    24079 
&& event!=    24082 
&& event!=    24114 
&& event!=    24128 
&& event!=    24133 
&& event!=    24134 
&& event!=    24172 
&& event!=    24180 
&& event!=    24199 
&& event!=    24214 
&& event!=    24241 
&& event!=    24248 
&& event!=    24268 
&& event!=    24296 
&& event!=    24316 
&& event!=    24347 
&& event!=    24350 
&& event!=    24458 
&& event!=    24468 
&& event!=    24478 
&& event!=    24491 
&& event!=    24506 
&& event!=    24512 
&& event!=    24528 
&& event!=    24529 
&& event!=    24533 
&& event!=    24551 
&& event!=    24557 
&& event!=    24579 
&& event!=    24581 
&& event!=    24610 
&& event!=    24620 
&& event!=    24634 
&& event!=    24637 
&& event!=    24638 
&& event!=    24641 
&& event!=    24670 
&& event!=    24694 
&& event!=    24801 
&& event!=    24831 
&& event!=    24848 
&& event!=    24853 
&& event!=    24875 
&& event!=    24938 
&& event!=    24950 
&& event!=    24962 
&& event!=    24966 
&& event!=    24998 
&& event!=    25009 
&& event!=    25024 
&& event!=    25035 
&& event!=    25042 
&& event!=    25050 
&& event!=    25065 
&& event!=    25081 
&& event!=    25087 
&& event!=    25103 
&& event!=    25112 
&& event!=    25132 
&& event!=    25144 
&& event!=    25161 
&& event!=    25167 
&& event!=    25175 
&& event!=    25204 
&& event!=    25218 
&& event!=    25227 
&& event!=    25232 
&& event!=    25234 
&& event!=    25236 
&& event!=    25266 
&& event!=    25267 
&& event!=    25268 
&& event!=    25285 
&& event!=    25306 
&& event!=    25318 
&& event!=    25360 
&& event!=    25390 
&& event!=    25407 
&& event!=    25468 
&& event!=    25489 
&& event!=    25500 
&& event!=    25538 
&& event!=    25589 
&& event!=    25614 
&& event!=    25618 
&& event!=    25668 
&& event!=    25683 
&& event!=    25695 
&& event!=    25708 
&& event!=    25709 
&& event!=    25719 
&& event!=    25749 
&& event!=    25793 
&& event!=    25836 
&& event!=    25837 
&& event!=    25857 
&& event!=    25859 
&& event!=    25891 
&& event!=    25892 
&& event!=    25907 
&& event!=    25914 
&& event!=    25916 
&& event!=    25928 
&& event!=    25950 
&& event!=    25995 
&& event!=    27504 
&& event!=    27505 
&& event!=    27518 
&& event!=    27531 
&& event!=    27552 
&& event!=    27573 
&& event!=    27582 
&& event!=    27593 
&& event!=    27607 
&& event!=    27615 
&& event!=    27634 
&& event!=    27643 
&& event!=    27647 
&& event!=    27648 
&& event!=    27657 
&& event!=    27666 
&& event!=    27716 
&& event!=    27737 
&& event!=    27738 
&& event!=    27740 
&& event!=    27767 
&& event!=    27777 
&& event!=    27796 
&& event!=    27813 
&& event!=    27837 
&& event!=    27842 
&& event!=    27860 
&& event!=    27864 
&& event!=    27961 
&& event!=    27966 
&& event!=    27977 
&& event!=    27999 
&& event!=    26007 
&& event!=    26018 
&& event!=    26020 
&& event!=    26022 
&& event!=    26026 
&& event!=    26046 
&& event!=    26050 
&& event!=    26053 
&& event!=    26074 
&& event!=    26076 
&& event!=    26099 
&& event!=    26104 
&& event!=    26145 
&& event!=    26162 
&& event!=    26166 
&& event!=    26170 
&& event!=    26180 
&& event!=    26186 
&& event!=    26193 
&& event!=    26227 
&& event!=    26248 
&& event!=    26252 
&& event!=    26347 
&& event!=    26376 
&& event!=    26385 
&& event!=    26397 
&& event!=    26416 
&& event!=    26426 
&& event!=    26453 
&& event!=    26485 
&& event!=    26487 
&& event!=    26498 
&& event!=    26502 
&& event!=    26513 
&& event!=    26538 
&& event!=    26564 
&& event!=    26581 
&& event!=    26585 
&& event!=    26592 
&& event!=    26598 
&& event!=    26618 
&& event!=    26635 
&& event!=    26642 
&& event!=    26645 
&& event!=    26661 
&& event!=    26684 
&& event!=    26688 
&& event!=    26725 
&& event!=    26737 
&& event!=    26769 
&& event!=    26824 
&& event!=    26832 
&& event!=    26836 
&& event!=    26845 
&& event!=    26849 
&& event!=    26860 
&& event!=    26869 
&& event!=    26892 
&& event!=    26905 
&& event!=    26908 
&& event!=    26941 
&& event!=    26965 
&& event!=    27026 
&& event!=    27046 
&& event!=    27056 
&& event!=    27086 
&& event!=    27132 
&& event!=    27148 
&& event!=    27155 
&& event!=    27156 
&& event!=    27177 
&& event!=    27189 
&& event!=    27193 
&& event!=    27218 
&& event!=    27242 
&& event!=    27250 
&& event!=    27265 
&& event!=    27268 
&& event!=    27269 
&& event!=    27271 
&& event!=    27276 
&& event!=    27296 
&& event!=    27317 
&& event!=    27341 
&& event!=    27359 
&& event!=    27364 
&& event!=    27369 
&& event!=    27456 
&& event!=    27469 
&& event!=    27480 
&& event!=    27498 
&& event!=    27499 
&& event!=    28008 
&& event!=    28066 
&& event!=    28079 
&& event!=    28080 
&& event!=    28082 
&& event!=    28094 
&& event!=    28099 
&& event!=    28110 
&& event!=    28115 
&& event!=    28126 
&& event!=    28177 
&& event!=    28179 
&& event!=    28186 
&& event!=    28253 
&& event!=    28311 
&& event!=    28312 
&& event!=    28367 
&& event!=    28418 
&& event!=    28436 
&& event!=    28466 
&& event!=    28484 
&& event!=    28504 
&& event!=    28506 
&& event!=    28526 
&& event!=    28530 
&& event!=    28553 
&& event!=    28608 
&& event!=    28624 
&& event!=    28625 
&& event!=    28667 
&& event!=    28681 
&& event!=    28701 
&& event!=    28717 
&& event!=    28723 
&& event!=    28738 
&& event!=    28754 
&& event!=    28822 
&& event!=    28906 
&& event!=    28912 
&& event!=    28945 
&& event!=    28988 
&& event!=    29003 
&& event!=    29006 
&& event!=    29065 
&& event!=    29096 
&& event!=    29140 
&& event!=    29154 
&& event!=    29209 
&& event!=    29213 
&& event!=    29219 
&& event!=    29228 
&& event!=    29248 
&& event!=    29256 
&& event!=    29264 
&& event!=    29269 
&& event!=    29274 
&& event!=    29276 
&& event!=    29329 
&& event!=    29330 
&& event!=    29341 
&& event!=    29346 
&& event!=    29359 
&& event!=    29373 
&& event!=    29394 
&& event!=    29439 
&& event!=    29442 
&& event!=    29468 
&& event!=    29472 
&& event!=    29493 
&& event!=    29525 
&& event!=    29559 
&& event!=    29632 
&& event!=    29666 
&& event!=    29673 
&& event!=    29680 
&& event!=    29698 
&& event!=    29702 
&& event!=    29708 
&& event!=    29711 
&& event!=    29714 
&& event!=    29781 
&& event!=    29799 
&& event!=    29823 
&& event!=    29829 
&& event!=    29896 
&& event!=    29910 
&& event!=    29951 
&& event!=    29960 
&& event!=    29968 
&& event!=    30002 
&& event!=    30027 
&& event!=    30042 
&& event!=    30059 
&& event!=    30064 
&& event!=    30086 
&& event!=    30132 
&& event!=    30203 
&& event!=    30218 
&& event!=    30253 
&& event!=    30264 
&& event!=    30299 
&& event!=    30309 
&& event!=    30404 
&& event!=    30415 
&& event!=    30424 
&& event!=    30428 
&& event!=    30459 
&& event!=    30467 
&& event!=    30477 
&& event!=    30498 
&& event!=    30502 
&& event!=    30519 
&& event!=    30533 
&& event!=    30537 
&& event!=    30557 
&& event!=    30563 
&& event!=    30579 
&& event!=    30593 
&& event!=    30600 
&& event!=    30621 
&& event!=    30657 
&& event!=    30706 
&& event!=    30716 
&& event!=    30727 
&& event!=    30763 
&& event!=    30765 
&& event!=    30779 
&& event!=    30783 
&& event!=    30792 
&& event!=    30809 
&& event!=    30839 
&& event!=    30885 
&& event!=    30925 
&& event!=    30971 
&& event!=    31026 
&& event!=    31048 
&& event!=    31061 
&& event!=    31067 
&& event!=    31070 
&& event!=    31084 
&& event!=    31107 
&& event!=    31139 
&& event!=    31176 
&& event!=    31190 
&& event!=    31214 
&& event!=    31276 
&& event!=    31282 
&& event!=    31292 
&& event!=    31313 
&& event!=    31315 
&& event!=    31338 
&& event!=    31372 
&& event!=    31375 
&& event!=    31384 
&& event!=    31410 
&& event!=    31449 
&& event!=    31450 
&& event!=    31457 
&& event!=    31516 
&& event!=    31600 
&& event!=    31620 
&& event!=    31650 
&& event!=    31654 
&& event!=    31761 
&& event!=    31767 
&& event!=    31815 
&& event!=    31826 
&& event!=    31828 
&& event!=    31833 
&& event!=    31837 
&& event!=    31859 
&& event!=    31890 
&& event!=    31893 
&& event!=    31895 
&& event!=    31898 
&& event!=    31900 
&& event!=    31912 
&& event!=    31955 
&& event!=    31961 
&& event!=    31982 
&& event!=    31997 
&& event!=    32060 
&& event!=    32106 
&& event!=    32126 
&& event!=    32139 
&& event!=    32157 
&& event!=    32158 
&& event!=    32161 
&& event!=    32174 
&& event!=    32260 
&& event!=    32310 
&& event!=    32331 
&& event!=    32356 
&& event!=    32367 
&& event!=    32389 
&& event!=    32406 
&& event!=    32436 
&& event!=    32476 
&& event!=    32482 
&& event!=    32483 
&& event!=    32494 
&& event!=    32538 
&& event!=    32603 
&& event!=    32626 
&& event!=    32633 
&& event!=    32636 
&& event!=    32659 
&& event!=    32669 
&& event!=    32681 
&& event!=    32690 
&& event!=    32695 
&& event!=    32709 
&& event!=    32718 
&& event!=    32719 
&& event!=    32735 
&& event!=    32749 
&& event!=    32763 
&& event!=    32821 
&& event!=    32829 
&& event!=    32857 
&& event!=    32867 
&& event!=    32900 
&& event!=    32903 
&& event!=    32911 
&& event!=    32914 
&& event!=    32945 
&& event!=    32949 
&& event!=    33048 
&& event!=    33074 
&& event!=    33076 
&& event!=    33120 
&& event!=    33121 
&& event!=    33130 
&& event!=    33148 
&& event!=    33154 
&& event!=    33158 
&& event!=    33174 
&& event!=    33177 
&& event!=    33180 
&& event!=    33186 
&& event!=    33189 
&& event!=    33216 
&& event!=    33248 
&& event!=    33254 
&& event!=    33352 
&& event!=    33353 
&& event!=    33356 
&& event!=    33499 
&& event!=    34003 
&& event!=    34027 
&& event!=    34047 
&& event!=    34056 
&& event!=    34079 
&& event!=    34089 
&& event!=    34157 
&& event!=    34198 
&& event!=    34201 
&& event!=    34212 
&& event!=    34221 
&& event!=    34229 
&& event!=    34247 
&& event!=    34269 
&& event!=    34298 
&& event!=    34314 
&& event!=    34330 
&& event!=    34342 
&& event!=    34352 
&& event!=    34356 
&& event!=    34368 
&& event!=    34371 
&& event!=    34375 
&& event!=    34432 
&& event!=    34500 
&& event!=    34501 
&& event!=    34516 
&& event!=    34545 
&& event!=    34549 
&& event!=    34563 
&& event!=    34584 
&& event!=    34604 
&& event!=    34613 
&& event!=    34623 
&& event!=    34629 
&& event!=    34647 
&& event!=    34666 
&& event!=    34699 
&& event!=    34723 
&& event!=    34745 
&& event!=    34890 
&& event!=    34897 
&& event!=    34922 
&& event!=    34926 
&& event!=    34935 
&& event!=    34972 
&& event!=    35000 
&& event!=    35002 
&& event!=    35004 
&& event!=    35015 
&& event!=    35031 
&& event!=    35038 
&& event!=    35041 
&& event!=    35050 
&& event!=    35059 
&& event!=    35082 
&& event!=    35086 
&& event!=    35100 
&& event!=    35104 
&& event!=    35106 
&& event!=    35112 
&& event!=    35114 
&& event!=    35136 
&& event!=    35153 
&& event!=    35166 
&& event!=    35203 
&& event!=    35209 
&& event!=    35236 
&& event!=    35243 
&& event!=    35278 
&& event!=    35292 
&& event!=    35293 
&& event!=    35295 
&& event!=    35306 
&& event!=    35337 
&& event!=    35362 
&& event!=    35374 
&& event!=    35397 
&& event!=    35407 
&& event!=    35409 
&& event!=    35431 
&& event!=    35465 
&& event!=    35483 
&& event!=    35485 
&& event!=    35515 
&& event!=    35519 
&& event!=    35583 
&& event!=    35595 
&& event!=    35621 
&& event!=    35651 
&& event!=    35666 
&& event!=    35675 
&& event!=    35690 
&& event!=    35691 
&& event!=    35733 
&& event!=    35751 
&& event!=    35757 
&& event!=    35776 
&& event!=    35780 
&& event!=    35792 
&& event!=    35824 
&& event!=    35865 
&& event!=    35912 
&& event!=    35971 
&& event!=    35987 
&& event!=    36026 
&& event!=    36028 
&& event!=    36034 
&& event!=    36087 
&& event!=    36091 
&& event!=    36105 
&& event!=    36160 
&& event!=    36201 
&& event!=    36216 
&& event!=    36217 
&& event!=    36259 
&& event!=    36290 
&& event!=    36294 
&& event!=    36306 
&& event!=    36308 
&& event!=    36313 
&& event!=    36339 
&& event!=    36386 
&& event!=    36438 
&& event!=    36462 
&& event!=    36466 
&& event!=     5017 
&& event!=     5024 
&& event!=     5041 
&& event!=     5058 
&& event!=     5061 
&& event!=     5069 
&& event!=     5074 
&& event!=     5087 
&& event!=     5097 
&& event!=     5109 
&& event!=     5113 
&& event!=     5124 
&& event!=     5144 
&& event!=     5146 
&& event!=     5147 
&& event!=     5154 
&& event!=     5187 
&& event!=     5228 
&& event!=     5234 
&& event!=     5241 
&& event!=     5269 
&& event!=     5280 
&& event!=     5295 
&& event!=     5317 
&& event!=     5357 
&& event!=     5366 
&& event!=     5388 
&& event!=     5406 
&& event!=     5408 
&& event!=     5418 
&& event!=     5443 
&& event!=     5444 
&& event!=     5456 
&& event!=     5457 
&& event!=     5474 
&& event!=     2506 
&& event!=     2508 
&& event!=     2518 
&& event!=     2553 
&& event!=     2571 
&& event!=     2578 
&& event!=     2608 
&& event!=     2624 
&& event!=     2632 
&& event!=     2654 
&& event!=     2669 
&& event!=     2676 
&& event!=     2708 
&& event!=     2709 
&& event!=     2774 
&& event!=     2837 
&& event!=     2850 
&& event!=     2862 
&& event!=     2906 
&& event!=     2913 
&& event!=     2914 
&& event!=     2930 
&& event!=     2934 
&& event!=     2935 
&& event!=     2946 
&& event!=     2985 
&& event!=     2997 
&& event!=     3518 
&& event!=     3519 
&& event!=     3551 
&& event!=     3563 
&& event!=     3585 
&& event!=     3598 
&& event!=     3603 
&& event!=     3604 
&& event!=     3716 
&& event!=     3721 
&& event!=     3738 
&& event!=     3745 
&& event!=     3747 
&& event!=     3760 
&& event!=     3770 
&& event!=     3775 
&& event!=     3801 
&& event!=     3805 
&& event!=     3849 
&& event!=     3862 
&& event!=     3881 
&& event!=     3900 
&& event!=     3918 
&& event!=     3942 
&& event!=     3953 
&& event!=     3979 
&& event!=     3982 
&& event!=     4022 
&& event!=     4038 
&& event!=     4039 
&& event!=     4044 
&& event!=     4058 
&& event!=     4083 
&& event!=     4160 
&& event!=     4174 
&& event!=     4179 
&& event!=     4221 
&& event!=     4226 
&& event!=     4229 
&& event!=     4266 
&& event!=     4284 
&& event!=     4314 
&& event!=     4315 
&& event!=     4318 
&& event!=     4430 
&& event!=     4441 
&& event!=     4466 
&& event!=     4500 
&& event!=     5515 
&& event!=     5537 
&& event!=     5542 
&& event!=     5544 
&& event!=     5559 
&& event!=     5560 
&& event!=     5570 
&& event!=     5581 
&& event!=     5585 
&& event!=     5615 
&& event!=     5621 
&& event!=     5639 
&& event!=     5675 
&& event!=     5690 
&& event!=     5697 
&& event!=     5737 
&& event!=     5739 
&& event!=     5740 
&& event!=     5761 
&& event!=     5765 
&& event!=     5768 
&& event!=     5772 
&& event!=     5776 
&& event!=     5827 
&& event!=     5842 
&& event!=     5854 
&& event!=     5876 
&& event!=     5889 
&& event!=     5912 
&& event!=     5916 
&& event!=     5925 
&& event!=     5929 
&& event!=     5973 
&& event!=     5975 
&& event!=     6045 
&& event!=     6058 
&& event!=     6116 
&& event!=     6121 
&& event!=     6123 
&& event!=     6138 
&& event!=     6151 
&& event!=     6161 
&& event!=     6163 
&& event!=     6201 
&& event!=     6235 
&& event!=     6276 
&& event!=     6287 
&& event!=     6336 
&& event!=     6340 
&& event!=     6353 
&& event!=     6371 
&& event!=     6389 
&& event!=     6395 
&& event!=     6455 
&& event!=     6471 
&& event!=     6501 
&& event!=     6509 
&& event!=     6512 
&& event!=     6515 
&& event!=     6528 
&& event!=     6539 
&& event!=     6570 
&& event!=     6571 
&& event!=     6647 
&& event!=     6656 
&& event!=     6750 
&& event!=     6751 
&& event!=     6770 
&& event!=     6779 
&& event!=     6851 
&& event!=     6853 
&& event!=     6868 
&& event!=     6869 
&& event!=     6917 
&& event!=     6920 
&& event!=     6924 
&& event!=     6931 
&& event!=     6933 
&& event!=     6934 
&& event!=     6962 
&& event!=     6964 
&& event!=     6971 
&& event!=     6979 
&& event!=     7043 
&& event!=     7047 
&& event!=     7058 
&& event!=     7125 
&& event!=     7136 
&& event!=     7142 
&& event!=     7148 
&& event!=     7150 
&& event!=     7175 
&& event!=     7204 
&& event!=     7220 
&& event!=     7256 
&& event!=     7287 
&& event!=     7296 
&& event!=     7312 
&& event!=     7337 
&& event!=     7351 
&& event!=     7373 
&& event!=     7380 
&& event!=     7403 
&& event!=     7441 
&& event!=     7446 
&& event!=     7456 
&& event!=     7472 
&& event!=     7473 
&& event!=     7475 
&& event!=     7476 
&& event!=     7490 
&& event!=     7500 )
) ; //continue;

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


        isLep_ele = false;

      if( electrons.size() + muons.size() > 1 ) {

        continue;

      } else if( electrons.size() > 0 ) {

        //if( electrons[0].DeltaR( phot1 ) < 0.5 || electrons[0].DeltaR( phot2 ) < 0.5 ) continue;

        isLeptonic_t = true;
        isMu_t = false;
        isLep_ele = true;
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

