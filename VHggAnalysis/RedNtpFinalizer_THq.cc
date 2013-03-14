#include "RedNtpFinalizer_THq.h"

#include "QGLikelihood/interface/QGLikelihoodCalculator.h"

#include "BTagSFUtil/interface/BTagSFUtil.h"




#define DEBUG_EVENT_NUMBER_ -1




RedNtpFinalizer_THq::RedNtpFinalizer_THq( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "THq", dataset ) 
{

  bTaggerType_ = bTaggerType;

  setSelectionType(selectionType);

  BLIND_ = true;

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
   TH1D*  h1_nbjets_loose_2jets = new TH1D("nbjets_loose_2jets", "", 11, -0.5, 10.5);
   h1_nbjets_loose_2jets->Sumw2();
   TH1D*  h1_nbjets_medium_2jets = new TH1D("nbjets_medium_2jets", "", 11, -0.5, 10.5);
   h1_nbjets_medium_2jets->Sumw2();


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

   TH1D* h1_mgg= new TH1D("mgg", "", 80, 100., 180.);
   h1_mgg->Sumw2();




   // higgs:
   TH1D*  h1_ptDiphot = new TH1D("ptDiphot", "", 100, 0., 500.);
   h1_ptDiphot->Sumw2();
   TH1D*  h1_ptRunDiphot = new TH1D("ptRunDiphot", "", 100, 0., 500.);
   h1_ptRunDiphot->Sumw2();

   TH1D*  h1_deltaPhi_top_higgs = new TH1D("deltaPhi_top_higgs", "", 100, 0., 10.);
   h1_deltaPhi_top_higgs->Sumw2();


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




   bool isLeptonic;

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
   float ptRungg_t;

   float Ht_t;

   int hasPassedSinglePhot_t;
   int hasPassedDoublePhot_t;

   float m3_t;
   float minv_lnu_t;
   float minv_blnu_t;
   float minv_bgg_t;

   float eventWeight = 1.;

   TTree* tree_passedEvents = new TTree();
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
   tree_passedEvents->Branch( "ptRungg", &ptRungg_t, "ptRungg_t/F" );
   tree_passedEvents->Branch("epfMet", &epfMet, "epfMet");
   tree_passedEvents->Branch( "hasPassedSinglePhot", &hasPassedSinglePhot_t, "hasPassedSinglePhot_t/I" );
   tree_passedEvents->Branch( "hasPassedDoublePhot", &hasPassedDoublePhot_t, "hasPassedDoublePhot_t/I" );

 
//   std::string qglFileName = "/afs/cern.ch/work/p/pandolf/CMSSW_5_2_5/src/UserCode/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1.root";
   //std::string qglFileName="/afs/cern.ch/user/m/micheli/public/ttH/QG_QCD_Pt-15to3000_TuneZ2_Flat_8TeV_pythia6_Summer12-PU_S7_START52_V9-v1.root";
   //std::string qglFileName="/afs/cern.ch/work/p/pandolf/CMSSW_5_3_6/src/QG/QGLikelihood/test/Histos_2012.root";
   std::string qglFileName="/afs/cern.ch/work/p/pandolf/public/ReducedHisto_2012.root";

   QGLikelihoodCalculator *qglikeli = new QGLikelihoodCalculator( qglFileName );


   BTagSFUtil* btsfutil = new BTagSFUtil(bTaggerType_, 13);



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









   float N_all = 0.;
   float N_qMatched = 0.;
   float N_tMatched = 0.;
   float N_hMatched = 0.;
   float N_bMatched = 0.;
   float N_WqMatched = 0.;
   float N_WqbarMatched = 0.;






   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      if( ientry % 20000 == 0 ) std::cout << "-> Entry " << ientry << " / " << nentries << std::endl;

      nb = tree_->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      bool isMC = ( run<5 );

     


      if( isMC ) {

        if( nGenEvents_<=0 ) {
          std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
          exit(99);
        }

        eventWeight = xSection_ / nGenEvents_ ;
 
        h1_nvertex->Fill( nvtx, eventWeight );

        // pu/pt reeventWeighting
        if(dopureeventWeight_) eventWeight *= pu_weight; 
        if(doptreeventWeight_) eventWeight *= pt_weight; 

      }



      TLorentzVector diphot;
      diphot.SetPtEtaPhiM( ptggnewvtx, etagg, phigg, massggnewvtx);



      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << std::endl << std::endl << "#### DEBUG LOG FOR EVENT " << DEBUG_EVENT_NUMBER_ << std::endl << std::endl;
        std::cout << std::endl << "Here are the photons: " << std::endl;
        std::cout << "   Phot1:  pt: " << ptphot1 << " eta: " << etascphot1 << std::endl;
        std::cout << "   Phot2:  pt: " << ptphot2 << " eta: " << etascphot2 << std::endl;
        std::cout << std::endl << "Diphoton pt: " << diphot.Pt() << " mass: " << diphot.M() << std::endl;
      }



      // analysis cuts:

      //if(npu>=30) continue;

      // if(massggnewvtx<90 || massggnewvtx>190) continue;
      //      if(diphot.M()<100 || diphot.M()>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
         || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed photon eta acceptance." << std::endl;
      }
      //     if(ptphot1<ptphot1cut) continue; //pt first photon

      bool ebeb = TMath::Abs(etascphot1)<1.4442 && TMath::Abs(etascphot2)<1.4442;

      float triggerThreshPhot1=33.;
      float triggerThreshPhot2=25.;

      if(ptphot1 < triggerThreshPhot1) continue; //pt first photon
      if(ptphot2 < triggerThreshPhot2) continue; //pt second photon

      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed photon trigger pt thresholds." << std::endl;
        std::cout << std::endl << "Checking analysis thresholds: " << std::endl;
        std::cout << ptphot1 << " > " << ptphot1cut_ << "*" << massggnewvtx << "/120. = " << ptphot1cut_*massggnewvtx/120. << " ? " << std::endl;
        std::cout << ptphot2 << " > " << ptphot2cut_ << "*" << massggnewvtx << "/120. = " << ptphot2cut_*massggnewvtx/120. << " ? " << std::endl;

      }

      if(ptphot1 < ptphot1cut_*massggnewvtx/120.) continue; //pt first photon
      if(ptphot2 < ptphot2cut_* massggnewvtx/120.) continue; //pt second photon

      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed photon analysis pt thresholds." << std::endl;
      }

      //if(ptgg_thresh_>0 && diphot.Pt()< pthiggsmincut_) continue; //pt higgs min


//    if(ptjet1cut>0 && (ptcorrjet1<ptjet1cut_ || TMath::Abs(etajet1)>4.7)) continue; //pt first jet
//    if(ptjet2cut>0 && (ptcorrjet2<ptjet2cut_ || TMath::Abs(etajet2)>4.7)) continue; //pt second jet


       // photon identification

       // r9 categories
       bool isr9phot1(0), isr9phot2(0);
 
       if(TMath::Abs(etascphot1)<1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)<1.4442 && r9phot2>.94) isr9phot2 = 1;
       if(TMath::Abs(etascphot1)>1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)>1.4442 && r9phot2>.94) isr9phot2 = 1;
 
       //if(r9cat_ == 1) {
       //  if(!isr9phot1 || !isr9phot2) continue;
       //} else if (r9cat_ == 0){
       //  if(isr9phot1 && isr9phot2) continue;
       //} 
 
       // photon id
       bool idphot1(0), idphot2(0), looseidphot1(0), looseidphot2(0), pxlphot1(1), pxlphot2(1);
 
 //       if(pixelseedcut) { 
 // 	pxlphot1 = !pid_haspixelseedphot1;
 // 	pxlphot2 = !pid_haspixelseedphot2;
 //       }

      if(!invert_photonCuts_){       
       idphot1 = (idcicpfphot1 >= photonID_thresh_);
       idphot2 = (idcicpfphot2 >= photonID_thresh_);
      }else{
	idphot1 = (idcicpfphot1 >= photonID_thresh_);
	idphot2 = (idcicpfphot2 < photonID_thresh_);
      }

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

      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed photon ID." << std::endl;
      }

      if(diphot.M()<100 || diphot.M()>180) continue;

      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed mgg 100-180 cut." << std::endl;
      }

      
      //       if( !isMC && BLIND_ && massggnewvtx>120. && massggnewvtx<130. ) continue;
      h1_nvertex_PUW->Fill( nvtx, eventWeight );
      
      
      if( !( !isMC && BLIND_ && massggnewvtx>120. && massggnewvtx<130.) )       h1_mgg_prepresel->Fill( massggnewvtx, eventWeight );
      
      
      
      // jets
      int njets_selected = 0;
      int njets_selected_btagloose = 0;
      int njets_selected_btagmedium = 0;
      std::vector<int> index_MatchedJet;
      std::vector<int> index_selected;
      std::vector<int> index_selected_btagloose;
      std::vector<int> index_selected_btagmedium;
      Ht_t = 0.;
      
      for( unsigned ijet=0; ijet<njets; ++ijet ) {
      	if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
	  std::cout << ijet << "/" << njets << " pt: " << ptcorrjet[ijet] << " eta: " << etajet[ijet] << std::endl;
	}
	if( ptcorrjet[ijet] < ptjetthresh_count_ ) continue;
	if( fabs(etajet[ijet]) > etajetthresh_count_ ) continue;
	
        //jet PU ID:
        bool passedPUID = true;
        if(use_PUID_){
          if((ijet+1)>=njets_PUID_thresh_){
            if(TMath::Abs(etajet[ijet]) < 2.5) {
	      if(betastarjet[ijet]<0)passedPUID=false;
	      if(betastarjet[ijet] > 0.2 * log( nvtx - PUID_betastar_thresh_ ) ) passedPUID = false;
              if(rmsjet[ijet] > 0.06) passedPUID = false;
            } else if(TMath::Abs(etajet[ijet]) < 3.){
              if(rmsjet[ijet] > 0.05) passedPUID = false;
            } else {
              if(rmsjet[ijet] > 0.055) passedPUID = false;
            }
          }
          if( !passedPUID )continue;
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
        if( isMC )
          btsfutil->modifyBTagsWithSF_fast(btagged_loose, btagged_medium, ptcorrjet[ijet], etajet[ijet], partPdgIDjet[ijet], "mean" );
      
      
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
      
      
        Ht_t += ptcorrjet[ijet];
      
      
        //AnalysisJet thisJet;
        //thisJet.set
        //selectedJets.push_back( thisJet );
      
      } //for jets
      

      if( index_selected_btagmedium.size()==0 ) continue;

      int bJet_index = index_selected_btagmedium[0];
      TLorentzVector bJet;
      bJet.SetPtEtaPhiE(ptcorrjet[bJet_index],etajet[bJet_index],phijet[bJet_index],ecorrjet[bJet_index]);


      isLeptonic = (ptele1>ptlep_ttH_leptonic_thresh_ || ptmu1>ptlep_ttH_leptonic_thresh_);//cut on pt lepton
      int isMu=ptmu1>ptlep_ttH_leptonic_thresh_;

      int index_jetW1=-1;
      int index_jetW2=-1;
      

      TLorentzVector top;
      TLorentzVector lept, neutrino;


      if( isLeptonic ) { // *** LEPTONIC CHANNEL

        if( isMu )
          lept.SetPtEtaPhiE( ptmu1, etamu1, phimu1, enemu1 );
        else
          lept.SetPtEtaPhiE( ptele1, etaele1, phiele1, eneele1 );

        neutrino.SetPtEtaPhiE( epfMet,0,phipfMet,epfMet );

        top = bJet + lept + neutrino;



      } else { // *** HADRONIC CHANNEL


        float maxDeltaPhi=0;
        
        // define top candidate as the one maximising deltaPhi(top, higgs):
        for( unsigned i=0; i<index_selected.size(); ++i ) {
          for( unsigned j=i+1; j<index_selected.size(); ++j ) {

            if( i==index_selected_btagmedium[0] || j==index_selected_btagmedium[0] ) continue;


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
          

        if( index_jetW1<0 && index_jetW2<0 ) continue;

        TLorentzVector jetW1,jetW2;
        jetW1.SetPtEtaPhiE(ptcorrjet[index_jetW1],etajet[index_jetW1],phijet[index_jetW1],ecorrjet[index_jetW1]);
        jetW2.SetPtEtaPhiE(ptcorrjet[index_jetW2],etajet[index_jetW2],phijet[index_jetW2],ecorrjet[index_jetW2]);

        top = bJet + jetW1 + jetW2;

      }      
      
      

      // now look for q-jet:
      int index_qJet=-1;
      TLorentzVector qJet;
      for( unsigned i=0; i<index_selected.size(); ++i ) {

        if( i==index_selected_btagmedium[0] ) continue;
        if( i==index_jetW1 ) continue;
        if( i==index_jetW2 ) continue;

        if( ptcorrjet[i]<20. ) continue;
        if( fabs(etajet[i])<1. ) continue;

        index_qJet=i;
        qJet.SetPtEtaPhiE(ptcorrjet[i],etajet[i],phijet[i],ecorrjet[i]);

        break;
    
      } 
      
      if( index_qJet<0 ) continue;



      
      if( isSignalMC ) {

        TLorentzVector h, t, q, b, Wq, Wqbar;
        h.SetPtEtaPhiE( pt_h, eta_h, phi_h, e_h );
        t.SetPtEtaPhiE( pt_t, eta_t, phi_t, e_t );
        b.SetPtEtaPhiE( pt_b, eta_b, phi_b, e_b );
        q.SetPtEtaPhiE( pt_q, eta_q, phi_q, e_q );
        Wq.SetPtEtaPhiE( pt_Wq, eta_Wq, phi_Wq, e_Wq );
        Wqbar.SetPtEtaPhiE( pt_Wqbar, eta_Wqbar, phi_Wqbar, e_Wqbar );


        N_all+=eventWeight;
 
        float deltaR_q = q.DeltaR(qJet);
        float deltaR_t = t.DeltaR(top);
        float deltaPhi_t = fabs(t.DeltaPhi(top));
        float deltaR_b = b.DeltaR(bJet);
        float deltaR_h = h.DeltaR(diphot);

        if( deltaR_q<0.5 )
          N_qMatched+=eventWeight;
        if( deltaR_t<0.5 || (isLeptonic && deltaPhi_t<0.5) )
          N_tMatched+=eventWeight;
        if( deltaR_b<0.5 )
          N_bMatched+=eventWeight;
        if( deltaR_h<0.5 )
          N_hMatched+=eventWeight;
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
      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed njets cut." << std::endl;
      }

      if(njets_selected>njets_upper_thresh_)continue;
      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed njets upper cut." << std::endl;
      }

      if(njets_selected_btagloose<nbtagloose_thresh_) continue;
      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed nbjets loose cut." << std::endl;
      }

      if(njets_selected_btagmedium<nbtagmedium_thresh_) continue;
      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed nbjets medium cut." << std::endl;
      }

      
      if(Ht_t<Ht_thresh_)continue;
      if( event == DEBUG_EVENT_NUMBER_ || DEBUG_EVENT_NUMBER_==-999 ) {
        std::cout << "-> Passed HT  cut." << std::endl;
      }

      
      
      
      
      
      
      h1_ptphot0->Fill(ptphot1, eventWeight);
      h1_ptphot1->Fill(ptphot2, eventWeight);
      
      h1_ptrunphot0->Fill(ptphot1*120/massggnewvtx, eventWeight);
      h1_ptrunphot1->Fill(ptphot2*120/massggnewvtx, eventWeight);
      
      h1_etaphot0->Fill(etaphot1, eventWeight);
      h1_etaphot1->Fill(etaphot2, eventWeight);
      
      
      h1_njets->Fill( njets_selected, eventWeight );
      h1_nbjets_loose->Fill( njets_selected_btagloose, eventWeight );
      h1_nbjets_medium->Fill( njets_selected_btagmedium, eventWeight );
      
      if( njets_selected>=2 ) {
        h1_nbjets_loose_2jets->Fill( njets_selected_btagloose, eventWeight );
        h1_nbjets_medium_2jets->Fill( njets_selected_btagmedium, eventWeight );
      }
      
      
      
      if( !( !isMC && BLIND_ && massggnewvtx>120. && massggnewvtx<130.) )       h1_mgg_presel->Fill( massggnewvtx, eventWeight );
      h1_njets_presel->Fill( njets_selected, eventWeight ); 
      
      

      float qgl_qJet = qglikeli->computeQGLikelihood2012( qJet.Pt(), qJet.Eta(), rhoPF, nChg_QCjet[index_qJet]+nNeutral_ptCutjet[index_qJet], ptD_QCjet[index_qJet], axis2_QCjet[index_qJet] );
      
      h1_pt_qJet->Fill( qJet.Pt(), eventWeight ); 
      h1_eta_qJet->Fill( qJet.Eta(), eventWeight ); 
      h1_qgl_qJet->Fill( qgl_qJet, eventWeight ); 
      
      h1_pt_bJet->Fill( bJet.Pt(), eventWeight ); 
      h1_eta_bJet->Fill( bJet.Eta(), eventWeight ); 
      
      
      
      h1_Ht->Fill(Ht_t,eventWeight);
      
      
      
      h1_ptDiphot->Fill( diphot.Pt(), eventWeight );
      h1_ptRunDiphot->Fill( diphot.Pt()*120./massggnewvtx, eventWeight );
      
      
      
      
      if( !( !isMC && BLIND_ && massggnewvtx>120. && massggnewvtx<130.) )       h1_mgg->Fill( massggnewvtx, eventWeight );
      
      
      
      if( isLeptonic ) {

        //leptons
        h1_pt_lept->Fill(lept.Pt(),eventWeight);
        h1_eta_lept->Fill(lept.Eta(),eventWeight);
      
        //MET
        h1_pfMet->Fill(epfMet,eventWeight);
        float deltaPhi_met_bJet = fabs(bJet.DeltaPhi(neutrino));
        h1_deltaPhi_met_bJet->Fill(deltaPhi_met_bJet,eventWeight);

      }
      

      h1_deltaPhi_top_higgs->Fill(fabs(top.DeltaPhi(diphot)),eventWeight);
      
      
      // set tree vars:
      njets_t  = njets_selected;
      nbjets_loose_t = njets_selected_btagloose;
      nbjets_medium_t = njets_selected_btagmedium;
      ptPhot1_t = ptphot1;
      ptPhot2_t = ptphot2;
      ptRunPhot1_t = ptphot1*120./massggnewvtx;
      ptRunPhot2_t = ptphot2*120./massggnewvtx;
      etaPhot1_t = etaphot1;
      etaPhot2_t = etaphot2;
      if( !( !isMC && BLIND_ && massggnewvtx>120. && massggnewvtx<130.) ){
        mgg_t = diphot.M();
      }else{
        mgg_t = -1;
      }
      ptgg_t = diphot.Pt();
      ptRungg_t = diphot.Pt()*120./massggnewvtx;
      
      
      
      tree_passedEvents->Fill();


   } //for entries


   if( isSignalMC ) {
     std::cout << "-> Matched Higgs in " << N_hMatched << "/" << N_all << " (" << (float)N_hMatched*100./N_all << "%%) of the cases." << std::endl;
     std::cout << "-> Matched bJet  in " << N_bMatched << "/" << N_all << " (" << (float)N_bMatched*100./N_all << "%%) of the cases." << std::endl;
     std::cout << "-> Matched top   in " << N_tMatched << "/" << N_all << " (" << (float)N_tMatched*100./N_all << "%%) of the cases." << std::endl;
     std::cout << "-> Matched qJet  in " << N_qMatched << "/" << N_all << " (" << (float)N_qMatched*100./N_all << "%%) of the cases." << std::endl;
   }


   outFile_->cd();


   tree_passedEvents->Write();





   h1_nGenEvents->Write();

   h1_nvertex->Write();
   h1_nvertex_PUW->Write();

   h1_njets_presel->Write();
   h1_njets->Write();
   h1_nbjets_loose->Write();
   h1_nbjets_medium->Write();
   h1_nbjets_loose_2jets->Write();
   h1_nbjets_medium_2jets->Write();

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


   h1_ptDiphot->Write();
   h1_ptRunDiphot->Write();


   h1_mgg_prepresel->Write();
   h1_mgg_presel->Write();
   h1_mgg->Write();


   h1_pt_lept->Write();
   h1_eta_lept->Write();

   h1_deltaPhi_met_bJet->Write();
   h1_deltaPhi_top_higgs->Write();

   h1_deltaR_q->Write();
   h1_deltaR_t->Write();
   h1_deltaR_h->Write();
   h1_deltaR_b->Write();
   h1_deltaR_Wq->Write();
   h1_deltaR_Wqbar->Write();
   h1_deltaPhi_t->Write();


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
  
  dopureeventWeight_ = true;
  doptreeventWeight_ = true;
  r9cat_ = 1;
  photonID_thresh_ = 4;
  cs_ = false;


  ptphot1cut_ = 30.;
  ptphot2cut_ = 25.;

  ptjetthresh_count_ = 20.;
  etajetthresh_count_ = 2.4;


  njets_thresh_=0;
  njets_upper_thresh_=1000;
  nbtagloose_thresh_=0;
  nbtagmedium_thresh_=0;


  Ht_thresh_=0;
  invert_photonCuts_=false;
  use_PUID_=true;

  njets_PUID_thresh_=0;
  ptlep_ttH_leptonic_thresh_=20;

  PUID_betastar_thresh_=0.67;
  chooseBtaggedJets_=true;//be careful with this


  if( selectionType=="presel" ) {

    // leave all cuts to default


  } else {
    std::cout << std::endl << std::endl << "Selection '" << selectionType << "' currently not implemented. Exiting." << std::endl;
    exit(12345);

  }


}


