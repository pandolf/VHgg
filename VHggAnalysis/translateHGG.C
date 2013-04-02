#define translateHGG_cxx
#include "translateHGG.h"
#include <iostream>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

void translateHGG::Loop()
{

   if (fChain == 0) return;


   TFile* outfile = TFile::Open("HGG_translated.root", "RECREATE");
   outfile->cd();

   TTree* tree = new TTree("tree_passedEvents", "");

   tree->Branch( "run", &run );
   float eventWeight =1.;
   tree->Branch( "eventWeight", &eventWeight );
   bool isLeptonic =true;
   tree->Branch( "isLeptonic", &isLeptonic );
   tree->Branch( "nvertex", &nvtx );
   tree->Branch( "pt_phot1", &ph1_pt );
   tree->Branch( "eta_phot1", &ph1_eta );
   tree->Branch( "phi_phot1", &ph1_phi );
   tree->Branch( "e_phot1", &ph1_e );
   tree->Branch( "pt_phot2", &ph2_pt );
   tree->Branch( "eta_phot2", &ph2_eta );
   tree->Branch( "phi_phot2", &ph2_phi );
   tree->Branch( "e_phot2", &ph2_e );
   tree->Branch( "mgg", &PhotonsMass );
   tree->Branch( "ptgg", &dipho_pt );

   float BDT_lept_t;
   tree->Branch( "BDT_lept", &BDT_lept_t );

   int nJets;
   tree->Branch( "nJets", &nJets );
   int nCentralJets;
   tree->Branch( "nCentralJets", &nCentralJets );
   int nBJetsLoose;
   tree->Branch( "nBJetsLoose", &nBJetsLoose );

   float pt_bJet;
   tree->Branch( "pt_bJet", &pt_bJet );
   float eta_bJet;
   tree->Branch( "eta_bJet", &eta_bJet );
   float phi_bJet;
   tree->Branch( "phi_bJet", &phi_bJet );

   float pt_qJet;
   tree->Branch( "pt_qJet", &pt_qJet );
   float eta_qJet;
   tree->Branch( "eta_qJet", &eta_qJet );
   float phi_qJet;
   tree->Branch( "phi_qJet", &phi_qJet );

   float pt_lept;
   tree->Branch( "pt_lept", &pt_lept );
   float eta_lept;
   tree->Branch( "eta_lept", &eta_lept );
   float phi_lept;
   tree->Branch( "phi_lept", &phi_lept );
   bool isMu;
   tree->Branch( "isMu", &isMu );
   int charge_lept;
   tree->Branch( "charge_lept", &charge_lept );
   float deltaR_lept_phot;
   tree->Branch( "deltaR_lept_phot", &deltaR_lept_phot );


   tree->Branch( "j1_betaStar", &j1_betaStar );
   tree->Branch( "j2_betaStar", &j2_betaStar );
   tree->Branch( "j3_betaStar", &j3_betaStar );
   tree->Branch( "j4_betaStar", &j4_betaStar );

   tree->Branch( "j1_eta", &j1_eta );
   tree->Branch( "j2_eta", &j2_eta );
   tree->Branch( "j3_eta", &j3_eta );
   tree->Branch( "j4_eta", &j4_eta );

   float deltaPhi_top_higgs;
   tree->Branch( "deltaPhi_top_higgs", &deltaPhi_top_higgs );

   float mt_top;
   tree->Branch( "mt_top", &mt_top );


   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" ); 

   float njets_bdt;
   float nbjets_loose_bdt;
   float deltaPhi_top_higgs_bdt;
   float mt_top_bdt;
   float eta_qJet_bdt;
   float charge_lept_bdt;
   float pt_bJet_bdt;

   reader->AddVariable( "njets", &njets_bdt );
   reader->AddVariable( "nbjets_loose", &nbjets_loose_bdt );
   reader->AddVariable( "deltaPhi_top_higgs", &deltaPhi_top_higgs_bdt );
   reader->AddVariable( "mt_top", &mt_top_bdt );
   reader->AddVariable( "eta_qJet", &eta_qJet_bdt );
   reader->AddVariable( "charge_lept", &charge_lept_bdt );
   reader->AddVariable( "pt_bJet", &pt_bJet_bdt );


   TString methodName = "BDTG method";
   //TString weightfile = "TMVA/weights/prova_nonjets_BDTG.weights.xml";
   TString weightfile = "TMVA/weights/prova_BDTG.weights.xml";

   std::cout << "-> Booking BDT" << std::endl;
   reader->BookMVA( methodName, weightfile ); 




   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      if( jentry % 5000 ==0 )
        std::cout << "Entry: " << jentry << "/" << nentries << std::endl;


      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);  
      

      if( category!=11 ) continue;
      if( PhotonsMass>115. && PhotonsMass<135. ) continue; //blinding

      TLorentzVector phot1, phot2;
      phot1.SetPtEtaPhiE( ph1_pt, ph1_eta, ph1_phi, ph1_e );
      phot2.SetPtEtaPhiE( ph2_pt, ph2_eta, ph2_phi, ph2_e );
    
      TLorentzVector diphot;
      diphot.SetPtEtaPhiE( dipho_pt, dipho_eta, dipho_phi, dipho_E );


      TLorentzVector bJet;
      int nbj=0;
      int i_bJet=0;
      if( j1_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j1_pt, j1_eta, j1_phi, j1_e );
        i_bJet = 1;
      }
      if( j2_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j2_pt, j2_eta, j2_phi, j2_e );
        i_bJet = 2;
      }
      if( j3_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j3_pt, j3_eta, j3_phi, j3_e );
        i_bJet = 3;
      }
      if( j4_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j4_pt, j4_eta, j4_phi, j4_e );
        i_bJet = 4;
      }
      if( j5_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j5_pt, j5_eta, j5_phi, j5_e );
        i_bJet = 5;
      }
      if( j6_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j6_pt, j6_eta, j6_phi, j6_e );
        i_bJet = 6;
      }
      if( j7_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j7_pt, j7_eta, j7_phi, j7_e );
        i_bJet = 7;
      }
      if( j8_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j8_pt, j8_eta, j8_phi, j8_e );
        i_bJet = 8;
      }
      if( j9_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j9_pt, j9_eta, j9_phi, j9_e );
        i_bJet = 9;
      }
      if( j10_algoPF1_csvBtag>0.679 ) { 
        nbj++;
        bJet.SetPtEtaPhiE( j10_pt, j10_eta, j10_phi, j10_e );
        i_bJet = 10;
      }

      if( nbj!=1 ) continue;
      if( bJet.Pt()< 20. ) {
        std::cout << "AAAAAAAAH!!! LE FOTTUTE SCIMMIE." << std::endl;
        exit(11);
      }


      nBJetsLoose=0;
      if( j1_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j2_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j3_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j4_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j5_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j6_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j7_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j8_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j9_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }
      if( j10_algoPF1_csvBtag>0.244 ) { 
        nBJetsLoose++;
      }


      TLorentzVector qJet;
      int i_qJet=0;
      if( fabs(j1_eta)>1. && i_bJet!=1 ) {
        qJet.SetPtEtaPhiE( j1_pt, j1_eta, j1_phi, j1_e );
        i_qJet=1;
      } else if( fabs(j2_eta)>1. && i_bJet!=2 ) {
        qJet.SetPtEtaPhiE( j2_pt, j2_eta, j2_phi, j2_e );
        i_qJet=2;
      } else if( fabs(j3_eta)>1. && i_bJet!=3 ) {
        qJet.SetPtEtaPhiE( j3_pt, j3_eta, j3_phi, j3_e );
        i_qJet=3;
      } else if( fabs(j4_eta)>1. && i_bJet!=4 ) {
        qJet.SetPtEtaPhiE( j4_pt, j4_eta, j4_phi, j4_e );
        i_qJet=4;
      } else if( fabs(j5_eta)>1. && i_bJet!=5 ) {
        qJet.SetPtEtaPhiE( j5_pt, j5_eta, j5_phi, j5_e );
        i_qJet=5;
      } else if( fabs(j6_eta)>1. && i_bJet!=6 ) {
        qJet.SetPtEtaPhiE( j6_pt, j6_eta, j6_phi, j6_e );
        i_qJet=6;
      } else if( fabs(j7_eta)>1. && i_bJet!=7 ) {
        qJet.SetPtEtaPhiE( j7_pt, j7_eta, j7_phi, j7_e );
        i_qJet=7;
      } else if( fabs(j8_eta)>1. && i_bJet!=8 ) {
        qJet.SetPtEtaPhiE( j8_pt, j8_eta, j8_phi, j8_e );
        i_qJet=8;
      } else if( fabs(j9_eta)>1. && i_bJet!=9 ) {
        qJet.SetPtEtaPhiE( j9_pt, j9_eta, j9_phi, j9_e );
        i_qJet=9;
      } else if( fabs(j10_eta)>1. && i_bJet!=10 ) {
        qJet.SetPtEtaPhiE( j10_pt, j10_eta, j10_phi, j10_e );
        i_qJet=10;
      }


      if( qJet.Pt()< 20. ) {
        std::cout << "LE FOTTUTE SCIMMIE. FOTTUTE FOTTUTISSIME SCIMMIE." << std::endl;
        std::cout << "event: " << event << std::endl;
        std::cout <<  j1_pt << " " <<  j1_eta << " " <<  j1_phi << " " <<  j1_e << " " <<  j1_algoPF1_csvBtag << std::endl;
        std::cout <<  j2_pt << " " <<  j2_eta << " " <<  j2_phi << " " <<  j2_e << " " <<  j2_algoPF1_csvBtag << std::endl;
        std::cout <<  j3_pt << " " <<  j3_eta << " " <<  j3_phi << " " <<  j3_e << " " <<  j3_algoPF1_csvBtag << std::endl;
        std::cout <<  j4_pt << " " <<  j4_eta << " " <<  j4_phi << " " <<  j4_e << " " <<  j4_algoPF1_csvBtag << std::endl;
        continue;
      }


      nCentralJets=0;
      if( i_qJet!=1 && i_bJet!=1 && fabs(j1_eta)<1. && j1_pt>20. ) nCentralJets++;
      if( i_qJet!=2 && i_bJet!=2 && fabs(j2_eta)<1. && j2_pt>20. ) nCentralJets++;
      if( i_qJet!=3 && i_bJet!=3 && fabs(j3_eta)<1. && j3_pt>20. ) nCentralJets++;
      if( i_qJet!=4 && i_bJet!=4 && fabs(j4_eta)<1. && j4_pt>20. ) nCentralJets++;
      if( i_qJet!=5 && i_bJet!=5 && fabs(j5_eta)<1. && j5_pt>20. ) nCentralJets++;
      if( i_qJet!=6 && i_bJet!=6 && fabs(j6_eta)<1. && j6_pt>20. ) nCentralJets++;
      if( i_qJet!=7 && i_bJet!=7 && fabs(j7_eta)<1. && j7_pt>20. ) nCentralJets++;
      if( i_qJet!=8 && i_bJet!=8 && fabs(j8_eta)<1. && j8_pt>20. ) nCentralJets++;
      if( i_qJet!=9 && i_bJet!=9 && fabs(j9_eta)<1. && j9_pt>20. ) nCentralJets++;
      if( i_qJet!=10 && i_bJet!=10 && fabs(j10_eta)<1. && j10_pt>20. ) nCentralJets++;
      if( nCentralJets>1 ) continue;
      //  std::cout << "SCIMMIE FOTTUTE." << std::endl;
      //  exit(11);
      //}

      nJets=0;
      if( j1_pt>20. ) nJets++;
      if( j2_pt>20. ) nJets++;
      if( j3_pt>20. ) nJets++;
      if( j4_pt>20. ) nJets++;
      if( j5_pt>20. ) nJets++;
      if( j6_pt>20. ) nJets++;
      if( j7_pt>20. ) nJets++;
      if( j8_pt>20. ) nJets++;
      if( j9_pt>20. ) nJets++;
      if( j10_pt>20. ) nJets++;
      if( nJets<3 ) continue;


      TLorentzVector lept;
      if( ptEle>10. && ptMu<=0. ) {
        isMu = false;
        charge_lept = chargeEle;
        lept.SetPtEtaPhiM( ptEle, etaEle, phiEle, 0. );
      } else if( ptEle<=10. && ptMu>0. ) {
        isMu = true;
        charge_lept = chargeMu;
        lept.SetPtEtaPhiM( ptMu, etaMu, phiMu, 0. );
      } else {
        std::cout << "ptEle: " <<  ptEle << " ptMu: " << ptMu << std::endl;
        continue;
      }


      pt_lept = lept.Pt();
      eta_lept = lept.Eta();
      phi_lept = lept.Phi();

      float deltaR_lept1 = lept.DeltaR(phot1);
      float deltaR_lept2 = lept.DeltaR(phot2);

      if( deltaR_lept1<0.1 ) {
        std::cout << "event: " << event << std::endl;
        std::cout << "lepton pt: " << lept.Pt() << " photon pt: " << phot1.Pt() << " deltaR: " << deltaR_lept1 << std::endl;
      }
      if( deltaR_lept2<0.1 ) {
        std::cout << "event: " << event << std::endl;
        std::cout << "lepton pt: " << lept.Pt() << " photon pt: " << phot2.Pt() << " deltaR: " << deltaR_lept2 << std::endl;
      }

      deltaR_lept_phot = TMath::Min(deltaR_lept1, deltaR_lept2);
      if( deltaR_lept_phot<0.5 ) continue;

      pt_bJet = bJet.Pt();
      eta_bJet = bJet.Eta();
      phi_bJet = bJet.Phi();

      pt_qJet = qJet.Pt();
      eta_qJet = qJet.Eta();
      phi_qJet = qJet.Phi();

      TLorentzVector neutrino;
      neutrino.SetPtEtaPhiE( met_pfmet, 0., met_phi_pfmet, met_pfmet );

      TLorentzVector top = neutrino + lept + bJet;

      

      njets_bdt              = nJets;
      nbjets_loose_bdt       = nBJetsLoose;
      deltaPhi_top_higgs_bdt = fabs(top.DeltaPhi(diphot));
      mt_top_bdt             = top.Mt();
      eta_qJet_bdt           = qJet.Eta();
      charge_lept_bdt        = charge_lept;
      pt_bJet_bdt            = bJet.Pt();
  
      BDT_lept_t = reader->EvaluateMVA( "BDTG method" );
      std::cout << "BDT: " << BDT_lept_t << std::endl;
      std::cout << "mgg: " << diphot.M() << std::endl;

      deltaPhi_top_higgs = fabs(deltaPhi_top_higgs_bdt);
      mt_top = top.Mt();

      tree->Fill();


   }


   std::cout << "Events passing selections: " << tree->GetEntries() << std::endl;

   outfile->cd();

   tree->Write();

   outfile->Close();

}
