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


   std::string outfilename = "HGG_translated_" + suffix_ + ".root";
   TFile* outfile = TFile::Open(outfilename.c_str(), "RECREATE");
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
   int nBJetsMedium;
   tree->Branch( "nBJetsMedium", &nBJetsMedium );

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
   float deltaEta_lept_qJet;
   tree->Branch( "deltaEta_lept_qJet", &deltaEta_lept_qJet );


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
   float deltaPhi_top_higgs_bdt;
   float mt_top_bdt;
   float eta_qJet_bdt;
   float charge_lept_bdt;
   float deltaEta_lept_qJet_bdt;
  

   reader->AddVariable( "njets", &njets_bdt );
   reader->AddVariable( "deltaPhi_top_higgs", &deltaPhi_top_higgs_bdt );
   reader->AddVariable( "mt_top", &mt_top_bdt );
   reader->AddVariable( "eta_qJet", &eta_qJet_bdt );
   reader->AddVariable( "charge_lept", &charge_lept_bdt );
   reader->AddVariable( "deltaEta_lept_qJet", &deltaEta_lept_qJet_bdt );


   TString methodName = "BDTG method";
   //TString weightfile = "TMVA/weights/prova_nonjets_BDTG.weights.xml";
   TString weightfile = "TMVA/weights/provaNEW_BDTG.weights.xml";

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
      if( PhotonsMass<100. || PhotonsMass>180. ) continue;
      if( PhotonsMass>115. && PhotonsMass<135. ) continue; //blinding

      TLorentzVector phot1, phot2;
      phot1.SetPtEtaPhiE( ph1_pt, ph1_eta, ph1_phi, ph1_e );
      phot2.SetPtEtaPhiE( ph2_pt, ph2_eta, ph2_phi, ph2_e );

      if( phot1.Pt() < 50.*PhotonsMass/120. ) continue;
      if( phot2.Pt() < 25. ) continue;

    
      TLorentzVector diphot;
      diphot.SetPtEtaPhiE( dipho_pt, dipho_eta, dipho_phi, dipho_E );


      TLorentzVector bJet;
      nBJetsMedium=0;
      int i_bJet=0;
      if( j1_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j1_pt, j1_eta, j1_phi, j1_e );
          i_bJet = 1;
        }
        nBJetsMedium++;
      }
      if( j2_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j2_pt, j2_eta, j2_phi, j2_e );
          i_bJet = 2;
        }
        nBJetsMedium++;
      }
      if( j3_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j3_pt, j3_eta, j3_phi, j3_e );
          i_bJet = 3;
        }
        nBJetsMedium++;
      }
      if( j4_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j4_pt, j4_eta, j4_phi, j4_e );
          i_bJet = 4;
        }
        nBJetsMedium++;
      }
      if( j5_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j5_pt, j5_eta, j5_phi, j5_e );
          i_bJet = 5;
        }
        nBJetsMedium++;
      }
      if( j6_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j6_pt, j6_eta, j6_phi, j6_e );
          i_bJet = 6;
        }
        nBJetsMedium++;
      }
      if( j7_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j7_pt, j7_eta, j7_phi, j7_e );
          i_bJet = 7;
        }
        nBJetsMedium++;
      }
      if( j8_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j8_pt, j8_eta, j8_phi, j8_e );
          i_bJet = 8;
        }
        nBJetsMedium++;
      }
      if( j9_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j9_pt, j9_eta, j9_phi, j9_e );
          i_bJet = 9;
        }
        nBJetsMedium++;
      }
      if( j10_algoPF1_csvBtag>0.679 ) { 
        if( nBJetsMedium==0 ) {
          bJet.SetPtEtaPhiE( j10_pt, j10_eta, j10_phi, j10_e );
          i_bJet = 10;
        }
        nBJetsMedium++;
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



      if( nBJetsMedium<1 ) continue;
      if( bJet.Pt()< 20. ) {
        std::cout << "AAAAAAAAH!!! LE FOTTUTE SCIMMIE." << std::endl;
        std::cout << "run: " << run << " event: " << event << std::endl;
        continue;
      }




      TLorentzVector lept;
      if( ptEle>10. && ptMu<=0. ) {
        isLeptonic = true;
        isMu = false;
        charge_lept = chargeEle;
        lept.SetPtEtaPhiM( ptEle, etaEle, phiEle, 0. );
      } else if( ptEle<=0. && ptMu>10. ) {
        isLeptonic = true;
        isMu = true;
        charge_lept = chargeMu;
        lept.SetPtEtaPhiM( ptMu, etaMu, phiMu, 0. );
      } else {
        isLeptonic = false;
      }


      std::vector<TLorentzVector> jets;
      TLorentzVector jet;
      if( j1_pt>20. ) {
        jet.SetPtEtaPhiE( j1_pt, j1_eta, j1_phi, j1_e );
        //jet.combinedSecondaryVertexBJetTag = j1_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j2_pt>20. ) {
        jet.SetPtEtaPhiE( j2_pt, j2_eta, j2_phi, j2_e );
        //jet.combinedSecondaryVertexBJetTag = j2_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j3_pt>20. ) {
        jet.SetPtEtaPhiE( j3_pt, j3_eta, j3_phi, j3_e );
        //jet.combinedSecondaryVertexBJetTag = j3_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j4_pt>20. ) {
        jet.SetPtEtaPhiE( j4_pt, j4_eta, j4_phi, j4_e );
        //jet.combinedSecondaryVertexBJetTag = j4_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j5_pt>20. ) {
        jet.SetPtEtaPhiE( j5_pt, j5_eta, j5_phi, j5_e );
        //jet.combinedSecondaryVertexBJetTag = j5_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j6_pt>20. ) {
        jet.SetPtEtaPhiE( j6_pt, j6_eta, j6_phi, j6_e );
        //jet.combinedSecondaryVertexBJetTag = j6_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j7_pt>20. ) {
        jet.SetPtEtaPhiE( j7_pt, j7_eta, j7_phi, j7_e );
        //jet.combinedSecondaryVertexBJetTag = j7_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j8_pt>20. ) {
        jet.SetPtEtaPhiE( j8_pt, j8_eta, j8_phi, j8_e );
        //jet.combinedSecondaryVertexBJetTag = j8_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j9_pt>20. ) {
        jet.SetPtEtaPhiE( j9_pt, j9_eta, j9_phi, j9_e );
        //jet.combinedSecondaryVertexBJetTag = j9_algoPF1_csvBtag;
        jets.push_back(jet);
      }
      if( j10_pt>20. ) {
        jet.SetPtEtaPhiE( j10_pt, j10_eta, j10_phi, j10_e );
        //jet.combinedSecondaryVertexBJetTag = j10_algoPF1_csvBtag;
        jets.push_back(jet);
      }


      nJets = jets.size();

      // now look for q-jet:
      int i_qJet=-1;
      TLorentzVector qJet;
      nCentralJets = 0;
      float hardestCentralJetPt=0.;
      float eta_thresh_qJet = (isLeptonic) ? 1. : 2.;
      for( unsigned ii=0; ii<jets.size(); ++ii ) {

        if( (ii+1)==i_bJet ) continue;

        if( jets[ii].Pt()<20. ) continue;

        if( isLeptonic ) {
          if( jets[ii].DeltaR(lept)<0.5 ) continue;
        }

        if( fabs(jets[ii].Eta())<eta_thresh_qJet ) {
          if( nCentralJets==0 ) hardestCentralJetPt = jets[ii].Pt(); //hardest central jet
          nCentralJets++;
        } else if( i_qJet<0 ) { //hardest forward jet
          i_qJet=ii+1;
          qJet = jets[ii];
        }

      }

      if( i_qJet<0 ) continue;


      if( isLeptonic ) {

std::cout << "**************************************************  isLeptonic" << std::endl;

        // *** LEPTONIC CHANNEL

        if( nCentralJets>1 ) continue;
        pt_lept = lept.Pt();
        eta_lept = lept.Eta();
        phi_lept = lept.Phi();

        if( fabs(qJet.Eta())<1. ) continue;
        if( qJet.Pt()<20. ) continue;


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

        deltaEta_lept_qJet = fabs( lept.Eta() - qJet.Eta() );

        njets_bdt              = nJets;
        deltaPhi_top_higgs_bdt = fabs(top.DeltaPhi(diphot));
        mt_top_bdt             = top.Mt();
        eta_qJet_bdt           = qJet.Eta();
        charge_lept_bdt        = charge_lept;
        deltaEta_lept_qJet_bdt = deltaEta_lept_qJet;
    
        BDT_lept_t = reader->EvaluateMVA( "BDTG method" );
        deltaPhi_top_higgs = fabs(deltaPhi_top_higgs_bdt);
        mt_top = top.Mt();

        std::cout << std::endl << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "*** Event passing LEPTONIC channel selection: run: " << run << " event: " << event << std::endl;
        std::cout << "+ phot1: pt: \t" << phot1.Pt() << " eta: \t" << phot1.Eta() << " phi: \t" << phot1.Phi() << std::endl;
        std::cout << "+ phot2: pt: \t" << phot2.Pt() << " eta: \t" << phot2.Eta() << " phi: \t" << phot2.Phi() << std::endl;
        std::cout << "+ Lepton: pt: \t" << lept.Pt() << " eta: \t" << lept.Eta() << " phi: \t" << lept.Phi() << std::endl;
        std::cout << "+ qJet: pt: \t" << qJet.Pt() << " eta: \t" << qJet.Eta() << " phi: \t" << qJet.Phi() << std::endl;
        std::cout << "njets             : " << njets_bdt              << std::endl;
        std::cout << "deltaPhi_top_higgs: " << deltaPhi_top_higgs_bdt << std::endl;
        std::cout << "mt_top            : " << mt_top_bdt             << std::endl;
        std::cout << "eta_qJet          : " << eta_qJet_bdt           << std::endl;
        std::cout << "charge_lept       : " << charge_lept_bdt        << std::endl;
        std::cout << "deltaEta_lept_qJet: " << deltaEta_lept_qJet_bdt << std::endl;
        std::cout << "BDT: " << BDT_lept_t << std::endl;
        std::cout << "mgg: " << diphot.M() << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << std::endl << std::endl;



      } else {


        // *** HADRONIC CHANNEL

        if( qJet.Pt()< 45. ) continue;

        if( nJets < 4 ) continue;


        float maxDeltaPhi=0;
        int i_jetW1=-1;
        int i_jetW2=-1;
        
        // define top candidate as the one maximising deltaPhi(top, higgs):
        for( unsigned i=0; i<jets.size(); ++i ) {
          for( unsigned j=i+1; j<jets.size(); ++j ) {

            if( (i+1)==i_bJet || (j+1)==i_bJet ) continue;
            if( (i+1)==i_qJet || (j+1)==i_qJet ) continue;


            TLorentzVector jetW1_tmp = jets[i];
            TLorentzVector jetW2_tmp = jets[j];

            TLorentzVector top_tmp = jetW1_tmp+jetW2_tmp+bJet;

            float thisDeltaPhi = fabs(top_tmp.DeltaPhi(diphot));
            if( thisDeltaPhi>maxDeltaPhi ) {
              maxDeltaPhi=thisDeltaPhi;
              i_jetW1 = i;
              i_jetW2 = j;
            }

          } //j
        } //i


        if( i_jetW1<0 || i_jetW2<0 ) continue;

        TLorentzVector jetW1 = jets[i_jetW1];
        TLorentzVector jetW2 = jets[i_jetW2];

        if( fabs(jetW1.Eta())<2. ) nCentralJets--;
        if( fabs(jetW2.Eta())<2. ) nCentralJets--;


        TLorentzVector W = jetW1 + jetW2;

        TLorentzVector top = bJet + W;

        float Wmass = 80.4;

        if( fabs(top.M()-172.5)>40. ) continue;


        std::cout << std::endl << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << "*** Event passing HADRONIC channel selection: run: " << run << " event: " << event << std::endl;
        std::cout << "+ phot1: pt: \t" << phot1.Pt() << " eta: \t" << phot1.Eta() << " phi: \t" << phot1.Phi() << std::endl;
        std::cout << "+ phot2: pt: \t" << phot2.Pt() << " eta: \t" << phot2.Eta() << " phi: \t" << phot2.Phi() << std::endl;
        std::cout << "+ qJet: pt: \t" << qJet.Pt() << " eta: \t" << qJet.Eta() << " phi: \t" << qJet.Phi() << std::endl;
        std::cout << "+ jetW1: pt: \t" << jetW1.Pt() << " eta: \t" << jetW1.Eta() << " phi: \t" << jetW1.Phi() << std::endl;
        std::cout << "+ jetW2: pt: \t" << jetW2.Pt() << " eta: \t" << jetW2.Eta() << " phi: \t" << jetW2.Phi() << std::endl;
        std::cout << "nJets             : " << nJets              << std::endl;
        std::cout << "nBJetsLoose       : " << nBJetsLoose          << std::endl;
        std::cout << "nBJetsMedium      : " << nBJetsMedium          << std::endl;
        std::cout << "nCentralJets      : " << nCentralJets              << std::endl;
        std::cout << "m_top             : " << top.M()             << std::endl;
        std::cout << "m_W               : " << W.M()             << std::endl;
        std::cout << "mgg: " << diphot.M() << std::endl;
        std::cout << "-------------------------------------------------------------------------------" << std::endl;
        std::cout << std::endl << std::endl;


      } // leptonic/hadronic channels


      tree->Fill();


   }


   std::cout << "Events passing selections: " <<std::endl;
   std::cout << " + Hadronic Channel: " << tree->GetEntries("!isLeptonic") << " (baseline)   \t" << tree->GetEntries("!isLeptonic && nCentralJets<2") << " (plus CJV)" << std::endl;
   std::cout << " + Leptonic Channel: " << tree->GetEntries("isLeptonic") << " (baseline)  \t" << tree->GetEntries("isLeptonic && BDT_lept>0.2") << " (plus BDT)" << std::endl;

   outfile->cd();

   tree->Write();

   outfile->Close();

}
