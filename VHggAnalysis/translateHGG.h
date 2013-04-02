//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 20 11:49:00 2013 by ROOT version 5.32/00
// from TChain Data/
//////////////////////////////////////////////////////////

#ifndef translateHGG_h
#define translateHGG_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class translateHGG {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           itype;
   Float_t         run;
   Int_t           event;
   Float_t         lumis;
   Float_t         weight;
   Float_t         evweight;
   Float_t         pu_weight;
   Float_t         pu_n;
   Float_t         nvtx;
   Float_t         rho;
   Int_t           category;
   Float_t         ph1_e;
   Float_t         ph2_e;
   Float_t         ph1_pt;
   Float_t         ph2_pt;
   Float_t         ph1_phi;
   Float_t         ph2_phi;
   Float_t         ph1_eta;
   Float_t         ph2_eta;
   Float_t         ph1_r9;
   Float_t         ph2_r9;
   Int_t           ph1_isPrompt;
   Int_t           ph2_isPrompt;
   Float_t         ph1_SCEta;
   Float_t         ph2_SCEta;
   Float_t         ph1_SCPhi;
   Float_t         ph2_SCPhi;
   Float_t         ph1_hoe;
   Float_t         ph2_hoe;
   Float_t         ph1_sieie;
   Float_t         ph2_sieie;
   Float_t         ph1_pfchargedisogood03;
   Float_t         ph2_pfchargedisogood03;
   Float_t         ph1_pfchargedisobad04;
   Float_t         ph2_pfchargedisobad04;
   Float_t         ph1_etawidth;
   Float_t         ph2_etawidth;
   Float_t         ph1_phiwidth;
   Float_t         ph2_phiwidth;
   Float_t         ph1_eseffssqrt;
   Float_t         ph2_eseffssqrt;
   Float_t         ph1_pfchargedisobad03;
   Float_t         ph2_pfchargedisobad03;
   Float_t         ph1_sieip;
   Float_t         ph2_sieip;
   Float_t         ph1_sipip;
   Float_t         ph2_sipip;
   Float_t         ph1_ecaliso;
   Float_t         ph2_ecaliso;
   Float_t         ph1_ecalisobad;
   Float_t         ph2_ecalisobad;
   Float_t         ph1_badvtx_Et;
   Float_t         ph2_badvtx_Et;
   Float_t         ph1_isconv;
   Float_t         ph2_isconv;
   Int_t           ph1_ciclevel;
   Int_t           ph2_ciclevel;
   Float_t         ph1_sigmaEoE;
   Float_t         ph2_sigmaEoE;
   Float_t         ph1_ptoM;
   Float_t         ph2_ptoM;
   Int_t           ph1_isEB;
   Int_t           ph2_isEB;
   Float_t         ph1_s4ratio;
   Float_t         ph2_s4ratio;
   Float_t         ph1_e3x3;
   Float_t         ph2_e3x3;
   Float_t         ph1_e5x5;
   Float_t         ph2_e5x5;
   Float_t         PhotonsMass;
   Float_t         dipho_E;
   Float_t         dipho_pt;
   Float_t         dipho_eta;
   Float_t         dipho_phi;
   Float_t         dipho_cosThetaStar_CS;
   Float_t         dipho_tanhYStar;
   Float_t         dipho_Y;
   Int_t           vtx_ind;
   Float_t         vtx_x;
   Float_t         vtx_y;
   Float_t         vtx_z;
   Float_t         vtx_mva;
   Float_t         vtx_mva_2;
   Float_t         vtx_mva_3;
   Float_t         vtx_ptbal;
   Float_t         vtx_ptasym;
   Float_t         vtx_logsumpt2;
   Float_t         vtx_pulltoconv;
   Float_t         vtx_prob;
   Int_t           njets_passing_kLooseID;
   Float_t         j1_e;
   Float_t         j1_pt;
   Float_t         j1_phi;
   Float_t         j1_eta;
   Float_t         j1_beta;
   Float_t         j1_betaStar;
   Float_t         j1_betaStarClassic;
   Float_t         j1_dR2Mean;
   Float_t         j1_algoPF1_csvBtag;
   Float_t         j2_e;
   Float_t         j2_pt;
   Float_t         j2_phi;
   Float_t         j2_eta;
   Float_t         j2_beta;
   Float_t         j2_betaStar;
   Float_t         j2_betaStarClassic;
   Float_t         j2_dR2Mean;
   Float_t         j2_algoPF1_csvBtag;
   Float_t         j3_e;
   Float_t         j3_pt;
   Float_t         j3_phi;
   Float_t         j3_eta;
   Float_t         j3_beta;
   Float_t         j3_betaStar;
   Float_t         j3_betaStarClassic;
   Float_t         j3_dR2Mean;
   Float_t         j3_algoPF1_csvBtag;
   Float_t         j4_e;
   Float_t         j4_pt;
   Float_t         j4_phi;
   Float_t         j4_eta;
   Float_t         j4_beta;
   Float_t         j4_betaStar;
   Float_t         j4_betaStarClassic;
   Float_t         j4_dR2Mean;
   Float_t         j4_algoPF1_csvBtag;
   Float_t         j5_e;
   Float_t         j5_pt;
   Float_t         j5_phi;
   Float_t         j5_eta;
   Float_t         j5_beta;
   Float_t         j5_betaStar;
   Float_t         j5_betaStarClassic;
   Float_t         j5_dR2Mean;
   Float_t         j5_algoPF1_csvBtag;
   Float_t         j6_e;
   Float_t         j6_pt;
   Float_t         j6_phi;
   Float_t         j6_eta;
   Float_t         j6_beta;
   Float_t         j6_betaStar;
   Float_t         j6_betaStarClassic;
   Float_t         j6_dR2Mean;
   Float_t         j6_algoPF1_csvBtag;
   Float_t         j7_e;
   Float_t         j7_pt;
   Float_t         j7_phi;
   Float_t         j7_eta;
   Float_t         j7_beta;
   Float_t         j7_betaStar;
   Float_t         j7_betaStarClassic;
   Float_t         j7_dR2Mean;
   Float_t         j7_algoPF1_csvBtag;
   Float_t         j8_e;
   Float_t         j8_pt;
   Float_t         j8_phi;
   Float_t         j8_eta;
   Float_t         j8_beta;
   Float_t         j8_betaStar;
   Float_t         j8_betaStarClassic;
   Float_t         j8_dR2Mean;
   Float_t         j8_algoPF1_csvBtag;
   Float_t         j9_e;
   Float_t         j9_pt;
   Float_t         j9_phi;
   Float_t         j9_eta;
   Float_t         j9_beta;
   Float_t         j9_betaStar;
   Float_t         j9_betaStarClassic;
   Float_t         j9_dR2Mean;
   Float_t         j9_algoPF1_csvBtag;
   Float_t         j10_e;
   Float_t         j10_pt;
   Float_t         j10_phi;
   Float_t         j10_eta;
   Float_t         j10_beta;
   Float_t         j10_betaStar;
   Float_t         j10_betaStarClassic;
   Float_t         j10_dR2Mean;
   Float_t         j10_algoPF1_csvBtag;
   Float_t         JetsMass;
   Float_t         dijet_E;
   Float_t         dijet_Pt;
   Float_t         dijet_Eta;
   Float_t         dijet_Phi;
   Int_t           nelectrons;
   Int_t           nmuons;
   Int_t           chargeEle;
   Int_t           chargeMu;
   Float_t         ptEle;
   Float_t         etaEle;
   Float_t         phiEle;
   Float_t         ptMu;
   Float_t         etaMu;
   Float_t         phiMu;
   Float_t         met_pfmet;
   Float_t         met_phi_pfmet;

   // List of branches
   TBranch        *b_itype;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_evweight;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_pu_n;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_category;   //!
   TBranch        *b_ph1_e;   //!
   TBranch        *b_ph2_e;   //!
   TBranch        *b_ph1_pt;   //!
   TBranch        *b_ph2_pt;   //!
   TBranch        *b_ph1_phi;   //!
   TBranch        *b_ph2_phi;   //!
   TBranch        *b_ph1_eta;   //!
   TBranch        *b_ph2_eta;   //!
   TBranch        *b_ph1_r9;   //!
   TBranch        *b_ph2_r9;   //!
   TBranch        *b_ph1_isPrompt;   //!
   TBranch        *b_ph2_isPrompt;   //!
   TBranch        *b_ph1_SCEta;   //!
   TBranch        *b_ph2_SCEta;   //!
   TBranch        *b_ph1_SCPhi;   //!
   TBranch        *b_ph2_SCPhi;   //!
   TBranch        *b_ph1_hoe;   //!
   TBranch        *b_ph2_hoe;   //!
   TBranch        *b_ph1_sieie;   //!
   TBranch        *b_ph2_sieie;   //!
   TBranch        *b_ph1_pfchargedisogood03;   //!
   TBranch        *b_ph2_pfchargedisogood03;   //!
   TBranch        *b_ph1_pfchargedisobad04;   //!
   TBranch        *b_ph2_pfchargedisobad04;   //!
   TBranch        *b_ph1_etawidth;   //!
   TBranch        *b_ph2_etawidth;   //!
   TBranch        *b_ph1_phiwidth;   //!
   TBranch        *b_ph2_phiwidth;   //!
   TBranch        *b_ph1_eseffssqrt;   //!
   TBranch        *b_ph2_eseffssqrt;   //!
   TBranch        *b_ph1_pfchargedisobad03;   //!
   TBranch        *b_ph2_pfchargedisobad03;   //!
   TBranch        *b_ph1_sieip;   //!
   TBranch        *b_ph2_sieip;   //!
   TBranch        *b_ph1_sipip;   //!
   TBranch        *b_ph2_sipip;   //!
   TBranch        *b_ph1_ecaliso;   //!
   TBranch        *b_ph2_ecaliso;   //!
   TBranch        *b_ph1_ecalisobad;   //!
   TBranch        *b_ph2_ecalisobad;   //!
   TBranch        *b_ph1_badvtx_Et;   //!
   TBranch        *b_ph2_badvtx_Et;   //!
   TBranch        *b_ph1_isconv;   //!
   TBranch        *b_ph2_isconv;   //!
   TBranch        *b_ph1_ciclevel;   //!
   TBranch        *b_ph2_ciclevel;   //!
   TBranch        *b_ph1_sigmaEoE;   //!
   TBranch        *b_ph2_sigmaEoE;   //!
   TBranch        *b_ph1_ptoM;   //!
   TBranch        *b_ph2_ptoM;   //!
   TBranch        *b_ph1_isEB;   //!
   TBranch        *b_ph2_isEB;   //!
   TBranch        *b_ph1_s4ratio;   //!
   TBranch        *b_ph2_s4ratio;   //!
   TBranch        *b_ph1_e3x3;   //!
   TBranch        *b_ph2_e3x3;   //!
   TBranch        *b_ph1_e5x5;   //!
   TBranch        *b_ph2_e5x5;   //!
   TBranch        *b_PhotonsMass;   //!
   TBranch        *b_dipho_E;   //!
   TBranch        *b_dipho_pt;   //!
   TBranch        *b_dipho_eta;   //!
   TBranch        *b_dipho_phi;   //!
   TBranch        *b_dipho_cosThetaStar_CS;   //!
   TBranch        *b_dipho_tanhYStar;   //!
   TBranch        *b_dipho_Y;   //!
   TBranch        *b_vtx_ind;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_mva;   //!
   TBranch        *b_vtx_mva_2;   //!
   TBranch        *b_vtx_mva_3;   //!
   TBranch        *b_vtx_ptbal;   //!
   TBranch        *b_vtx_ptasym;   //!
   TBranch        *b_vtx_logsumpt2;   //!
   TBranch        *b_vtx_pulltoconv;   //!
   TBranch        *b_vtx_prob;   //!
   TBranch        *b_njets_passing_kLooseID;   //!
   TBranch        *b_j1_e;   //!
   TBranch        *b_j1_pt;   //!
   TBranch        *b_j1_phi;   //!
   TBranch        *b_j1_eta;   //!
   TBranch        *b_j1_beta;   //!
   TBranch        *b_j1_betaStar;   //!
   TBranch        *b_j1_betaStarClassic;   //!
   TBranch        *b_j1_dR2Mean;   //!
   TBranch        *b_j1_algoPF1_csvBtag;   //!
   TBranch        *b_j2_e;   //!
   TBranch        *b_j2_pt;   //!
   TBranch        *b_j2_phi;   //!
   TBranch        *b_j2_eta;   //!
   TBranch        *b_j2_beta;   //!
   TBranch        *b_j2_betaStar;   //!
   TBranch        *b_j2_betaStarClassic;   //!
   TBranch        *b_j2_dR2Mean;   //!
   TBranch        *b_j2_algoPF1_csvBtag;   //!
   TBranch        *b_j3_e;   //!
   TBranch        *b_j3_pt;   //!
   TBranch        *b_j3_phi;   //!
   TBranch        *b_j3_eta;   //!
   TBranch        *b_j3_beta;   //!
   TBranch        *b_j3_betaStar;   //!
   TBranch        *b_j3_betaStarClassic;   //!
   TBranch        *b_j3_dR2Mean;   //!
   TBranch        *b_j3_algoPF1_csvBtag;   //!
   TBranch        *b_j4_e;   //!
   TBranch        *b_j4_pt;   //!
   TBranch        *b_j4_phi;   //!
   TBranch        *b_j4_eta;   //!
   TBranch        *b_j4_beta;   //!
   TBranch        *b_j4_betaStar;   //!
   TBranch        *b_j4_betaStarClassic;   //!
   TBranch        *b_j4_dR2Mean;   //!
   TBranch        *b_j4_algoPF1_csvBtag;   //!
   TBranch        *b_j5_e;   //!
   TBranch        *b_j5_pt;   //!
   TBranch        *b_j5_phi;   //!
   TBranch        *b_j5_eta;   //!
   TBranch        *b_j5_beta;   //!
   TBranch        *b_j5_betaStar;   //!
   TBranch        *b_j5_betaStarClassic;   //!
   TBranch        *b_j5_dR2Mean;   //!
   TBranch        *b_j5_algoPF1_csvBtag;   //!
   TBranch        *b_j6_e;   //!
   TBranch        *b_j6_pt;   //!
   TBranch        *b_j6_phi;   //!
   TBranch        *b_j6_eta;   //!
   TBranch        *b_j6_beta;   //!
   TBranch        *b_j6_betaStar;   //!
   TBranch        *b_j6_betaStarClassic;   //!
   TBranch        *b_j6_dR2Mean;   //!
   TBranch        *b_j6_algoPF1_csvBtag;   //!
   TBranch        *b_j7_e;   //!
   TBranch        *b_j7_pt;   //!
   TBranch        *b_j7_phi;   //!
   TBranch        *b_j7_eta;   //!
   TBranch        *b_j7_beta;   //!
   TBranch        *b_j7_betaStar;   //!
   TBranch        *b_j7_betaStarClassic;   //!
   TBranch        *b_j7_dR2Mean;   //!
   TBranch        *b_j7_algoPF1_csvBtag;   //!
   TBranch        *b_j8_e;   //!
   TBranch        *b_j8_pt;   //!
   TBranch        *b_j8_phi;   //!
   TBranch        *b_j8_eta;   //!
   TBranch        *b_j8_beta;   //!
   TBranch        *b_j8_betaStar;   //!
   TBranch        *b_j8_betaStarClassic;   //!
   TBranch        *b_j8_dR2Mean;   //!
   TBranch        *b_j8_algoPF1_csvBtag;   //!
   TBranch        *b_j9_e;   //!
   TBranch        *b_j9_pt;   //!
   TBranch        *b_j9_phi;   //!
   TBranch        *b_j9_eta;   //!
   TBranch        *b_j9_beta;   //!
   TBranch        *b_j9_betaStar;   //!
   TBranch        *b_j9_betaStarClassic;   //!
   TBranch        *b_j9_dR2Mean;   //!
   TBranch        *b_j9_algoPF1_csvBtag;   //!
   TBranch        *b_j10_e;   //!
   TBranch        *b_j10_pt;   //!
   TBranch        *b_j10_phi;   //!
   TBranch        *b_j10_eta;   //!
   TBranch        *b_j10_beta;   //!
   TBranch        *b_j10_betaStar;   //!
   TBranch        *b_j10_betaStarClassic;   //!
   TBranch        *b_j10_dR2Mean;   //!
   TBranch        *b_j10_algoPF1_csvBtag;   //!
   TBranch        *b_JetsMass;   //!
   TBranch        *b_dijet_E;   //!
   TBranch        *b_dijet_Pt;   //!
   TBranch        *b_dijet_Eta;   //!
   TBranch        *b_dijet_Phi;   //!
   TBranch        *b_nelectrons;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_chargeEle;   //!
   TBranch        *b_chargeMu;   //!
   TBranch        *b_ptEle;   //!
   TBranch        *b_etaEle;   //!
   TBranch        *b_phiEle;   //!
   TBranch        *b_ptMu;   //!
   TBranch        *b_etaMu;   //!
   TBranch        *b_phiMu;   //!
   TBranch        *b_met_pfmet;   //!
   TBranch        *b_met_phi_pfmet;   //!

   translateHGG(TTree *tree=0);
   virtual ~translateHGG();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef translateHGG_cxx
translateHGG::translateHGG(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/p/pandolf/public/histograms_CMS-HGG_OLD.root");
      TFile *f ;
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/work/p/pandolf/public/histograms_CMS-HGG_NEW.root");
      }
      f->GetObject("Data",tree);

//#else // SINGLE_TREE
//
//      // The following code should be used if you want this class to access a chain
//      // of trees.
//      TChain * chain = new TChain("Data","");
//      chain->Add("/afs/cern.ch/work/p/pandolf/public/histograms_CMS-HGG.root/Data");
//      tree = chain;
//#endif // SINGLE_TREE

   }
   Init(tree);
}

translateHGG::~translateHGG()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t translateHGG::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t translateHGG::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void translateHGG::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("itype", &itype, &b_itype);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("evweight", &evweight, &b_evweight);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("pu_n", &pu_n, &b_pu_n);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("ph1_e", &ph1_e, &b_ph1_e);
   fChain->SetBranchAddress("ph2_e", &ph2_e, &b_ph2_e);
   fChain->SetBranchAddress("ph1_pt", &ph1_pt, &b_ph1_pt);
   fChain->SetBranchAddress("ph2_pt", &ph2_pt, &b_ph2_pt);
   fChain->SetBranchAddress("ph1_phi", &ph1_phi, &b_ph1_phi);
   fChain->SetBranchAddress("ph2_phi", &ph2_phi, &b_ph2_phi);
   fChain->SetBranchAddress("ph1_eta", &ph1_eta, &b_ph1_eta);
   fChain->SetBranchAddress("ph2_eta", &ph2_eta, &b_ph2_eta);
   fChain->SetBranchAddress("ph1_r9", &ph1_r9, &b_ph1_r9);
   fChain->SetBranchAddress("ph2_r9", &ph2_r9, &b_ph2_r9);
   fChain->SetBranchAddress("ph1_isPrompt", &ph1_isPrompt, &b_ph1_isPrompt);
   fChain->SetBranchAddress("ph2_isPrompt", &ph2_isPrompt, &b_ph2_isPrompt);
   fChain->SetBranchAddress("ph1_SCEta", &ph1_SCEta, &b_ph1_SCEta);
   fChain->SetBranchAddress("ph2_SCEta", &ph2_SCEta, &b_ph2_SCEta);
   fChain->SetBranchAddress("ph1_SCPhi", &ph1_SCPhi, &b_ph1_SCPhi);
   fChain->SetBranchAddress("ph2_SCPhi", &ph2_SCPhi, &b_ph2_SCPhi);
   fChain->SetBranchAddress("ph1_hoe", &ph1_hoe, &b_ph1_hoe);
   fChain->SetBranchAddress("ph2_hoe", &ph2_hoe, &b_ph2_hoe);
   fChain->SetBranchAddress("ph1_sieie", &ph1_sieie, &b_ph1_sieie);
   fChain->SetBranchAddress("ph2_sieie", &ph2_sieie, &b_ph2_sieie);
   fChain->SetBranchAddress("ph1_pfchargedisogood03", &ph1_pfchargedisogood03, &b_ph1_pfchargedisogood03);
   fChain->SetBranchAddress("ph2_pfchargedisogood03", &ph2_pfchargedisogood03, &b_ph2_pfchargedisogood03);
   fChain->SetBranchAddress("ph1_pfchargedisobad04", &ph1_pfchargedisobad04, &b_ph1_pfchargedisobad04);
   fChain->SetBranchAddress("ph2_pfchargedisobad04", &ph2_pfchargedisobad04, &b_ph2_pfchargedisobad04);
   fChain->SetBranchAddress("ph1_etawidth", &ph1_etawidth, &b_ph1_etawidth);
   fChain->SetBranchAddress("ph2_etawidth", &ph2_etawidth, &b_ph2_etawidth);
   fChain->SetBranchAddress("ph1_phiwidth", &ph1_phiwidth, &b_ph1_phiwidth);
   fChain->SetBranchAddress("ph2_phiwidth", &ph2_phiwidth, &b_ph2_phiwidth);
   fChain->SetBranchAddress("ph1_eseffssqrt", &ph1_eseffssqrt, &b_ph1_eseffssqrt);
   fChain->SetBranchAddress("ph2_eseffssqrt", &ph2_eseffssqrt, &b_ph2_eseffssqrt);
   fChain->SetBranchAddress("ph1_pfchargedisobad03", &ph1_pfchargedisobad03, &b_ph1_pfchargedisobad03);
   fChain->SetBranchAddress("ph2_pfchargedisobad03", &ph2_pfchargedisobad03, &b_ph2_pfchargedisobad03);
   fChain->SetBranchAddress("ph1_sieip", &ph1_sieip, &b_ph1_sieip);
   fChain->SetBranchAddress("ph2_sieip", &ph2_sieip, &b_ph2_sieip);
   fChain->SetBranchAddress("ph1_sipip", &ph1_sipip, &b_ph1_sipip);
   fChain->SetBranchAddress("ph2_sipip", &ph2_sipip, &b_ph2_sipip);
   fChain->SetBranchAddress("ph1_ecaliso", &ph1_ecaliso, &b_ph1_ecaliso);
   fChain->SetBranchAddress("ph2_ecaliso", &ph2_ecaliso, &b_ph2_ecaliso);
   fChain->SetBranchAddress("ph1_ecalisobad", &ph1_ecalisobad, &b_ph1_ecalisobad);
   fChain->SetBranchAddress("ph2_ecalisobad", &ph2_ecalisobad, &b_ph2_ecalisobad);
   fChain->SetBranchAddress("ph1_badvtx_Et", &ph1_badvtx_Et, &b_ph1_badvtx_Et);
   fChain->SetBranchAddress("ph2_badvtx_Et", &ph2_badvtx_Et, &b_ph2_badvtx_Et);
   fChain->SetBranchAddress("ph1_isconv", &ph1_isconv, &b_ph1_isconv);
   fChain->SetBranchAddress("ph2_isconv", &ph2_isconv, &b_ph2_isconv);
   fChain->SetBranchAddress("ph1_ciclevel", &ph1_ciclevel, &b_ph1_ciclevel);
   fChain->SetBranchAddress("ph2_ciclevel", &ph2_ciclevel, &b_ph2_ciclevel);
   fChain->SetBranchAddress("ph1_sigmaEoE", &ph1_sigmaEoE, &b_ph1_sigmaEoE);
   fChain->SetBranchAddress("ph2_sigmaEoE", &ph2_sigmaEoE, &b_ph2_sigmaEoE);
   fChain->SetBranchAddress("ph1_ptoM", &ph1_ptoM, &b_ph1_ptoM);
   fChain->SetBranchAddress("ph2_ptoM", &ph2_ptoM, &b_ph2_ptoM);
   fChain->SetBranchAddress("ph1_isEB", &ph1_isEB, &b_ph1_isEB);
   fChain->SetBranchAddress("ph2_isEB", &ph2_isEB, &b_ph2_isEB);
   fChain->SetBranchAddress("ph1_s4ratio", &ph1_s4ratio, &b_ph1_s4ratio);
   fChain->SetBranchAddress("ph2_s4ratio", &ph2_s4ratio, &b_ph2_s4ratio);
   fChain->SetBranchAddress("ph1_e3x3", &ph1_e3x3, &b_ph1_e3x3);
   fChain->SetBranchAddress("ph2_e3x3", &ph2_e3x3, &b_ph2_e3x3);
   fChain->SetBranchAddress("ph1_e5x5", &ph1_e5x5, &b_ph1_e5x5);
   fChain->SetBranchAddress("ph2_e5x5", &ph2_e5x5, &b_ph2_e5x5);
   fChain->SetBranchAddress("PhotonsMass", &PhotonsMass, &b_PhotonsMass);
   fChain->SetBranchAddress("dipho_E", &dipho_E, &b_dipho_E);
   fChain->SetBranchAddress("dipho_pt", &dipho_pt, &b_dipho_pt);
   fChain->SetBranchAddress("dipho_eta", &dipho_eta, &b_dipho_eta);
   fChain->SetBranchAddress("dipho_phi", &dipho_phi, &b_dipho_phi);
   fChain->SetBranchAddress("dipho_cosThetaStar_CS", &dipho_cosThetaStar_CS, &b_dipho_cosThetaStar_CS);
   fChain->SetBranchAddress("dipho_tanhYStar", &dipho_tanhYStar, &b_dipho_tanhYStar);
   fChain->SetBranchAddress("dipho_Y", &dipho_Y, &b_dipho_Y);
   fChain->SetBranchAddress("vtx_ind", &vtx_ind, &b_vtx_ind);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_mva", &vtx_mva, &b_vtx_mva);
   fChain->SetBranchAddress("vtx_mva_2", &vtx_mva_2, &b_vtx_mva_2);
   fChain->SetBranchAddress("vtx_mva_3", &vtx_mva_3, &b_vtx_mva_3);
   fChain->SetBranchAddress("vtx_ptbal", &vtx_ptbal, &b_vtx_ptbal);
   fChain->SetBranchAddress("vtx_ptasym", &vtx_ptasym, &b_vtx_ptasym);
   fChain->SetBranchAddress("vtx_logsumpt2", &vtx_logsumpt2, &b_vtx_logsumpt2);
   fChain->SetBranchAddress("vtx_pulltoconv", &vtx_pulltoconv, &b_vtx_pulltoconv);
   fChain->SetBranchAddress("vtx_prob", &vtx_prob, &b_vtx_prob);
   fChain->SetBranchAddress("njets_passing_kLooseID", &njets_passing_kLooseID, &b_njets_passing_kLooseID);
   fChain->SetBranchAddress("j1_e", &j1_e, &b_j1_e);
   fChain->SetBranchAddress("j1_pt", &j1_pt, &b_j1_pt);
   fChain->SetBranchAddress("j1_phi", &j1_phi, &b_j1_phi);
   fChain->SetBranchAddress("j1_eta", &j1_eta, &b_j1_eta);
   fChain->SetBranchAddress("j1_beta", &j1_beta, &b_j1_beta);
   fChain->SetBranchAddress("j1_betaStar", &j1_betaStar, &b_j1_betaStar);
   fChain->SetBranchAddress("j1_betaStarClassic", &j1_betaStarClassic, &b_j1_betaStarClassic);
   fChain->SetBranchAddress("j1_dR2Mean", &j1_dR2Mean, &b_j1_dR2Mean);
   fChain->SetBranchAddress("j1_algoPF1_csvBtag", &j1_algoPF1_csvBtag, &b_j1_algoPF1_csvBtag);
   fChain->SetBranchAddress("j2_e", &j2_e, &b_j2_e);
   fChain->SetBranchAddress("j2_pt", &j2_pt, &b_j2_pt);
   fChain->SetBranchAddress("j2_phi", &j2_phi, &b_j2_phi);
   fChain->SetBranchAddress("j2_eta", &j2_eta, &b_j2_eta);
   fChain->SetBranchAddress("j2_beta", &j2_beta, &b_j2_beta);
   fChain->SetBranchAddress("j2_betaStar", &j2_betaStar, &b_j2_betaStar);
   fChain->SetBranchAddress("j2_betaStarClassic", &j2_betaStarClassic, &b_j2_betaStarClassic);
   fChain->SetBranchAddress("j2_dR2Mean", &j2_dR2Mean, &b_j2_dR2Mean);
   fChain->SetBranchAddress("j2_algoPF1_csvBtag", &j2_algoPF1_csvBtag, &b_j2_algoPF1_csvBtag);
   fChain->SetBranchAddress("j3_e", &j3_e, &b_j3_e);
   fChain->SetBranchAddress("j3_pt", &j3_pt, &b_j3_pt);
   fChain->SetBranchAddress("j3_phi", &j3_phi, &b_j3_phi);
   fChain->SetBranchAddress("j3_eta", &j3_eta, &b_j3_eta);
   fChain->SetBranchAddress("j3_beta", &j3_beta, &b_j3_beta);
   fChain->SetBranchAddress("j3_betaStar", &j3_betaStar, &b_j3_betaStar);
   fChain->SetBranchAddress("j3_betaStarClassic", &j3_betaStarClassic, &b_j3_betaStarClassic);
   fChain->SetBranchAddress("j3_dR2Mean", &j3_dR2Mean, &b_j3_dR2Mean);
   fChain->SetBranchAddress("j3_algoPF1_csvBtag", &j3_algoPF1_csvBtag, &b_j3_algoPF1_csvBtag);
   fChain->SetBranchAddress("j4_e", &j4_e, &b_j4_e);
   fChain->SetBranchAddress("j4_pt", &j4_pt, &b_j4_pt);
   fChain->SetBranchAddress("j4_phi", &j4_phi, &b_j4_phi);
   fChain->SetBranchAddress("j4_eta", &j4_eta, &b_j4_eta);
   fChain->SetBranchAddress("j4_beta", &j4_beta, &b_j4_beta);
   fChain->SetBranchAddress("j4_betaStar", &j4_betaStar, &b_j4_betaStar);
   fChain->SetBranchAddress("j4_betaStarClassic", &j4_betaStarClassic, &b_j4_betaStarClassic);
   fChain->SetBranchAddress("j4_dR2Mean", &j4_dR2Mean, &b_j4_dR2Mean);
   fChain->SetBranchAddress("j4_algoPF1_csvBtag", &j4_algoPF1_csvBtag, &b_j4_algoPF1_csvBtag);
   fChain->SetBranchAddress("j5_e", &j5_e, &b_j5_e);
   fChain->SetBranchAddress("j5_pt", &j5_pt, &b_j5_pt);
   fChain->SetBranchAddress("j5_phi", &j5_phi, &b_j5_phi);
   fChain->SetBranchAddress("j5_eta", &j5_eta, &b_j5_eta);
   fChain->SetBranchAddress("j5_beta", &j5_beta, &b_j5_beta);
   fChain->SetBranchAddress("j5_betaStar", &j5_betaStar, &b_j5_betaStar);
   fChain->SetBranchAddress("j5_betaStarClassic", &j5_betaStarClassic, &b_j5_betaStarClassic);
   fChain->SetBranchAddress("j5_dR2Mean", &j5_dR2Mean, &b_j5_dR2Mean);
   fChain->SetBranchAddress("j5_algoPF1_csvBtag", &j5_algoPF1_csvBtag, &b_j5_algoPF1_csvBtag);
   fChain->SetBranchAddress("j6_e", &j6_e, &b_j6_e);
   fChain->SetBranchAddress("j6_pt", &j6_pt, &b_j6_pt);
   fChain->SetBranchAddress("j6_phi", &j6_phi, &b_j6_phi);
   fChain->SetBranchAddress("j6_eta", &j6_eta, &b_j6_eta);
   fChain->SetBranchAddress("j6_beta", &j6_beta, &b_j6_beta);
   fChain->SetBranchAddress("j6_betaStar", &j6_betaStar, &b_j6_betaStar);
   fChain->SetBranchAddress("j6_betaStarClassic", &j6_betaStarClassic, &b_j6_betaStarClassic);
   fChain->SetBranchAddress("j6_dR2Mean", &j6_dR2Mean, &b_j6_dR2Mean);
   fChain->SetBranchAddress("j6_algoPF1_csvBtag", &j6_algoPF1_csvBtag, &b_j6_algoPF1_csvBtag);
   fChain->SetBranchAddress("j7_e", &j7_e, &b_j7_e);
   fChain->SetBranchAddress("j7_pt", &j7_pt, &b_j7_pt);
   fChain->SetBranchAddress("j7_phi", &j7_phi, &b_j7_phi);
   fChain->SetBranchAddress("j7_eta", &j7_eta, &b_j7_eta);
   fChain->SetBranchAddress("j7_beta", &j7_beta, &b_j7_beta);
   fChain->SetBranchAddress("j7_betaStar", &j7_betaStar, &b_j7_betaStar);
   fChain->SetBranchAddress("j7_betaStarClassic", &j7_betaStarClassic, &b_j7_betaStarClassic);
   fChain->SetBranchAddress("j7_dR2Mean", &j7_dR2Mean, &b_j7_dR2Mean);
   fChain->SetBranchAddress("j7_algoPF1_csvBtag", &j7_algoPF1_csvBtag, &b_j7_algoPF1_csvBtag);
   fChain->SetBranchAddress("j8_e", &j8_e, &b_j8_e);
   fChain->SetBranchAddress("j8_pt", &j8_pt, &b_j8_pt);
   fChain->SetBranchAddress("j8_phi", &j8_phi, &b_j8_phi);
   fChain->SetBranchAddress("j8_eta", &j8_eta, &b_j8_eta);
   fChain->SetBranchAddress("j8_beta", &j8_beta, &b_j8_beta);
   fChain->SetBranchAddress("j8_betaStar", &j8_betaStar, &b_j8_betaStar);
   fChain->SetBranchAddress("j8_betaStarClassic", &j8_betaStarClassic, &b_j8_betaStarClassic);
   fChain->SetBranchAddress("j8_dR2Mean", &j8_dR2Mean, &b_j8_dR2Mean);
   fChain->SetBranchAddress("j8_algoPF1_csvBtag", &j8_algoPF1_csvBtag, &b_j8_algoPF1_csvBtag);
   fChain->SetBranchAddress("j9_e", &j9_e, &b_j9_e);
   fChain->SetBranchAddress("j9_pt", &j9_pt, &b_j9_pt);
   fChain->SetBranchAddress("j9_phi", &j9_phi, &b_j9_phi);
   fChain->SetBranchAddress("j9_eta", &j9_eta, &b_j9_eta);
   fChain->SetBranchAddress("j9_beta", &j9_beta, &b_j9_beta);
   fChain->SetBranchAddress("j9_betaStar", &j9_betaStar, &b_j9_betaStar);
   fChain->SetBranchAddress("j9_betaStarClassic", &j9_betaStarClassic, &b_j9_betaStarClassic);
   fChain->SetBranchAddress("j9_dR2Mean", &j9_dR2Mean, &b_j9_dR2Mean);
   fChain->SetBranchAddress("j9_algoPF1_csvBtag", &j9_algoPF1_csvBtag, &b_j9_algoPF1_csvBtag);
   fChain->SetBranchAddress("j10_e", &j10_e, &b_j10_e);
   fChain->SetBranchAddress("j10_pt", &j10_pt, &b_j10_pt);
   fChain->SetBranchAddress("j10_phi", &j10_phi, &b_j10_phi);
   fChain->SetBranchAddress("j10_eta", &j10_eta, &b_j10_eta);
   fChain->SetBranchAddress("j10_beta", &j10_beta, &b_j10_beta);
   fChain->SetBranchAddress("j10_betaStar", &j10_betaStar, &b_j10_betaStar);
   fChain->SetBranchAddress("j10_betaStarClassic", &j10_betaStarClassic, &b_j10_betaStarClassic);
   fChain->SetBranchAddress("j10_dR2Mean", &j10_dR2Mean, &b_j10_dR2Mean);
   fChain->SetBranchAddress("j10_algoPF1_csvBtag", &j10_algoPF1_csvBtag, &b_j10_algoPF1_csvBtag);
   fChain->SetBranchAddress("JetsMass", &JetsMass, &b_JetsMass);
   fChain->SetBranchAddress("dijet_E", &dijet_E, &b_dijet_E);
   fChain->SetBranchAddress("dijet_Pt", &dijet_Pt, &b_dijet_Pt);
   fChain->SetBranchAddress("dijet_Eta", &dijet_Eta, &b_dijet_Eta);
   fChain->SetBranchAddress("dijet_Phi", &dijet_Phi, &b_dijet_Phi);
   fChain->SetBranchAddress("nelectrons", &nelectrons, &b_nelectrons);
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("chargeEle", &chargeEle, &b_chargeEle);
   fChain->SetBranchAddress("chargeMu", &chargeMu, &b_chargeMu);
   fChain->SetBranchAddress("ptEle", &ptEle, &b_ptEle);
   fChain->SetBranchAddress("etaEle", &etaEle, &b_etaEle);
   fChain->SetBranchAddress("phiEle", &phiEle, &b_phiEle);
   fChain->SetBranchAddress("ptMu", &ptMu, &b_ptMu);
   fChain->SetBranchAddress("etaMu", &etaMu, &b_etaMu);
   fChain->SetBranchAddress("phiMu", &phiMu, &b_phiMu);
   fChain->SetBranchAddress("met_pfmet", &met_pfmet, &b_met_pfmet);
   fChain->SetBranchAddress("met_phi_pfmet", &met_phi_pfmet, &b_met_phi_pfmet);
   Notify();
}

Bool_t translateHGG::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void translateHGG::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t translateHGG::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef translateHGG_cxx
