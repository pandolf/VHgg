#ifndef RedNtpFinalizer_VHgg_cxx
#define RedNtpFinalizer_VHgg_cxx

#include "RedNtpFinalizer_VHgg.h"
#include <TH2.h>





RedNtpFinalizer_VHgg::RedNtpFinalizer_VHgg( const std::string& dataset, const std::string& selectionType, const std::string& bTaggerType ) : RedNtpFinalizer( "VHgg", dataset ) {
{

  bTaggerType_ = bTaggerType;

  setSelectionType(selectionType);

}

RedNtpFinalizer::~RedNtpFinalizer()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}





void RedNtpFinalizer_VHgg::finalize()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      bool isMC = ( run<5 );

      
      float eventWeight = 1.;


      if( isMC ) {

        if( nGenEvents_<=0 ) {
          std::cout << std::endl << "-> There must be a problem: nGenEvents is not positive!! Exiting" << std::endl;
          exit(99);
        }

        eventWeight = xSection_ / nGenEvents_ ;
 
        // pu/pt reeventWeighting
        if(dopureeventWeight) eventWeight *= pu_eventWeight; 
        if(doptreeventWeight) eventWeight *= pt_eventWeight; 

      }



      // analysis cuts:

      if(npu>=30) continue;

      // if(massggnewvtx<90 || massggnewvtx>190) continue;
      if(massggnewvtx<100 || massggnewvtx>180) continue;

      if((TMath::Abs(etascphot1)>1.4442&&TMath::Abs(etascphot1)<1.566)||(TMath::Abs(etascphot2)>1.4442&&TMath::Abs(etascphot2)<1.566)
         || TMath::Abs(etascphot1)>2.5 || TMath::Abs(etascphot2)>2.5) continue;  // acceptance

      //     if(ptphot1<ptphot1cut) continue; //pt first photon
      if(ptphot2<ptphot2cut) continue; //pt second photon


      if(ptphot1<ptphot1cut* massggnewvtx/120.) continue; //pt first photon
//       if(ptphot2<ptphot2cut* massggnewvtx/120.) continue; //pt second photon

      if(pthiggsmincut>0 && ptgg<pthiggsmincut) continue; //pt higgs min
      if(pthiggsmaxcut>0 && ptgg>=pthiggsmaxcut) continue; //pt higgs max


      if(ptjet1cut>0 && (ptcorrjet1<ptjet1cut || TMath::Abs(etajet1)>4.7)) continue; //pt first jet
      if(ptjet2cut>0 && (ptcorrjet2<ptjet2cut || TMath::Abs(etajet2)>4.7)) continue; //pt second jet


       // photon identification

       // r9 categories
       bool isr9phot1(0), isr9phot2(0);
 
       if(TMath::Abs(etascphot1)<1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)<1.4442 && r9phot2>.94) isr9phot2 = 1;
       if(TMath::Abs(etascphot1)>1.4442 && r9phot1>.94) isr9phot1 = 1;
       if(TMath::Abs(etascphot2)>1.4442 && r9phot2>.94) isr9phot2 = 1;
 
       if(r9cat == 1) {
         if(!isr9phot1 || !isr9phot2) continue;
       } else if (r9cat == 0){
         if(isr9phot1 && isr9phot2) continue;
       } 
 
       // photon id
       bool idphot1(0), idphot2(0), looseidphot1(0), looseidphot2(0), pxlphot1(1), pxlphot2(1);
 
 //       if(pixelseedcut) { 
 // 	pxlphot1 = !pid_haspixelseedphot1;
 // 	pxlphot2 = !pid_haspixelseedphot2;
 //       }
       
       idphot1 = (idcicphot1 >= cicselection);
       idphot2 = (idcicphot2 >= cicselection);
 
       if(!cs){ // photon id no control sample
 
       if(cicselection>0) {
         if(!(idphot1)) continue;
         if(!(idphot2)) continue;
       }else{
         if(!(idphot1 && pxlphot1)) continue;
         if(!(idphot2 && pxlphot2)) continue;
       }
       
       }else{ // photon id for control sample
        
         looseidphot1 = (idcicphot1 > 0 );
         looseidphot2 = (idcicphot2 > 0 );
         //	  if( !( (idphot1 && looseidphot2 && !idphot2) || (idphot2 && looseidphot1 && !idphot1) ) ) continue; 
         // Not perfect should be using the same electronVeto wrt CiC selection (now using matchedPromptEle veto)
         if( !( (idphot1 && !idphot2 && !pid_hasMatchedPromptElephot2) || (idphot2 && !idphot1 && !pid_hasMatchedPromptElephot1) ) ) continue; 
 
       }





   } //for entries



} //finalize


















Int_t RedNtpFinalizer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RedNtpFinalizer::LoadTree(Long64_t entry)
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

void RedNtpFinalizer::Init(TTree *tree)
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

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("rhoPF", &rhoPF, &b_rhoPF);
   fChain->SetBranchAddress("massgg", &massgg, &b_massgg);
   fChain->SetBranchAddress("ptgg", &ptgg, &b_ptgg);
   fChain->SetBranchAddress("ptggnewvtx", &ptggnewvtx, &b_ptggnewvtx);
   fChain->SetBranchAddress("phigg", &phigg, &b_phigg);
   fChain->SetBranchAddress("etagg", &etagg, &b_etagg);
   fChain->SetBranchAddress("massggnewvtx", &massggnewvtx, &b_massggnewvtx);
   fChain->SetBranchAddress("ptphot1", &ptphot1, &b_ptphot1);
   fChain->SetBranchAddress("ptphot2", &ptphot2, &b_ptphot2);
   fChain->SetBranchAddress("deltaRToTrackphot1", &deltaRToTrackphot1, &b_deltaRToTrackphot1);
   fChain->SetBranchAddress("deltaRToTrackphot2", &deltaRToTrackphot2, &b_deltaRToTrackphot2);
   fChain->SetBranchAddress("timephot1", &timephot1, &b_timephot1);
   fChain->SetBranchAddress("timephot2", &timephot2, &b_timephot2);
   fChain->SetBranchAddress("etaphot1", &etaphot1, &b_etaphot1);
   fChain->SetBranchAddress("etaphot2", &etaphot2, &b_etaphot2);
   fChain->SetBranchAddress("phiphot1", &phiphot1, &b_phiphot1);
   fChain->SetBranchAddress("phiphot2", &phiphot2, &b_phiphot2);
   fChain->SetBranchAddress("etascphot1", &etascphot1, &b_etascphot1);
   fChain->SetBranchAddress("etascphot2", &etascphot2, &b_etascphot2);
   fChain->SetBranchAddress("phiscphot1", &phiscphot1, &b_phiscphot1);
   fChain->SetBranchAddress("phiscphot2", &phiscphot2, &b_phiscphot2);
   fChain->SetBranchAddress("E1phot1", &E1phot1, &b_E1phot1);
   fChain->SetBranchAddress("E1phot2", &E1phot2, &b_E1phot2);
   fChain->SetBranchAddress("E9phot1", &E9phot1, &b_E9phot1);
   fChain->SetBranchAddress("E9phot2", &E9phot2, &b_E9phot2);
   fChain->SetBranchAddress("energyErrphot1", &energyErrphot1, &b_energyErrphot1);
   fChain->SetBranchAddress("energyErrphot2", &energyErrphot2, &b_energyErrphot2);
   fChain->SetBranchAddress("energySmearingphot1", &energySmearingphot1, &b_energySmearingphot1);
   fChain->SetBranchAddress("energySmearingphot2", &energySmearingphot2, &b_energySmearingphot2);
   fChain->SetBranchAddress("r9phot1", &r9phot1, &b_r9phot1);
   fChain->SetBranchAddress("r9phot2", &r9phot2, &b_r9phot2);
   fChain->SetBranchAddress("isemEGphot1", &isemEGphot1, &b_isemEGphot1);
   fChain->SetBranchAddress("isemEGphot2", &isemEGphot2, &b_isemEGphot2);
   fChain->SetBranchAddress("promptGamma", &promptGamma, &b_promptGamma);
   fChain->SetBranchAddress("LOGamma", &LOGamma, &b_LOGamma);
   fChain->SetBranchAddress("ISRGamma", &ISRGamma, &b_ISRGamma);
   fChain->SetBranchAddress("FSRGamma", &FSRGamma, &b_FSRGamma);
   fChain->SetBranchAddress("idmvaphot1", &idmvaphot1, &b_idmvaphot1);
   fChain->SetBranchAddress("idmvaphot2", &idmvaphot2, &b_idmvaphot2);
   fChain->SetBranchAddress("idcicphot1", &idcicphot1, &b_idcicphot1);
   fChain->SetBranchAddress("idcicphot2", &idcicphot2, &b_idcicphot2);
   fChain->SetBranchAddress("idcicnoelvetophot1", &idcicnoelvetophot1, &b_idcicnoelvetophot1);
   fChain->SetBranchAddress("idcicnoelvetophot2", &idcicnoelvetophot2, &b_idcicnoelvetophot2);
   fChain->SetBranchAddress("idcicpfphot1", &idcicpfphot1, &b_idcicpfphot1);
   fChain->SetBranchAddress("idcicpfphot2", &idcicpfphot2, &b_idcicpfphot2);
   fChain->SetBranchAddress("idcicpfnoelvetophot1", &idcicpfnoelvetophot1, &b_idcicpfnoelvetophot1);
   fChain->SetBranchAddress("idcicpfnoelvetophot2", &idcicpfnoelvetophot2, &b_idcicpfnoelvetophot2);
   fChain->SetBranchAddress("idelephot1", &idelephot1, &b_idelephot1);
   fChain->SetBranchAddress("idelephot2", &idelephot2, &b_idelephot2);
   fChain->SetBranchAddress("pid_isEMphot1", &pid_isEMphot1, &b_pid_isEMphot1);
   fChain->SetBranchAddress("pid_isEMphot2", &pid_isEMphot2, &b_pid_isEMphot2);
   fChain->SetBranchAddress("pid_haspixelseedphot1", &pid_haspixelseedphot1, &b_pid_haspixelseedphot1);
   fChain->SetBranchAddress("pid_haspixelseedphot2", &pid_haspixelseedphot2, &b_pid_haspixelseedphot2);
   fChain->SetBranchAddress("pid_jurECALphot1", &pid_jurECALphot1, &b_pid_jurECALphot1);
   fChain->SetBranchAddress("pid_jurECALphot2", &pid_jurECALphot2, &b_pid_jurECALphot2);
   fChain->SetBranchAddress("pid_twrHCALphot1", &pid_twrHCALphot1, &b_pid_twrHCALphot1);
   fChain->SetBranchAddress("pid_twrHCALphot2", &pid_twrHCALphot2, &b_pid_twrHCALphot2);
   fChain->SetBranchAddress("pid_HoverEphot1", &pid_HoverEphot1, &b_pid_HoverEphot1);
   fChain->SetBranchAddress("pid_HoverEphot2", &pid_HoverEphot2, &b_pid_HoverEphot2);
   fChain->SetBranchAddress("pid_hlwTrackphot1", &pid_hlwTrackphot1, &b_pid_hlwTrackphot1);
   fChain->SetBranchAddress("pid_hlwTrackphot2", &pid_hlwTrackphot2, &b_pid_hlwTrackphot2);
   fChain->SetBranchAddress("pid_etawidphot1", &pid_etawidphot1, &b_pid_etawidphot1);
   fChain->SetBranchAddress("pid_etawidphot2", &pid_etawidphot2, &b_pid_etawidphot2);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot1", &pid_hlwTrackNoDzphot1, &b_pid_hlwTrackNoDzphot1);
   fChain->SetBranchAddress("pid_hlwTrackNoDzphot2", &pid_hlwTrackNoDzphot2, &b_pid_hlwTrackNoDzphot2);
   fChain->SetBranchAddress("pid_hasMatchedConvphot1", &pid_hasMatchedConvphot1, &b_pid_hasMatchedConvphot1);
   fChain->SetBranchAddress("pid_hasMatchedConvphot2", &pid_hasMatchedConvphot2, &b_pid_hasMatchedConvphot2);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot1", &pid_hasMatchedPromptElephot1, &b_pid_hasMatchedPromptElephot1);
   fChain->SetBranchAddress("pid_hasMatchedPromptElephot2", &pid_hasMatchedPromptElephot2, &b_pid_hasMatchedPromptElephot2);
   fChain->SetBranchAddress("pid_sminphot1", &pid_sminphot1, &b_pid_sminphot1);
   fChain->SetBranchAddress("pid_sminphot2", &pid_sminphot2, &b_pid_sminphot2);
   fChain->SetBranchAddress("pid_smajphot1", &pid_smajphot1, &b_pid_smajphot1);
   fChain->SetBranchAddress("pid_smajphot2", &pid_smajphot2, &b_pid_smajphot2);
   fChain->SetBranchAddress("pid_ntrkphot1", &pid_ntrkphot1, &b_pid_ntrkphot1);
   fChain->SetBranchAddress("pid_ntrkphot2", &pid_ntrkphot2, &b_pid_ntrkphot2);
   fChain->SetBranchAddress("pid_ptisophot1", &pid_ptisophot1, &b_pid_ptisophot1);
   fChain->SetBranchAddress("pid_ptisophot2", &pid_ptisophot2, &b_pid_ptisophot2);
   fChain->SetBranchAddress("pid_ntrkcsphot1", &pid_ntrkcsphot1, &b_pid_ntrkcsphot1);
   fChain->SetBranchAddress("pid_ntrkcsphot2", &pid_ntrkcsphot2, &b_pid_ntrkcsphot2);
   fChain->SetBranchAddress("pid_ptisocsphot1", &pid_ptisocsphot1, &b_pid_ptisocsphot1);
   fChain->SetBranchAddress("pid_ptisocsphot2", &pid_ptisocsphot2, &b_pid_ptisocsphot2);
   fChain->SetBranchAddress("pid_ecalisophot1", &pid_ecalisophot1, &b_pid_ecalisophot1);
   fChain->SetBranchAddress("pid_ecalisophot2", &pid_ecalisophot2, &b_pid_ecalisophot2);
   fChain->SetBranchAddress("pid_hcalisophot1", &pid_hcalisophot1, &b_pid_hcalisophot1);
   fChain->SetBranchAddress("pid_hcalisophot2", &pid_hcalisophot2, &b_pid_hcalisophot2);
   fChain->SetBranchAddress("ptjet1", &ptjet1, &b_ptjet1);
   fChain->SetBranchAddress("ptjet2", &ptjet2, &b_ptjet2);
   fChain->SetBranchAddress("ptjet3", &ptjet3, &b_ptjet3);
   fChain->SetBranchAddress("ptjet4", &ptjet4, &b_ptjet4);
   fChain->SetBranchAddress("ptcorrjet1", &ptcorrjet1, &b_ptcorrjet1);
   fChain->SetBranchAddress("ptcorrjet2", &ptcorrjet2, &b_ptcorrjet2);
   fChain->SetBranchAddress("ptcorrjet3", &ptcorrjet3, &b_ptcorrjet3);
   fChain->SetBranchAddress("ptcorrjet4", &ptcorrjet4, &b_ptcorrjet4);
   fChain->SetBranchAddress("ecorrjet1", &ecorrjet1, &b_ecorrjet1);
   fChain->SetBranchAddress("ecorrjet2", &ecorrjet2, &b_ecorrjet2);
   fChain->SetBranchAddress("ecorrjet3", &ecorrjet3, &b_ecorrjet3);
   fChain->SetBranchAddress("ecorrjet4", &ecorrjet4, &b_ecorrjet4);
   fChain->SetBranchAddress("etajet1", &etajet1, &b_etajet1);
   fChain->SetBranchAddress("etajet2", &etajet2, &b_etajet2);
   fChain->SetBranchAddress("etajet3", &etajet3, &b_etajet3);
   fChain->SetBranchAddress("etajet4", &etajet4, &b_etajet4);
   fChain->SetBranchAddress("phijet1", &phijet1, &b_phijet1);
   fChain->SetBranchAddress("phijet2", &phijet2, &b_phijet2);
   fChain->SetBranchAddress("phijet3", &phijet3, &b_phijet3);
   fChain->SetBranchAddress("phijet4", &phijet4, &b_phijet4);
   fChain->SetBranchAddress("betajet1", &betajet1, &b_betajet1);
   fChain->SetBranchAddress("betajet2", &betajet2, &b_betajet2);
   fChain->SetBranchAddress("betastarjet1", &betastarjet1, &b_betastarjet1);
   fChain->SetBranchAddress("betastarjet2", &betastarjet2, &b_betastarjet2);
   fChain->SetBranchAddress("btagvtxjet1", &btagvtxjet1, &b_btagvtxjet1);
   fChain->SetBranchAddress("btagvtxjet2", &btagvtxjet2, &b_btagvtxjet2);
   fChain->SetBranchAddress("btagvtxjet3", &btagvtxjet3, &b_btagvtxjet3);
   fChain->SetBranchAddress("btagvtxjet4", &btagvtxjet4, &b_btagvtxjet4);
   fChain->SetBranchAddress("btagtrkjet1", &btagtrkjet1, &b_btagtrkjet1);
   fChain->SetBranchAddress("btagtrkjet2", &btagtrkjet2, &b_btagtrkjet2);
   fChain->SetBranchAddress("btagtrkjet3", &btagtrkjet3, &b_btagtrkjet3);
   fChain->SetBranchAddress("btagtrkjet4", &btagtrkjet4, &b_btagtrkjet4);
   fChain->SetBranchAddress("btagjprobjet1", &btagjprobjet1, &b_btagjprobjet1);
   fChain->SetBranchAddress("btagjprobjet2", &btagjprobjet2, &b_btagjprobjet2);
   fChain->SetBranchAddress("btagjprobjet3", &btagjprobjet3, &b_btagjprobjet3);
   fChain->SetBranchAddress("btagjprobjet4", &btagjprobjet4, &b_btagjprobjet4);
   fChain->SetBranchAddress("ptDjet1", &ptDjet1, &b_ptDjet1);
   fChain->SetBranchAddress("rmsjet1", &rmsjet1, &b_rmsjet1);
   fChain->SetBranchAddress("ntrkjet1", &ntrkjet1, &b_ntrkjet1);
   fChain->SetBranchAddress("nneutjet1", &nneutjet1, &b_nneutjet1);
   fChain->SetBranchAddress("jetIdSimple_mvajet1", &jetIdSimple_mvajet1, &b_jetIdSimple_mvajet1);
   fChain->SetBranchAddress("jetIdFull_mvajet1", &jetIdFull_mvajet1, &b_jetIdFull_mvajet1);
   fChain->SetBranchAddress("jetId_dR2Meanjet1", &jetId_dR2Meanjet1, &b_jetId_dR2Meanjet1);
   fChain->SetBranchAddress("jetId_betaStarClassicjet1", &jetId_betaStarClassicjet1, &b_jetId_betaStarClassicjet1);
   fChain->SetBranchAddress("jetId_frac01jet1", &jetId_frac01jet1, &b_jetId_frac01jet1);
   fChain->SetBranchAddress("jetId_frac02jet1", &jetId_frac02jet1, &b_jetId_frac02jet1);
   fChain->SetBranchAddress("jetId_frac03jet1", &jetId_frac03jet1, &b_jetId_frac03jet1);
   fChain->SetBranchAddress("jetId_frac04jet1", &jetId_frac04jet1, &b_jetId_frac04jet1);
   fChain->SetBranchAddress("jetId_frac05jet1", &jetId_frac05jet1, &b_jetId_frac05jet1);
   fChain->SetBranchAddress("jetId_betajet1", &jetId_betajet1, &b_jetId_betajet1);
   fChain->SetBranchAddress("jetId_betaStarjet1", &jetId_betaStarjet1, &b_jetId_betaStarjet1);
   fChain->SetBranchAddress("jetIdCutBased_wpjet1", &jetIdCutBased_wpjet1, &b_jetIdCutBased_wpjet1);
   fChain->SetBranchAddress("jetIdSimple_wpjet1", &jetIdSimple_wpjet1, &b_jetIdSimple_wpjet1);
   fChain->SetBranchAddress("jetIdFull_wpjet1", &jetIdFull_wpjet1, &b_jetIdFull_wpjet1);
   fChain->SetBranchAddress("ptDjet2", &ptDjet2, &b_ptDjet2);
   fChain->SetBranchAddress("rmsjet2", &rmsjet2, &b_rmsjet2);
   fChain->SetBranchAddress("ntrkjet2", &ntrkjet2, &b_ntrkjet2);
   fChain->SetBranchAddress("nneutjet2", &nneutjet2, &b_nneutjet2);
   fChain->SetBranchAddress("jetIdSimple_mvajet2", &jetIdSimple_mvajet2, &b_jetIdSimple_mvajet2);
   fChain->SetBranchAddress("jetIdFull_mvajet2", &jetIdFull_mvajet2, &b_jetIdFull_mvajet2);
   fChain->SetBranchAddress("jetId_dR2Meanjet2", &jetId_dR2Meanjet2, &b_jetId_dR2Meanjet2);
   fChain->SetBranchAddress("jetId_betaStarClassicjet2", &jetId_betaStarClassicjet2, &b_jetId_betaStarClassicjet2);
   fChain->SetBranchAddress("jetIdCutBased_wpjet2", &jetIdCutBased_wpjet2, &b_jetIdCutBased_wpjet2);
   fChain->SetBranchAddress("jetIdSimple_wpjet2", &jetIdSimple_wpjet2, &b_jetIdSimple_wpjet2);
   fChain->SetBranchAddress("jetIdFull_wpjet2", &jetIdFull_wpjet2, &b_jetIdFull_wpjet2);
   fChain->SetBranchAddress("jetId_frac01jet2", &jetId_frac01jet2, &b_jetId_frac01jet2);
   fChain->SetBranchAddress("jetId_frac02jet2", &jetId_frac02jet2, &b_jetId_frac02jet2);
   fChain->SetBranchAddress("jetId_frac03jet2", &jetId_frac03jet2, &b_jetId_frac03jet2);
   fChain->SetBranchAddress("jetId_frac04jet2", &jetId_frac04jet2, &b_jetId_frac04jet2);
   fChain->SetBranchAddress("jetId_frac05jet2", &jetId_frac05jet2, &b_jetId_frac05jet2);
   fChain->SetBranchAddress("jetId_betajet2", &jetId_betajet2, &b_jetId_betajet2);
   fChain->SetBranchAddress("jetId_betaStarjet2", &jetId_betaStarjet2, &b_jetId_betaStarjet2);
   fChain->SetBranchAddress("ptDjet3", &ptDjet3, &b_ptDjet3);
   fChain->SetBranchAddress("rmsjet3", &rmsjet3, &b_rmsjet3);
   fChain->SetBranchAddress("ntrkjet3", &ntrkjet3, &b_ntrkjet3);
   fChain->SetBranchAddress("nneutjet3", &nneutjet3, &b_nneutjet3);
   fChain->SetBranchAddress("ptDjet4", &ptDjet4, &b_ptDjet4);
   fChain->SetBranchAddress("rmsjet4", &rmsjet4, &b_rmsjet4);
   fChain->SetBranchAddress("ntrkjet4", &ntrkjet4, &b_ntrkjet4);
   fChain->SetBranchAddress("nneutjet4", &nneutjet4, &b_nneutjet4);
   fChain->SetBranchAddress("assjet1", &assjet1, &b_assjet1);
   fChain->SetBranchAddress("assjet2", &assjet2, &b_assjet2);
   fChain->SetBranchAddress("deltaeta", &deltaeta, &b_deltaeta);
   fChain->SetBranchAddress("zeppenjet", &zeppenjet, &b_zeppenjet);
   fChain->SetBranchAddress("deltaphi", &deltaphi, &b_deltaphi);
   fChain->SetBranchAddress("deltaphinewvtx", &deltaphinewvtx, &b_deltaphinewvtx);
   fChain->SetBranchAddress("deltaphigg", &deltaphigg, &b_deltaphigg);
   fChain->SetBranchAddress("invmassjet", &invmassjet, &b_invmassjet);
   fChain->SetBranchAddress("invmass2g1j", &invmass2g1j, &b_invmass2g1j);
   fChain->SetBranchAddress("invmass2g2j", &invmass2g2j, &b_invmass2g2j);
   fChain->SetBranchAddress("pt2g2j", &pt2g2j, &b_pt2g2j);
   fChain->SetBranchAddress("eta2j", &eta2j, &b_eta2j);
   fChain->SetBranchAddress("phi2j", &phi2j, &b_phi2j);
   fChain->SetBranchAddress("pt2j", &pt2j, &b_pt2j);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("vtxId", &vtxId, &b_vtxId);
   fChain->SetBranchAddress("vtxPos_x", &vtxPos_x, &b_vtxPos_x);
   fChain->SetBranchAddress("vtxPos_y", &vtxPos_y, &b_vtxPos_y);
   fChain->SetBranchAddress("vtxPos_z", &vtxPos_z, &b_vtxPos_z);
   fChain->SetBranchAddress("vtxIdMVA", &vtxIdMVA, &b_vtxIdMVA);
   fChain->SetBranchAddress("vtxIdEvtProb", &vtxIdEvtProb, &b_vtxIdEvtProb);
   fChain->SetBranchAddress("sMet", &sMet, &b_sMet);
   fChain->SetBranchAddress("eMet", &eMet, &b_eMet);
   fChain->SetBranchAddress("phiMet", &phiMet, &b_phiMet);
   fChain->SetBranchAddress("signifMet", &signifMet, &b_signifMet);
   fChain->SetBranchAddress("eSmearedMet", &eSmearedMet, &b_eSmearedMet);
   fChain->SetBranchAddress("phiSmearedMet", &phiSmearedMet, &b_phiSmearedMet);
   fChain->SetBranchAddress("eShiftedMet", &eShiftedMet, &b_eShiftedMet);
   fChain->SetBranchAddress("phiShiftedMet", &phiShiftedMet, &b_phiShiftedMet);
   fChain->SetBranchAddress("eShiftedScaledMet", &eShiftedScaledMet, &b_eShiftedScaledMet);
   fChain->SetBranchAddress("phiShiftedScaledMet", &phiShiftedScaledMet, &b_phiShiftedScaledMet);
   fChain->SetBranchAddress("eSmearedShiftedMet", &eSmearedShiftedMet, &b_eSmearedShiftedMet);
   fChain->SetBranchAddress("phiSmearedShiftedMet", &phiSmearedShiftedMet, &b_phiSmearedShiftedMet);
   fChain->SetBranchAddress("eShiftedScaledMetPUcorr", &eShiftedScaledMetPUcorr, &b_eShiftedScaledMetPUcorr);
   fChain->SetBranchAddress("phiShiftedScaledMetPUcorr", &phiShiftedScaledMetPUcorr, &b_phiShiftedScaledMetPUcorr);
   fChain->SetBranchAddress("eSmearedShiftedMePUcorrt", &eSmearedShiftedMePUcorrt, &b_eSmearedShiftedMetPUcorr);
   fChain->SetBranchAddress("phiSmearedShiftedMetPUcorr", &phiSmearedShiftedMetPUcorr, &b_phiSmearedShiftedMetPUcorr);
   fChain->SetBranchAddress("sCorrMet", &sCorrMet, &b_sCorrMet);
   fChain->SetBranchAddress("eCorrMet", &eCorrMet, &b_eCorrMet);
   fChain->SetBranchAddress("phiCorrMet", &phiCorrMet, &b_phiCorrMet);
   fChain->SetBranchAddress("signifCorrMet", &signifCorrMet, &b_signifCorrMet);
   fChain->SetBranchAddress("smuCorrMet", &smuCorrMet, &b_smuCorrMet);
   fChain->SetBranchAddress("emuCorrMet", &emuCorrMet, &b_emuCorrMet);
   fChain->SetBranchAddress("phimuCorrMet", &phimuCorrMet, &b_phimuCorrMet);
   fChain->SetBranchAddress("signifmuCorrMet", &signifmuCorrMet, &b_signifmuCorrMet);
   fChain->SetBranchAddress("sNoHFMet", &sNoHFMet, &b_sNoHFMet);
   fChain->SetBranchAddress("eNoHFMet", &eNoHFMet, &b_eNoHFMet);
   fChain->SetBranchAddress("phiNoHFMet", &phiNoHFMet, &b_phiNoHFMet);
   fChain->SetBranchAddress("signifNoHFMet", &signifNoHFMet, &b_signifNoHFMet);
   fChain->SetBranchAddress("stcMet", &stcMet, &b_stcMet);
   fChain->SetBranchAddress("etcMet", &etcMet, &b_etcMet);
   fChain->SetBranchAddress("phitcMet", &phitcMet, &b_phitcMet);
   fChain->SetBranchAddress("signiftcMet", &signiftcMet, &b_signiftcMet);
   fChain->SetBranchAddress("sglobalPfMet", &sglobalPfMet, &b_sglobalPfMet);
   fChain->SetBranchAddress("eglobalPfMet", &eglobalPfMet, &b_eglobalPfMet);
   fChain->SetBranchAddress("phiglobalPfMet", &phiglobalPfMet, &b_phiglobalPfMet);
   fChain->SetBranchAddress("signifglobalPfMet", &signifglobalPfMet, &b_signifglobalPfMet);
   fChain->SetBranchAddress("scentralPfMet", &scentralPfMet, &b_scentralPfMet);
   fChain->SetBranchAddress("ecentralPfMet", &ecentralPfMet, &b_ecentralPfMet);
   fChain->SetBranchAddress("phicentralPfMet", &phicentralPfMet, &b_phicentralPfMet);
   fChain->SetBranchAddress("signifcentralPfMet", &signifcentralPfMet, &b_signifcentralPfMet);
   fChain->SetBranchAddress("eassocPfMet", &eassocPfMet, &b_eassocPfMet);
   fChain->SetBranchAddress("phiassocPfMet", &phiassocPfMet, &b_phiassocPfMet);
   fChain->SetBranchAddress("signifassocPfMet", &signifassocPfMet, &b_signifassocPfMet);
   fChain->SetBranchAddress("eassocOtherVtxPfMet", &eassocOtherVtxPfMet, &b_eassocOtherVtxPfMet);
   fChain->SetBranchAddress("phiassocOtherVtxPfMet", &phiassocOtherVtxPfMet, &b_phiassocOtherVtxPfMet);
   fChain->SetBranchAddress("signifassocOtherVtxPfMet", &signifassocOtherVtxPfMet, &b_signifassocOtherVtxPfMet);
   fChain->SetBranchAddress("etrkPfMet", &etrkPfMet, &b_etrkPfMet);
   fChain->SetBranchAddress("phitrkPfMet", &phitrkPfMet, &b_phitrkPfMet);
   fChain->SetBranchAddress("signiftrkPfMet", &signiftrkPfMet, &b_signiftrkPfMet);
   fChain->SetBranchAddress("ecleanPfMet", &ecleanPfMet, &b_ecleanPfMet);
   fChain->SetBranchAddress("phicleanPfMet", &phicleanPfMet, &b_phicleanPfMet);
   fChain->SetBranchAddress("signifcleanPfMet", &signifcleanPfMet, &b_signifcleanPfMet);
   fChain->SetBranchAddress("ecleanedSaclayPfMet", &ecleanedSaclayPfMet, &b_ecleanedSaclayPfMet);
   fChain->SetBranchAddress("phicleanedSaclayPfMet", &phicleanedSaclayPfMet, &b_phicleanedSaclayPfMet);
   fChain->SetBranchAddress("signifcleanedSaclayPfMet", &signifcleanedSaclayPfMet, &b_signifcleanedSaclayPfMet);
   fChain->SetBranchAddress("eminTypeICleanSaclayPfMet", &eminTypeICleanSaclayPfMet, &b_eminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("phiminTypeICleanSaclayPfMet", &phiminTypeICleanSaclayPfMet, &b_phiminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("signifminTypeICleanSaclayPfMet", &signifminTypeICleanSaclayPfMet, &b_signifminTypeICleanSaclayPfMet);
   fChain->SetBranchAddress("globalPfSums", &globalPfSums, &b_globalPfSums);
   fChain->SetBranchAddress("spfMet", &spfMet, &b_spfMet);
   fChain->SetBranchAddress("epfMet", &epfMet, &b_epfMet);
   fChain->SetBranchAddress("phipfMet", &phipfMet, &b_phipfMet);
   fChain->SetBranchAddress("signifpfMet", &signifpfMet, &b_signifpfMet);
   fChain->SetBranchAddress("spfMetType1", &spfMetType1, &b_spfMetType1);
   fChain->SetBranchAddress("epfMetType1", &epfMetType1, &b_epfMetType1);
   fChain->SetBranchAddress("phipfMetType1", &phipfMetType1, &b_phipfMetType1);
   fChain->SetBranchAddress("signifpfMetType1", &signifpfMetType1, &b_signifpfMetType1);
   fChain->SetBranchAddress("sMetGen", &sMetGen, &b_sMetGen);
   fChain->SetBranchAddress("eMetGen", &eMetGen, &b_eMetGen);
   fChain->SetBranchAddress("phiMetGen", &phiMetGen, &b_phiMetGen);
   fChain->SetBranchAddress("signifMetGen", &signifMetGen, &b_signifMetGen);
   fChain->SetBranchAddress("sMetGen2", &sMetGen2, &b_sMetGen2);
   fChain->SetBranchAddress("eMetGen2", &eMetGen2, &b_eMetGen2);
   fChain->SetBranchAddress("phiMetGen2", &phiMetGen2, &b_phiMetGen2);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("NtotEvents", &NtotEvents, &b_NtotEvents);
   fChain->SetBranchAddress("xsection", &xsection, &b_xsection);
   fChain->SetBranchAddress("EquivLumi", &EquivLumi, &b_EquivLumi);
   fChain->SetBranchAddress("SampleID", &SampleID, &b_SampleID);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("pt_weight", &pt_weight, &b_pt_weight);
   fChain->SetBranchAddress("gen_custom_processId", &gen_custom_processId, &b_gen_custom_processId);
   fChain->SetBranchAddress("gen_pt_gamma1", &gen_pt_gamma1, &b_gen_pt_gamma1);
   fChain->SetBranchAddress("gen_pt_gamma2", &gen_pt_gamma2, &b_gen_pt_gamma2);
   fChain->SetBranchAddress("gen_eta_gamma1", &gen_eta_gamma1, &b_gen_eta_gamma1);
   fChain->SetBranchAddress("gen_eta_gamma2", &gen_eta_gamma2, &b_gen_eta_gamma2);
   fChain->SetBranchAddress("gen_phi_gamma1", &gen_phi_gamma1, &b_gen_phi_gamma1);
   fChain->SetBranchAddress("gen_phi_gamma2", &gen_phi_gamma2, &b_gen_phi_gamma2);
   fChain->SetBranchAddress("gen_pt_genjet1", &gen_pt_genjet1, &b_gen_pt_genjet1);
   fChain->SetBranchAddress("gen_pt_genjet2", &gen_pt_genjet2, &b_gen_pt_genjet2);
   fChain->SetBranchAddress("gen_eta_genjet1", &gen_eta_genjet1, &b_gen_eta_genjet1);
   fChain->SetBranchAddress("gen_eta_genjet2", &gen_eta_genjet2, &b_gen_eta_genjet2);
   fChain->SetBranchAddress("gen_phi_genjet1", &gen_phi_genjet1, &b_gen_phi_genjet1);
   fChain->SetBranchAddress("gen_phi_genjet2", &gen_phi_genjet2, &b_gen_phi_genjet2);
   fChain->SetBranchAddress("gen_mass_diphoton", &gen_mass_diphoton, &b_gen_mass_diphoton);
   fChain->SetBranchAddress("gen_pt_diphoton", &gen_pt_diphoton, &b_gen_pt_diphoton);
   fChain->SetBranchAddress("gen_eta_diphoton", &gen_eta_diphoton, &b_gen_eta_diphoton);
   fChain->SetBranchAddress("gen_phi_diphoton", &gen_phi_diphoton, &b_gen_phi_diphoton);
   fChain->SetBranchAddress("gen_mass_dijet", &gen_mass_dijet, &b_gen_mass_dijet);
   fChain->SetBranchAddress("gen_pt_dijet", &gen_pt_dijet, &b_gen_pt_dijet);
   fChain->SetBranchAddress("gen_eta_dijet", &gen_eta_dijet, &b_gen_eta_dijet);
   fChain->SetBranchAddress("gen_phi_dijet", &gen_phi_dijet, &b_gen_phi_dijet);
   fChain->SetBranchAddress("gen_zeppenfeld", &gen_zeppenfeld, &b_gen_zeppenfeld);
   fChain->SetBranchAddress("gen_pt_lep1", &gen_pt_lep1, &b_gen_pt_lep1);
   fChain->SetBranchAddress("gen_pt_lep2", &gen_pt_lep2, &b_gen_pt_lep2);
   fChain->SetBranchAddress("gen_eta_lep1", &gen_eta_lep1, &b_gen_eta_lep1);
   fChain->SetBranchAddress("gen_eta_lep2", &gen_eta_lep2, &b_gen_eta_lep2);
   fChain->SetBranchAddress("gen_phi_lep1", &gen_phi_lep1, &b_gen_phi_lep1);
   fChain->SetBranchAddress("gen_phi_lep2", &gen_phi_lep2, &b_gen_phi_lep2);
   fChain->SetBranchAddress("gen_pid_lep1", &gen_pid_lep1, &b_gen_pid_lep1);
   fChain->SetBranchAddress("gen_pid_lep2", &gen_pid_lep2, &b_gen_pid_lep2);
   fChain->SetBranchAddress("ptele1", &ptele1, &b_ptele1);
   fChain->SetBranchAddress("ptele2", &ptele2, &b_ptele2);
   fChain->SetBranchAddress("etaele1", &etaele1, &b_etaele1);
   fChain->SetBranchAddress("etaele2", &etaele2, &b_etaele2);
   fChain->SetBranchAddress("phiele1", &phiele1, &b_phiele1);
   fChain->SetBranchAddress("phiele2", &phiele2, &b_phiele2);
   fChain->SetBranchAddress("eneele1", &eneele1, &b_eneele1);
   fChain->SetBranchAddress("eneele2", &eneele2, &b_eneele2);
   fChain->SetBranchAddress("sIeIeele1", &sIeIeele1, &b_sIeIeele1);
   fChain->SetBranchAddress("sIeIeele2", &sIeIeele2, &b_sIeIeele2);
   fChain->SetBranchAddress("dphiele1", &dphiele1, &b_dphiele1);
   fChain->SetBranchAddress("dphiele2", &dphiele2, &b_dphiele2);
   fChain->SetBranchAddress("detaele1", &detaele1, &b_detaele1);
   fChain->SetBranchAddress("detaele2", &detaele2, &b_detaele2);
   fChain->SetBranchAddress("mhitsele1", &mhitsele1, &b_mhitsele1);
   fChain->SetBranchAddress("mhitsele2", &mhitsele2, &b_mhitsele2);
   fChain->SetBranchAddress("d0ele1", &d0ele1, &b_d0ele1);
   fChain->SetBranchAddress("d0ele2", &d0ele2, &b_d0ele2);
   fChain->SetBranchAddress("dzele1", &dzele1, &b_dzele1);
   fChain->SetBranchAddress("dzele2", &dzele2, &b_dzele2);
   fChain->SetBranchAddress("invMassele1g1", &invMassele1g1, &b_invMassele1g1);
   fChain->SetBranchAddress("invMassele1g2", &invMassele1g2, &b_invMassele1g2);
   fChain->SetBranchAddress("invMassele2g1", &invMassele2g1, &b_invMassele2g1);
   fChain->SetBranchAddress("invMassele2g2", &invMassele2g2, &b_invMassele2g2);
   fChain->SetBranchAddress("oEmoPele1", &oEmoPele1, &b_oEmoPele1);
   fChain->SetBranchAddress("oEmoPele2", &oEmoPele2, &b_oEmoPele2);
   fChain->SetBranchAddress("mvanotrigele1", &mvanotrigele1, &b_mvanotrigele1);
   fChain->SetBranchAddress("mvanotrigele2", &mvanotrigele2, &b_mvanotrigele2);
   fChain->SetBranchAddress("mvatrigele1", &mvatrigele1, &b_mvatrigele1);
   fChain->SetBranchAddress("mvatrigele2", &mvatrigele2, &b_mvatrigele2);
   fChain->SetBranchAddress("matchconvele1", &matchconvele1, &b_matchconvele1);
   fChain->SetBranchAddress("matchconvele2", &matchconvele2, &b_matchconvele2);
   fChain->SetBranchAddress("chHadIso03ele1", &chHadIso03ele1, &b_chHadIso03ele1);
   fChain->SetBranchAddress("chHadIso03ele2", &chHadIso03ele2, &b_chHadIso03ele2);
   fChain->SetBranchAddress("nHadIso03ele1", &nHadIso03ele1, &b_nHadIso03ele1);
   fChain->SetBranchAddress("nHadIso03ele2", &nHadIso03ele2, &b_nHadIso03ele2);
   fChain->SetBranchAddress("photIso03ele1", &photIso03ele1, &b_photIso03ele1);
   fChain->SetBranchAddress("photIso03ele2", &photIso03ele2, &b_photIso03ele2);
   fChain->SetBranchAddress("pteleloose1", &pteleloose1, &b_pteleloose1);
   fChain->SetBranchAddress("pteleloose2", &pteleloose2, &b_pteleloose2);
   fChain->SetBranchAddress("etaeleloose1", &etaeleloose1, &b_etaeleloose1);
   fChain->SetBranchAddress("etaeleloose2", &etaeleloose2, &b_etaeleloose2);
   fChain->SetBranchAddress("phieleloose1", &phieleloose1, &b_phieleloose1);
   fChain->SetBranchAddress("phieleloose2", &phieleloose2, &b_phieleloose2);
   fChain->SetBranchAddress("eneeleloose1", &eneeleloose1, &b_eneeleloose1);
   fChain->SetBranchAddress("eneeleloose2", &eneeleloose2, &b_eneeleloose2);
   fChain->SetBranchAddress("sIeIeeleloose1", &sIeIeeleloose1, &b_sIeIeeleloose1);
   fChain->SetBranchAddress("sIeIeeleloose2", &sIeIeeleloose2, &b_sIeIeeleloose2);
   fChain->SetBranchAddress("dphieleloose1", &dphieleloose1, &b_dphieleloose1);
   fChain->SetBranchAddress("dphieleloose2", &dphieleloose2, &b_dphieleloose2);
   fChain->SetBranchAddress("detaeleloose1", &detaeleloose1, &b_detaeleloose1);
   fChain->SetBranchAddress("detaeleloose2", &detaeleloose2, &b_detaeleloose2);
   fChain->SetBranchAddress("mhitseleloose1", &mhitseleloose1, &b_mhitseleloose1);
   fChain->SetBranchAddress("mhitseleloose2", &mhitseleloose2, &b_mhitseleloose2);
   fChain->SetBranchAddress("d0eleloose1", &d0eleloose1, &b_d0eleloose1);
   fChain->SetBranchAddress("d0eleloose2", &d0eleloose2, &b_d0eleloose2);
   fChain->SetBranchAddress("dzeleloose1", &dzeleloose1, &b_dzeleloose1);
   fChain->SetBranchAddress("dzeleloose2", &dzeleloose2, &b_dzeleloose2);
   fChain->SetBranchAddress("invMasseleloose1g1", &invMasseleloose1g1, &b_invMasseleloose1g1);
   fChain->SetBranchAddress("invMasseleloose1g2", &invMasseleloose1g2, &b_invMasseleloose1g2);
   fChain->SetBranchAddress("invMasseleloose2g1", &invMasseleloose2g1, &b_invMasseleloose2g1);
   fChain->SetBranchAddress("invMasseleloose2g2", &invMasseleloose2g2, &b_invMasseleloose2g2);
   fChain->SetBranchAddress("oEmoPeleloose1", &oEmoPeleloose1, &b_oEmoPeleloose1);
   fChain->SetBranchAddress("oEmoPeleloose2", &oEmoPeleloose2, &b_oEmoPeleloose2);
   fChain->SetBranchAddress("mvanotrigeleloose1", &mvanotrigeleloose1, &b_mvanotrigeleloose1);
   fChain->SetBranchAddress("mvanotrigeleloose2", &mvanotrigeleloose2, &b_mvanotrigeleloose2);
   fChain->SetBranchAddress("mvatrigeleloose1", &mvatrigeleloose1, &b_mvatrigeleloose1);
   fChain->SetBranchAddress("mvatrigeleloose2", &mvatrigeleloose2, &b_mvatrigeleloose2);
   fChain->SetBranchAddress("matchconveleloose1", &matchconveleloose1, &b_matchconveleloose1);
   fChain->SetBranchAddress("matchconveleloose2", &matchconveleloose2, &b_matchconveleloose2);
   fChain->SetBranchAddress("chHadIso03eleloose1", &chHadIso03eleloose1, &b_chHadIso03eleloose1);
   fChain->SetBranchAddress("chHadIso03eleloose2", &chHadIso03eleloose2, &b_chHadIso03eleloose2);
   fChain->SetBranchAddress("nHadIso03eleloose1", &nHadIso03eleloose1, &b_nHadIso03eleloose1);
   fChain->SetBranchAddress("nHadIso03eleloose2", &nHadIso03eleloose2, &b_nHadIso03eleloose2);
   fChain->SetBranchAddress("photIso03eleloose1", &photIso03eleloose1, &b_photIso03eleloose1);
   fChain->SetBranchAddress("photIso03eleloose2", &photIso03eleloose2, &b_photIso03eleloose2);
   fChain->SetBranchAddress("ptmu1", &ptmu1, &b_ptmu1);
   fChain->SetBranchAddress("ptmu2", &ptmu2, &b_ptmu2);
   fChain->SetBranchAddress("etamu1", &etamu1, &b_etamu1);
   fChain->SetBranchAddress("etamu2", &etamu2, &b_etamu2);
   fChain->SetBranchAddress("phimu1", &phimu1, &b_phimu1);
   fChain->SetBranchAddress("phimu2", &phimu2, &b_phimu2);
   fChain->SetBranchAddress("enemu1", &enemu1, &b_enemu1);
   fChain->SetBranchAddress("enemu2", &enemu2, &b_enemu2);
   fChain->SetBranchAddress("pixhitsmu1", &pixhitsmu1, &b_pixhitsmu1);
   fChain->SetBranchAddress("pixhitsmu2", &pixhitsmu2, &b_pixhitsmu2);
   fChain->SetBranchAddress("trkhitsmu1", &trkhitsmu1, &b_trkhitsmu1);
   fChain->SetBranchAddress("trkhitsmu2", &trkhitsmu2, &b_trkhitsmu2);
   fChain->SetBranchAddress("hitsmu1", &hitsmu1, &b_hitsmu1);
   fChain->SetBranchAddress("hitsmu2", &hitsmu2, &b_hitsmu2);
   fChain->SetBranchAddress("chi2mu1", &chi2mu1, &b_chi2mu1);
   fChain->SetBranchAddress("chi2mu2", &chi2mu2, &b_chi2mu2);
   fChain->SetBranchAddress("matchmu1", &matchmu1, &b_matchmu1);
   fChain->SetBranchAddress("matchmu2", &matchmu2, &b_matchmu2);
   fChain->SetBranchAddress("d0mu1", &d0mu1, &b_d0mu1);
   fChain->SetBranchAddress("d0mu2", &d0mu2, &b_d0mu2);
   fChain->SetBranchAddress("dzmu1", &dzmu1, &b_dzmu1);
   fChain->SetBranchAddress("dzmu2", &dzmu2, &b_dzmu2);
   fChain->SetBranchAddress("chHadmu1", &chHadmu1, &b_chHadmu1);
   fChain->SetBranchAddress("chHadmu2", &chHadmu2, &b_chHadmu2);
   fChain->SetBranchAddress("nHadmu1", &nHadmu1, &b_nHadmu1);
   fChain->SetBranchAddress("nHadmu2", &nHadmu2, &b_nHadmu2);
   fChain->SetBranchAddress("photmu1", &photmu1, &b_photmu1);
   fChain->SetBranchAddress("photmu2", &photmu2, &b_photmu2);
   fChain->SetBranchAddress("puptmu1", &puptmu1, &b_puptmu1);
   fChain->SetBranchAddress("puptmu2", &puptmu2, &b_puptmu2);
   fChain->SetBranchAddress("ptmuloose1", &ptmuloose1, &b_ptmuloose1);
   fChain->SetBranchAddress("ptmuloose2", &ptmuloose2, &b_ptmuloose2);
   fChain->SetBranchAddress("etamuloose1", &etamuloose1, &b_etamuloose1);
   fChain->SetBranchAddress("etamuloose2", &etamuloose2, &b_etamuloose2);
   fChain->SetBranchAddress("phimuloose1", &phimuloose1, &b_phimuloose1);
   fChain->SetBranchAddress("phimuloose2", &phimuloose2, &b_phimuloose2);
   fChain->SetBranchAddress("enemuloose1", &enemuloose1, &b_enemuloose1);
   fChain->SetBranchAddress("enemuloose2", &enemuloose2, &b_enemuloose2);
   fChain->SetBranchAddress("pixhitsmuloose1", &pixhitsmuloose1, &b_pixhitsmuloose1);
   fChain->SetBranchAddress("pixhitsmuloose2", &pixhitsmuloose2, &b_pixhitsmuloose2);
   fChain->SetBranchAddress("trkhitsmuloose1", &trkhitsmuloose1, &b_trkhitsmuloose1);
   fChain->SetBranchAddress("trkhitsmuloose2", &trkhitsmuloose2, &b_trkhitsmuloose2);
   fChain->SetBranchAddress("hitsmuloose1", &hitsmuloose1, &b_hitsmuloose1);
   fChain->SetBranchAddress("hitsmuloose2", &hitsmuloose2, &b_hitsmuloose2);
   fChain->SetBranchAddress("chi2muloose1", &chi2muloose1, &b_chi2muloose1);
   fChain->SetBranchAddress("chi2muloose2", &chi2muloose2, &b_chi2muloose2);
   fChain->SetBranchAddress("matchmuloose1", &matchmuloose1, &b_matchmuloose1);
   fChain->SetBranchAddress("matchmuloose2", &matchmuloose2, &b_matchmuloose2);
   fChain->SetBranchAddress("d0muloose1", &d0muloose1, &b_d0muloose1);
   fChain->SetBranchAddress("d0muloose2", &d0muloose2, &b_d0muloose2);
   fChain->SetBranchAddress("dzmuloose1", &dzmuloose1, &b_dzmuloose1);
   fChain->SetBranchAddress("dzmuloose2", &dzmuloose2, &b_dzmuloose2);
   fChain->SetBranchAddress("chHadmuloose1", &chHadmuloose1, &b_chHadmuloose1);
   fChain->SetBranchAddress("chHadmuloose2", &chHadmuloose2, &b_chHadmuloose2);
   fChain->SetBranchAddress("nHadmuloose1", &nHadmuloose1, &b_nHadmuloose1);
   fChain->SetBranchAddress("nHadmuloose2", &nHadmuloose2, &b_nHadmuloose2);
   fChain->SetBranchAddress("photmuloose1", &photmuloose1, &b_photmuloose1);
   fChain->SetBranchAddress("photmuloose2", &photmuloose2, &b_photmuloose2);
   fChain->SetBranchAddress("puptmuloose1", &puptmuloose1, &b_puptmuloose1);
   fChain->SetBranchAddress("puptmuloose2", &puptmuloose2, &b_puptmuloose2);
   Notify();
}




void RedNtpFinalizer_VHgg::setSelectionType( const std::string& selectionType ) {

  selectionType_ = selectionType;


  if( selectionType=="sel0" ) {



  } else {

    std::cout << std::endl << std::endl << "Selection '" << selectionType << "' currently not implemented. Exiting." << std::endl;
    exit(12345);

  }


}
