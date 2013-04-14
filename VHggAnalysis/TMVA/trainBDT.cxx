void trainBDT(TString name) {

    TFile *fOut = new TFile(name+".root","RECREATE");
    TMVA::Factory *factory = new TMVA::Factory(name, fOut, "!V:!Color");

    TFile *fSig = TFile::Open("THq_tHqLeptonic_mH125_8TeV_testtest_presel_isLeptonic_CSV.root");
    TTree *tSig = (TTree *) fSig->Get("tree_passedEvents");
    //tSig.AddFriend("tree_passedEvents", "THq_tHqLeptonic_mH125_8TeV_testtest_00_presel_CSV_friend.root");
    factory->AddSignalTree(tSig, 1.0);

    TCut cuts = "isLeptonic";

    //TFile *fBkg1 = TFile::Open("THq_TTbarGG_0Jet_Summer12-PU_S7_START52_V9-v1_presel_CSV.root");
    //TTree *tBkg1 = (TTree *) fBkg1->Get("tree_passedEvents");
    //tBkg1.AddFriend("BDT_lept/t", "THq_TTbarGG_0Jet_Summer12-PU_S7_START52_V9-v1_presel_CSV_plusBDT.root");
    //factory->AddBackgroundTree(tBkg1, 1.0);
    TFile *fBkg2 = TFile::Open("THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_isLeptonic_CSV.root");
    TTree *tBkg2 = (TTree *) fBkg2->Get("tree_passedEvents");
    //tBkg2.AddFriend("tree_passedEvents", "THq_TTH_HToGG_M-125_8TeV-pythia6_Summer12-PU_S7_START52_V9-v2_presel_CSV_friend.root");
    factory->AddBackgroundTree(tBkg2, 1.0);


    //VARIABLES
    factory->AddVariable("njets", 'I');
    //factory->AddVariable("nbjets_loose", 'I');
    factory->AddVariable("deltaPhi_top_higgs", 'F');
    factory->AddVariable("mt_top", 'F');
    factory->AddVariable("eta_qJet", 'F');
    factory->AddVariable("charge_lept", 'I');
    factory->AddVariable("deltaEta_lept_qJet", 'F');
    //factory->AddVariable("pt_bJet", 'F');


    factory->SetWeightExpression("1");
    factory->PrepareTrainingAndTestTree( cuts, cuts, "SplitMode=Random" );

    //factory->BookMethod( TMVA::Types::kLD, "LD", "!H:!V:VarTransform=None" );

    TString BDTGopt = "!H:!V:NTrees=500:BoostType=Grad:Shrinkage=0.05:!UseBaggedGrad:nCuts=200:nEventsMin=100:NNodesMax=5";

    BDTGopt += ":CreateMVAPdfs"; // Create Rarity distribution
    factory->BookMethod( TMVA::Types::kBDT, "BDTG", BDTGopt);


    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    fOut->Close();
}

