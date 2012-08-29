# usage: source finalize_all.sh [selectionType]


 cp finalize_VHgg finalize_VHgg_tmp
 ./finalize_VHgg_tmp  HToGG_M-125_8TeV-pythia6 $1
 ./finalize_VHgg_tmp  WH_ZH_HToGG_M-125_8TeV-pythia6 $1
 ./finalize_VHgg_tmp  DiPhoton_8TeV-pythia6 $1
 ./finalize_VHgg_tmp  V_8TeV $1
 ./finalize_VHgg_tmp  GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $1
 ./finalize_VHgg_tmp  QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $1
 ./finalize_VHgg_tmp  TT_8TeV $1
 ./finalize_VHgg_tmp  VV_8TeV $1
 ./finalize_VHgg_tmp  VGG_8TeV $1
 rm finalize_VHgg_tmp
