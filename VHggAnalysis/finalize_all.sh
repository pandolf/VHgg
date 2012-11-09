# usage: source finalize_all.sh [redntpProdVersion] [selectionType] [use_inverted_photID=\"false\"] [bTaggerType=\"JP\"]"


 cp finalize_TTVHgg finalize_TTVHgg_tmp
 ./finalize_TTVHgg_tmp  $1 DATA_Run2012ABC $2 $3
 ./finalize_TTVHgg_tmp  $1 HToGG_M-125_8TeV-pythia6 $2 $3
 ./finalize_TTVHgg_tmp  $1 DiPhoton_8TeV-pythia6 $2 $3
 ./finalize_TTVHgg_tmp  $1 V_8TeV $2 $3
 ./finalize_TTVHgg_tmp  $1 GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 ./finalize_TTVHgg_tmp  $1 QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2 $3
 ./finalize_TTVHgg_tmp  $1 TT_8TeV $2 $3 
 ./finalize_TTVHgg_tmp  $1 VV_8TeV $2 $3
 ./finalize_TTVHgg_tmp  $1 VGG_8TeV $2 $3
 rm finalize_TTVHgg_tmp
