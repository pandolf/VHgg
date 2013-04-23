#!/bin/sh
#########################
#
# Driver script for Toy Monte Carlo submission with CRAB 
#
# author: Luca Lista, INFN
#                      
#########################

if [ -e outputToy ]; then 
  rm -rf outputToy 
fi
mkdir outputToy

i="$1"
if [ "$i" == "help" ]; then
  echo "usage: combine_crab.sh <job index> [<max events>]"
  exit 0;
fi
if [ "$i" = "" ]; then
  echo "Error: missing job index"
  exit 1;
fi
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents" 
if [ "$n" = "" ]; then
  n="$2"
fi
if [ "$n" = "" ]; then
  echo "Error: missing number of experiments"
  exit 2;
fi

# first, link locally the executable:
# ln -s ../../../../bin/slc5_amd64_gcc434/combine .

echo "job number: seed #$i"
#./combine datacard_VH_cat6-7_tth_wggH.dat -m125 -s$i -M HybridNew --freq  --saveToys --saveHybridResult -H ProfileLikelihood >& log.txt
./combine -M HybridNew --frequentist datacard_VH_cat6-7_tth_wggH.dat --expectedFromGrid=0.5 --saveToys --fullBToys --saveHybridResult -m125 -T400 -s$i -H ProfileLikelihood >& log.txt
#./combine d.txt -m125 -s$i -M HybridNew --freq  --saveToys --saveHybridResult -H ProfileLikelihood >& log.txt
mv higgsCombineTest.HybridNew.mH12*.root outputToy/
mv log.txt outputToy/
echo "pack the results"
tar cvfz outputToy.tgz outputToy/
