./combine -M HybridNew datacard_0001.root --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 5000 -i 20 -s13 >&! logHybrid_Signif_0001.log
./combine -M HybridNew datacard_001.root  --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 1000 -i 10 -s13 >&! logHybrid_Signif_001.log
./combine -M HybridNew datacard_01.root   --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 200 -i 10 -s13 >&! logHybrid_Signif_01.log
./combine -M HybridNew datacard_1.root    --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 200 -i 10 -s13 >&! logHybrid_Signif_1.log
./combine -M HybridNew datacard_2.root    --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 200 -i 10 -s13 >&! logHybrid_Signif_2.log
./combine -M HybridNew datacard_10.root   --significance  --expectedFromGrid=0.5  --saveToys --fullBToys --saveHybridResult -T 200 -i 10 -s13 >&! logHybrid_Signif_10.log

