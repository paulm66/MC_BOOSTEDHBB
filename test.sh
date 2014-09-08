run-pythia \
    -n $1 \
    -o hepmc.fifo \
    --collision-energy 13000 \
    -c "HiggsSM:ffbar2HZ=on" \
    -c "HiggsSM:ffbar2HW=on" \
    -c "25:onMode=off" \
    -c "25:onIfMatch=5 -5" &

rivet -a MC_BOOSTEDHBB -H MC_BOOSTEDHBB.yoda \
    --pwd hepmc.fifo
