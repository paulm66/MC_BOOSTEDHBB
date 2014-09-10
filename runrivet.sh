for f in $@
do
    echo "rivet -a MC_BOOSTEDHBB -H ${f/hepmc/yoda} --pwd $f >& ${f/hepmc/rivet.log} &"
    rivet -a MC_BOOSTEDHBB -H ${f/hepmc/yoda} --pwd $f >& ${f/hepmc/rivet.log} &
done

echo "all jobs submitted in runrivet.sh. waiting for completion."

wait
