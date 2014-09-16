for f in $@
do
    echo "rivet --pwd -a MC_BOOSTEDHBB -H ${f/hepmc/yoda} $f >& ${f/hepmc/rivet.log} &"
    rivet --verbose --pwd -a MC_BOOSTEDHBB -H ${f/hepmc/yoda} $f >& ${f/hepmc/rivet.log} &
done

echo "all jobs submitted in runrivet.sh. waiting for completion."

wait
