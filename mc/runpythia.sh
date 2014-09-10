for f in $@
do
    gunzip -c $f/Events/run_01/events.lhe.gz > $f/events.lhe

    # TODO
    # why is this necessary?!
    perl -p -i -e "s/\"3.0\"/\"1.0\"/" $f/events.lhe

    run-pythia -n 999999999 -l $f/events.lhe -o $f/events.hepmc >& $f/pythia.log &

done
