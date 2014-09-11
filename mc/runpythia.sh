for f in $@
do
    echo "gunzip -c $f/Events/run_01/events.lhe.gz > $f/$f.lhe"
    gunzip -c $f/Events/run_01/events.lhe.gz > $f/$f.lhe

    # TODO
    # why is this necessary?!
    echo "perl -p -i -e "s/\"3.0\"/\"1.0\"/" $f/$f.lhe"
    perl -p -i -e "s/\"3.0\"/\"1.0\"/" $f/$f.lhe

    # 999999999 -> shower all events.
    echo "run-pythia -n 999999999 -l $f/$f.lhe -o $f/$f.hepmc >& $f/pythia.log &"
    run-pythia -n 999999999 -l $f/$f.lhe -o $f/$f.hepmc >& $f/pythia.log &

done

echo "all jobs submitted in runpythia.sh. waiting for completion."

wait
