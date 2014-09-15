for f in $@
do
    echo "gunzip -c $f/Events/run_01/events.lhe.gz > $f/$f.lhe"
    gunzip -c $f/Events/run_01/events.lhe.gz > $f/$f.lhe

    # TODO
    # why is this necessary?!
    echo "perl -p -i -e "s/\"3.0\"/\"1.0\"/" $f/$f.lhe"
    perl -p -i -e "s/\"3.0\"/\"1.0\"/" $f/$f.lhe

    python splitLHE.py -i $f/$f.lhe -n 25000

    for f1 in `ls $f/${f}_*.lhe`
    do

        # 999999999 -> shower all events.
        echo "run-pythia -n 999999999 -l $f1 -o ${f1/lhe/hepmc} >& ${f1/lhe/pythia.log} &"
        run-pythia -n 999999999 -l $f1 -o ${f1/lhe/hepmc} >& ${f1/lhe/pythia.log} &

    done

done

echo "all jobs submitted in runpythia.sh. waiting for completion."

wait
