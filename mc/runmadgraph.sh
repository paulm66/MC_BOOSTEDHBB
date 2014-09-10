for f in $@
    do mg5_aMC $f
    done

rm py.py

for f in $@
    do echo "00" | ./${f/.mg5/}/bin/generate_events >& ${f/.mg5/}/run.log &
    done

wait
