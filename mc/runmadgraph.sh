for f in $@
do
    mg5_aMC $f
done

rm py.py

for f in $@
do
    echo "echo \"00\" | ./${f/.mg5/}/bin/generate_events >& ${f/.mg5/}/madgraph.log &"
    echo "00" | ./${f/.mg5/}/bin/generate_events >& ${f/.mg5/}/madgraph.log &
done

echo "all jobs submitted in runmadgraph.sh. waiting for completion."

wait
