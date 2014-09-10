# $1 -> input madgraph folder
# $2 -> output hepmc file

gunzip -c $1/Events/run_01/events.lhe.gz > $1/events.lhe

# TODO
# why is this necessary?!
perl -p -i -e "s/\"3.0\"/\"1.0\"/" $1/events.lhe

run-pythia -n 999999999 -l $1/events.lhe -o $2 >& ${2/.hepmc/.log}
