# $1 -> number of events
# $2 -> input cmnd file
# $3 -> output folder

run-pythia -n $1 -i $2 -o $3/${2/cmnd/}13TeV.hepmc --collision-energy 13000 >& $3/${2/cmnd/}13TeV.log
