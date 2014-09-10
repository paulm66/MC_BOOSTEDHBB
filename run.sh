INFILE=`echo $1 | sed -e 's/.*\///'`
OUTFILE=run/${INFILE/hepmc/yoda}
LOGFILE=${OUTFILE/yoda/log}

rivet -a MC_BOOSTEDHBB --pwd -H $OUTFILE $1 >& $LOGFILE &
