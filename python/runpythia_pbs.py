#!/usr/bin/env python

from subprocess import Popen
from sys import argv, stdout
import splitLHE
from glob import glob
from pbssubmit import pbssubmit
from time import sleep

def runpythia_pbs(folder, nevt=25000):
    flhegz = "%s/Events/run_01/events.lhe.gz" % folder

    fnamebase = folder.split('/')[-1]
    fnamelhe = "%s/%s.lhe" % (folder, fnamebase)
    flhe = open(fnamelhe, 'w')

    cmd = ["gunzip", "-c", flhegz]
    print " ".join(cmd)

    p = Popen(cmd, stdout=flhe)
    p.wait()
    flhe.close()

    flhe = open(fnamelhe, 'r')
    lines = flhe.readlines()
    flhe.close()

    fin = fnamelhe
    fout = fin.replace("lhe", "hepmc")
    flog = fin.replace("lhe", "pythia.log")

    cmd = "run-pythia -n 999999999 -l %s -o %s" % \
            (fin, fout)

    p = pbssubmit("pythia.%s" % fnamebase, cmd, outfile=flog)
    p.wait()
    print p.stdout.read()
    print p.stderr.read()
    stdout.flush()

    return


def main(args):
    if len(argv) < 2:
        print "no madgraph output specified"
        exit()

    for f in args[1:]:
        runpythia_pbs(f.rstrip("/"))
        continue

    return 0


if __name__ == "__main__":
    main(argv)
