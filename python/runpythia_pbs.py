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

    # pythia doesn't like LHE version 3.0 for some reason...
    flhe = open(fnamelhe, 'w')
    lines[0] = lines[0].replace("3.0", "1.0")
    flhe.writelines(lines)
    flhe.close()

    # TODO
    # splitting turns out to be a bad idea.....
    # split the lhe file into smaller files with nevt events each
    # splitLHE.main(["-i", fnamelhe, "-n", str(nevt)])

    # flhenames = glob("%s/%s_*lhe" % (folder, folder))

    # running = []
    # for i in range(len(flhenames)):
    fin = fnamelhe
    fout = fin.replace("lhe", "hepmc")
    flog = fin.replace("lhe", "pythia.log")

    cmd = "run-pythia -n 999999999 -l %s -o %s" % \
            (fin, fout)

    # running.append(pbssubmit("pythia.%s.%d" % (folder, i), cmd,
        # outfile=flog))

    p = pbssubmit("pythia.%s" % fnamebase, cmd, outfile=flog)
    p.wait()
    print p.stdout.read()
    print p.stderr.read()
    stdout.flush()

    # map(lambda p: p.wait(), running)
    # for p in running:
        # print p.stdout.read()
        # print p.stderr.read()
    
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
