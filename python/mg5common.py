# mg5common.py

from subprocess import Popen, PIPE
from pbssubmit import pbssubmit
from sys import stdout
import re

defnevents = 25000

defruncarddict = {
    "nevents": defnevents,
    "lhe_version": "1.0",
}

class mg5proc:
    def __init__(self, name, cmd, runcarddict=defruncarddict):

        self.name = name
        self.cmd = cmd
        self.runcarddict = runcarddict

        return


    def initialize(self):

        mg5 = Popen("mg5_aMC", stdin=PIPE, stdout=PIPE, stderr=PIPE)

        print "sending madgraph process commands:"
        print self.cmd; stdout.flush()
        mg5.stdin.write(self.cmd)
        mg5.stdin.write("\n")

        print "output madevent %s -f" % self.name; stdout.flush()
        mg5.stdin.write("output madevent %s -f\n" % self.name)

        print "quit"; stdout.flush()
        mg5.stdin.write("quit\n")

        mg5.wait()

        print "fixing run_card.dat"; stdout.flush()

        runcard = open("%s/Cards/run_card.dat" % self.name, "r")
        lines = runcard.readlines()
        runcard.close()

        runcard = open("%s/Cards/run_card.dat" % self.name, "w")
        for line in lines:
            haskey = False
            for key, v in self.runcarddict.iteritems():
                if " %s " % key in line:
                    runcard.write(" %s = %s\n" % (v, key))
                    haskey = True
                    break

                continue

            if not haskey:
                runcard.write(line)

            continue

        runcard.close()

        print "fixing me5_configuration.txt"; stdout.flush()
        meconfig = open("%s/Cards/me5_configuration.txt" % self.name, "w")
        meconfig.write("cluster_type=pbs\n")
        meconfig.write("cluster_queue=short6\n")
        meconfig.write("cluster_status_update=60 15\n")
        meconfig.close()


        print "%s initialized\n\n" % self.name
        stdout.flush()

        return


    def generate_events(self, opts=["--nb_core=1"]):
        print "starting event generation for %s." % self.name
        stdout.flush()

        cmd = " ".join(["./%s/bin/generate_events" % self.name] + opts + ["-f"])
        print cmd; print; stdout.flush()

        p = pbssubmit("madgraph.%s" % self.name, cmd,
                outfile="%s/generate_events.log")

        p.wait()
        print p.stdout.read()
        print p.stderr.read()
        stdout.flush()

        return


def mg5split(proc, minvar, maxvar, vals):
    procs = []
    for i in range(1, len(vals)):
        valmin = vals[i-1]
        valmax = vals[i]

        name = "%s_%s%.2f_%s%.2f" % \
                (proc.name, minvar, valmin, maxvar, valmax)

        # need at least a shallow copy here.
        d = dict(proc.runcarddict)

        d[minvar] = valmin
        d[maxvar] = valmax

        procs.append(mg5proc(name, proc.cmd, d))

        continue

    return procs
