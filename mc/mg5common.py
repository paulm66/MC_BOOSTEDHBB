# mg5common.py

from subprocess import Popen, PIPE
from sys import stdout
import re

smheader = """import model sm
define l+ e+ mu+ ta+
define l- e- mu- ta-
define l l+ l-
define v ve vm vt
define v~ ve~ vm~ vt~
define vall v v~
"""

# split files in ht of all partons
defhts = [0, 400, 800, 1600, 3200]

class mg5proc:
    def __init__(self, name, cmd,
            runcarddict={"nevents": 10000}):

        self.name = name
        self.cmd = cmd
        self.runcarddict = runcarddict

        return


    def nevents(self, n=-1):
        if n < 0:
            return self.runcarddict["nevents"]
        else:
            self.runcarddict["nevents"] = n
            return n


    def run(self):

        mg5 = Popen("mg5_aMC", stdin=PIPE, stdout=PIPE, stderr=PIPE)

        print "sending madgraph process command."
        print self.cmd; stdout.flush()
        mg5.stdin.write(self.cmd)

        mg5.stdin.write("output madevent %s -f" % self.name)
        mg5.communicate("")

        print "fixing run_card.dat"; stdout.flush()

        runcard = open("%s/Cards/run_card.dat" % self.name, "r")
        lines = runcard.readlines()
        runcard.close()

        runcard = open("%s/Cards/run_card.dat" % self.name, "w")
        for line in lines:
            haskey = False
            for key, v in self.runcarddict.iteritems():
                if " %s " % key in line:
                    val = ("%s" % v).rjust(8)
                    runcard.write("%s = %s\n" % (val, key))
                    haskey = True
                    break

                continue

            if not haskey:
                runcard.write(line)

        runcard.close()

        return


def ihtsplit(proc, ihts):
    procs = []
    for i in range(len(ihts)):
        ihtmin = ihts[i]
        try:
            ihtmax = ihts[i+1]
        except IndexError:
            ihtmax = -1

        name = "%s_htmin%04d%s"  % (proc.name, htmin,
                "_htmax%04d" % htmax if htmax != -1 else "")

        # need at least a shallow copy here.
        d = dict(proc.runcarddict)

        d["ihtmin"] = ihtmin
        d["ihtmax"] = ihtmax

        procs.append(mg5proc(name, proc.cmd, d))

        continue

    return procs
