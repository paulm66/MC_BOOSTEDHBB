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


    def initialize(self):

        mg5 = Popen("mg5_aMC", stdin=PIPE, stdout=PIPE, stderr=PIPE)

        print "sending madgraph process command."
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
                    val = ("%s" % v).rjust(8)
                    runcard.write("%s = %s\n" % (val, key))
                    haskey = True
                    break

                continue

            if not haskey:
                runcard.write(line)

            continue

        runcard.close()

        return


    def generate_events(self):
        print "starting event generation for %s." % self.name
        stdout.flush()

        outf = open("%s/generate_events.log" % proc.name, 'w')
        evgenproc = Popen("%s/bin/generate_events" % proc.name,
                stdin=PIPE, stdout=outf, stderr=outf)

        evgenproc.stdin.write("0\n0\n")

        return evgenproc
