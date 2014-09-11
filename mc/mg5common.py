# mg5common.py

from subprocess import Popen, PIPE
from sys import stdout

mg5header = """import model sm
define l+ e+ mu+ ta+
define l- e- mu- ta-
define l l+ l-
define v ve vm vt
define v~ ve~ vm~ vt~
define vall v v~
"""

# split files in ht of all partons
defhts = [0, 400, 800, 1600, 3200]

def runmg5(procname, proccmd, nevents, hts=[0]):
    for i in range(len(hts)):
        htmin = hts[i]
        try:
            htmax = hts[i+1]
        except IndexError:
            htmax = -1

        mg5 = Popen("mg5_aMC", stdin=PIPE, stdout=PIPE, stderr=PIPE)

        print "sending madgraph sm information:"
        print mg5header; stdout.flush()
        mg5.stdin.write(mg5header)

        print "sending madgraph process information:"
        print proccmd; stdout.flush()
        mg5.stdin.write(proccmd)

        foldername = "%s_htmin%04d%s"  % (procname, htmin,
                "_htmax%04d" % htmax if htmax != -1 else "")

        mg5.stdin.write("output madevent %s -f" % foldername)
        mg5.communicate("")

        print "fixing run_card.dat"; stdout.flush()
        runcard = open("%s/Cards/run_card.dat" % foldername, "r")
        lines = runcard.readlines()
        runcard.close()

        for i in range(len(lines)):
            line = lines[i]
            if "nevents" in line:
                lines[i] = "%04d = nevents ! Number of unweighted events requested\n" % nevents
            elif "= ihtmin" in line:
                lines[i] = "%04d = ihtmin !inclusive Ht for all partons (including b)\n" % htmin
            elif "= ihtmax" in line:
                lines[i] = "%04d = ihtmax !inclusive Ht for all partons (including b)\n" % htmax

            continue

        runcard = open("%s/Cards/run_card.dat" % foldername, "w")
        runcard.writelines(lines)
        runcard.close()

        continue

    return
