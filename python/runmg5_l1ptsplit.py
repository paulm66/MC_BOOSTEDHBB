#!/usr/bin/env python

from mg5common import mg5proc

defptl1s = [0, 100, 200, 400]
ptl1splitdict = {
        "WH": defptl1s,
        "ZH": defptl1s,
        "Wbb": defptl1s,
        "Zbb": defptl1s,
        "ttbar": [0, 25, 35, 50, 75, 100, 150, 200, 250, 300, 400, 500]
        }

defnevents = 25000

def ptl1split(proc, ptl1s):
    procs = []
    for i in range(len(ptl1s)):
        ptl1min = ptl1s[i]
        try:
            ptl1max = ptl1s[i+1]
        except IndexError:
            ptl1max = -1

        name = "%s_ptl1min%04d%s"  % (proc.name, ptl1min,
                "_ptl1max%04d" % ptl1max if ptl1max > 0 else "")

        # need at least a shallow copy here.
        d = dict(proc.runcarddict)

        d["ptl1min"] = ptl1min
        d["ptl1max"] = ptl1max

        procs.append(mg5proc(name, proc.cmd, d))

        continue

    return procs


if __name__ == "__main__":
    from mg5procs import procdict
    from subprocess import Popen
    from time import sleep
    from sys import argv, stdout

    if len(argv) < 2:
        print "no processes specified."
        exit()

    elif "all" in argv:
        procs = procdict.itervalues()

    else:
        procs = map(procdict.get, argv[1:])

    # list of list of ptl1split procs
    ptl1splitprocs = map(lambda p: ptl1split(p, ptl1splitdict[p.name]), procs)

    # list of all procs to run
    procs = sum(ptl1splitprocs[1:], ptl1splitprocs[0])

    for p in procs:
        p.nevents(defnevents)
        p.initialize()

    # remove spurious leftover files
    Popen(["rm", "py.py"]).wait()

    running = []
    for p in procs:
        # run on PBS cluster.
        running.append(p.generate_events(["--cluster"]))

    nprocs = len(running)
    while nprocs > 0:
        nprocs = sum(map(lambda p: bool(p.poll()), running))
        print "%d processes still running." % nprocs
        stdout.flush()
        sleep(30)

    print "all event generation complete."
