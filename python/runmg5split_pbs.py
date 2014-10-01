#!/usr/bin/env python

from mg5common import mg5proc, mg5split

# -1 -> no max value
defdrbbs = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0, -1]
defnevents = 25000

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

    for p in procs:
        p.runcarddict["xptb"] = 15
        p.runcarddict["etab"] = 4
        p.runcarddict["xptl"] = 15

    # list of list of split procs
    splitprocs = map(lambda p: mg5split(p, "drbb", "drbbmax", defdrbbs), procs)

    # list of all procs to run
    procs = sum(splitprocs[1:], splitprocs[0])

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
