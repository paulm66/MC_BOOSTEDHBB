#!/usr/bin/env python

from mg5common import mg5proc, mg5split

# -1 -> no max value
defptbs = [0, 15, 25, 30, 35, 40, 50, 60, 75, 100, 125, 150, 175, 200,
        250, 300, 350, 400, 500, -1]

if __name__ == "__main__":
    from mg5procs import procdict
    from subprocess import Popen
    from time import sleep
    from sys import argv, stdout

    if len(argv) < 2:
        print "no processes specified."
        exit()

    elif "all" in argv:
        procs = procdict.values()

    else:
        procs = map(procdict.get, argv[1:])


    # list of list of split procs
    splitprocs = map(lambda p: mg5split(p, "ptj1min", "ptj1max", defptbs), procs)

    # list of all procs to run
    procs = sum(splitprocs[1:], splitprocs[0])

    for p in procs:
        p.initialize()

    # remove spurious leftover files
    Popen(["rm", "py.py"]).wait()

    map(mg5proc.generate_events, procs)

    print "all event generation jobs submitted."
