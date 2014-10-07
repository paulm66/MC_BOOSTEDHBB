#!/usr/bin/env python

from mg5common import mg5proc, mg5split

# -1 -> no max value
defptbs = [15, 25, 30, 40, 50, 75, 100, 125, 150, 200, 300, 500, 750, -1]

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
    splitprocs = map(lambda p: mg5split(p, "ptb", "ptbmax", defptbs), procs)

    # list of all procs to run
    procs = sum(splitprocs[1:], splitprocs[0])

    for p in procs:
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
