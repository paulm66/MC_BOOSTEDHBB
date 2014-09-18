#!/usr/bin/env python

defnevents = 1000000

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
        p.nevents(defnevents)
        p.initialize()

    # remove spurious leftover files
    Popen(["rm", "py.py"]).wait()

    running = []
    for p in procs:
        opts = ["--nb-cores=1"]
        # run on PBS cluster.
        opts = ["--cluster"]
        running.append(p.generate_events(opts))

    nprocs = len(running)
    while nprocs > 0:
        print "%d processes still running." % nprocs
        stdout.flush()
        sleep(30)
        nprocs = sum(map(lambda p: bool(p.poll()), running))

    print "all event generation complete."
