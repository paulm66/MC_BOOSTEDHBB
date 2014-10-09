#!/usr/bin/env python

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
        p.initialize()

    # remove spurious leftover files
    Popen(["rm", "py.py"]).wait()

    map(mg5proc.generate_events, procs)


    print "all event generation jobs submitted."
