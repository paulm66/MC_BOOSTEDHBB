#!/usr/bin/env python

defnevents = 1000000

def mg5split(minvar, maxvar, vals):
    procs = []
    for i in range(len(vals)):
        valmin = vals[i]
        try:
            valmax = vals[i+1]
        except IndexError:
            valmax = -1

        name = "%s_%s%04d%s"  % (proc.name, minvar, valmin,
                "_%s%04d" % valmax if valmax > 0 else "")

        # need at least a shallow copy here.
        d = dict(proc.runcarddict)

        d[minvar] = valmin
        d[maxvar] = valmax

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
