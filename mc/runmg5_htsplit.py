from mg5common import mg5proc
from mg5procs import procdict
from subprocess import Popen

defihts = [0, 400, 800, 1600, 3200]
defnevents = 100000

def ihtsplit(proc, ihts):
    procs = []
    for i in range(len(ihts)):
        ihtmin = ihts[i]
        try:
            ihtmax = ihts[i+1]
        except IndexError:
            ihtmax = -1

        name = "%s_ihtmin%04d%s"  % (proc.name, ihtmin,
                "_ihtmax%04d" % ihtmax if ihtmax != -1 else "")

        # need at least a shallow copy here.
        d = dict(proc.runcarddict)

        d["ihtmin"] = ihtmin
        d["ihtmax"] = ihtmax

        procs.append(mg5proc(name, proc.cmd, d))

        continue

    return procs


if __name__ == "__main__":
    from sys import argv
    if len(argv) < 2:
        print "no processes specified."
        exit()

    elif "all" in argv:
        procs = procdict.itervalues()

    else:
        procs = map(procdict.get, argv[1:])

    # list of list of ihtsplit procs
    ihtsplitprocs = map(lambda p: ihtsplit(p, defihts), procs)

    # list of all procs to run
    procs = sum(ihtsplitprocs[1:], ihtsplitprocs[0])

    for p in procs:
        p.nevents(defnevents)
        p.initialize()

    # remove spurious leftover files
    Popen(["rm", "py.py"]).wait()

    running = []
    for p in procs:
        # run on PBS cluster.
        # call p.generate_events() to run locally.
        running.append(p.generate_events(["--cluster"]))

    map(Popen.wait, running)

    print "all event generation complete."
