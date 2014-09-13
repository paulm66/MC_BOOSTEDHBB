import yoda
from sys import stdout

def haddFiles(lfin, fout):
    hists = {}
    for fin in lfin:
        dhists = yoda.read(fin)
        for path, hist in dhists.iteritems():
            if path in hists:
                print "hadding histogram %s" % path
                print "two integrals: %.2e, %.2e" % \
                        (hists[path].integral(), hist.integral())
                stdout.flush()
                hists[path] += hist
                print "integral after: %.2e" % hists[path].integral()
                stdout.flush()
            else:
                hists[path] = hist

            continue

        continue

    return yoda.writeYODA(hists, fout)

if __name__ == "__main__":
    from sys import argv
    haddFiles(argv[1:-1], argv[-1])
