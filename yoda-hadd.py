import yoda
from sys import stdout

def haddFiles(lfin, fout):
    hists = {}
    for fin in lfin:
        dhists = yoda.read(fin)
        for path, hist in dhists.iteritems():
            if path in hists:
                hists[path] += hist
            else:
                hists[path] = hist

            continue

        continue

    return yoda.writeYODA(dhists, fout)

if __name__ == "__main__":
    from sys import argv
    haddFiles(argv[1:-1], argv[-1])
