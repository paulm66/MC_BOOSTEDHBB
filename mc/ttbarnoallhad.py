from mg5common import mg5proc, smheader

cmd = \
"""
generate p p > t t~, t > l+ v b
add process p p > t t~, t~ > l- v~ b~
"""

proc = mg5proc("ttbarnoallhad", smheader + cmd)

if __name__ == "__main__":
    from sys import argv
    if len(argv) > 2:
        proc.nevents(int(argv[-1]))
    proc.run()
