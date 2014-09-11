from sys import argv
from mg5common import runmg5, defhts

procname = "ttbar"
proccmd = """generate p p > t t~, t > l+ v b
add process p p > t t~, t~ > l- v~ b~
"""

if len(argv) > 2:
    nevents = int(argv[-1])
else:
    nevents = 10000

runmg5(procname, proccmd, nevents, defhts)
