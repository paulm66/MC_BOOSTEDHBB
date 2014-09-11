from sys import argv
from mg5common import mg5proc, smheader

ttbarnoallhadcmd = smheader + 
"""
generate p p > t t~, t > l+ v b
add process p p > t t~, t~ > l- v~ b~
"""

ttbarnoallhadproc = mg5proc("ttbarnoallhad", ttbarnoallhadcmd)
