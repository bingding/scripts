#!/usr/bin/python
import sys,os
from sets import Set
f=sys.argv[1]
seen=Set()
for line in open(f):
    line=line.rstrip()
    tags=line.split()
    ID=tags[0]
    pos=tags[2]
    if pos in seen:
        print pos
    else:
        seen.add(pos)
    
