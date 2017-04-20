#!/usr/bin/python
"""
map ccds ID to hugo gene name 
not used any more
"""
import os,sys

# it is protein-coding_gene.txt
fn=sys.argv[1]
# it is ccdsGene.txt.gz 

line=""
ccds2hugo={}
lines=open(fn).readlines()

for line in lines[1:]:
    line=line.rstrip()
    tags=line.split("\t")
    if len(tags) >25:
        geneSymbol=tags[1]
        geneName=tags[2]
        ccdsID=tags[24]
        #refID=tags[23]
        if ccdsID.rfind("|") !=-1:
            ccdsID=ccdsID.replace('\"','')
            for c in ccdsID.split("|"):
                ccds2hugo[c]=geneSymbol
        else:
            if ccdsID.startswith("CCDS"):
                ccds2hugo[ccdsID]=geneSymbol
#print len(ccds2hugo)
ccdsIDs=list(ccds2hugo.keys())
ccdsIDs.sort()
for c in ccdsIDs:
    print "%s.1\t%s"%(c,ccds2hugo[c])
