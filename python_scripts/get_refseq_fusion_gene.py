#!/usr/bin/python
"""
get gene cooridates for fusion detect

"""
import os,sys

# it is hg19_refSeq_genes.txt
fn=sys.argv[1]

line=""
posgene={}

for line in open(fn):
    line=line.rstrip()
    tags=line.split("\t")
    chrom=tags[0]
    start=tags[2]
    end=tags[3]
    geneName=tags[9]
    pos_gene="%s\t%s\t%s\t%s"%(chrom,start,end,geneName)
    posgene[pos_gene]=1

#print len(pos2gene)

for p in posgene.keys():
    print "%s"%(p)

