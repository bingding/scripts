#!/usr/bin/python

import pysam
import os,sys
from sets import Set

a=sys.argv[1]
b=sys.argv[2]


samA=pysam.AlignmentFile(a,"rb")
samB=pysam.AlignmentFile(b,"rb")

Areads=Set()
Breads=Set()

read=pysam.AlignedSegment()
for read in samA.fetch(until_eof=True):
    r1="read1"
    if read.is_read2:
        r1="read2"
    
    Areads.add("%s-%s"%(read.query_name,r1))
read=pysam.AlignedSegment()
#pysam.AlignedSegment.mate_is_unmapped
dif_num=0
for read in samB.fetch(until_eof=True):
    r1="read1"
    if read.is_read2:
        r1="read2"

    read_name="%s-%s"%(read.query_name,r1)
    Breads.add("%s-%s"%(read.query_name,r1))

    """
    if read_name not in Areads:
        #if read.mapq >0 and read.pos !=read.mpos:
        if read.mapq >0 and (not read.mate_is_unmapped):
            #print read
            #print read.mpos
            #print  ""
            dif_num +=1
    """

#print "diff num ",dif_num



    
print "total reads in ",a 
print len(Areads)
print "total reads in ",b  
print len(Breads)

print ""
print "A to B"
print "overlap: ", len(Areads.intersection(Breads))
print "difference: ", len(Areads.difference(Breads))

print "B to A"
print "overlap: ", len(Breads.intersection(Areads))
print "difference: ", len(Breads.difference(Areads))


