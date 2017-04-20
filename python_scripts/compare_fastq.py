#!/usr/bin/python
#import HTSeq
#import Bio.SeqIO

import os,sys
from sets import Set
from parse_fastqorfasta import * 
f1=sys.argv[1]
f2=sys.argv[2]

results = []
seen = Set
## for s in Bio.SeqIO.parse(f1,"fastq"):
##     seen.add(str(s.seq))
##     print s
    
    

## for s in Bio.SeqIO.parse(f2,"fastq"):
##     if str(s.seq) not in seen:
##         results.append(x)

## output_handle = open("DIFF.fastq","w")
## SeqIO.write(results,output_handle,"fastq")
## output_handle.close()

# use htseq to read fastq files 
#reads1=HTSeq.FastqReader(f1,"phred")
#reads2=HTSeq.FastqReader(f2,"phred")

seqs1={}
seqs2={}

#reads1 = ParseFastQ(f1)
#reads2 = ParseFastQ(f2)
dscs_1=0
dscs_2=0
for name, seq, qual in readfq(f1):
     if  name.startswith("c"):
          seqs1[name.split("/")[0]]=[seq,qual]
          if name.count("#") == 5:
               dscs_1 +=1
     
for name, seq, qual in readfq(f2):
     if  name.startswith("c"):
          seqs2[name.split("/")[0]]=[seq,qual]
          if name.count("#") == 5:
               dscs_2 +=1

#    read_name=name.split("/")[0]
#    if not read_name in seqs1:
#        print read_name

     
#print len(reads1)
#print len(reads2)
#for read in reads1:
    # delete /1,2 
#    seqs1.add(read[0].split("/")[0])

#for read in reads2:
#    seqs2.add(read[0].split("/")[0])


print len(seqs1.keys()),len(seqs2.keys())
s1=Set(seqs1.keys())
s2=Set(seqs2.keys())

#for r in s2.difference(s1):
#     print r 

print "overlapped reads: ",len(s1.intersection(s2))
print "different reads: ",len(s1.difference(s2)),len(s2.difference(s1))
print "dscs reads 1: ",dscs_1
print "dscs reads 2: ",dscs_2 
          
diff=0
same=0
for s in s1.intersection(s2):
    if seqs1[s][0] == seqs2[s][0]:
         same +=1
         #print seqs1[s]
         #print seqs2[s]
    else:
         diff +=1
         #print s
         #print seqs1[s][0]
         #print seqs2[s][0]
         #print
print "Same sequence ", same," diff seq ",diff




        
