#!/usr/bin/python
"""
pre-process the cosmic mutant tab export tsv file to generate cosmic.txt file
in the anonotation table folder
currently not used. 
"""

import os,sys,gzip


fn=sys.argv[1]
output=sys.argv[2]
fw=open(output,"w")
def buildname2col(line):
    line=line.rstrip()
    colnames = line.split("\t")
    colname2colnum ={}
    colnum2colname={}
    for i,colname in enumerate(colnames):
        colname2colnum[colname]=i
        colnum2colname[i]=colname
    return colname2colnum,colnum2colname

line=gzip.open(fn).readline()

colname2colnum,colnum2colname=buildname2col(line)
#print colname2colnum
proteinChange2num={}
proteinChange2type={}
proteinChange2line={}
for line in gzip.open(fn):
    line=line.rstrip()
    if not line.startswith("Gene"): # skip the first header line
        tags=line.split("\t")
        geneName=tags[colname2colnum["Gene name"]]
        proteinChange=tags[colname2colnum["Mutation AA"]]
        cosmicID=tags[colname2colnum["Mutation ID"]]
        mutant="%s:%s"%(geneName,proteinChange)
        mutantion_type=tags[colname2colnum["Mutation Description"]]
        mutantion_type=mutantion_type.replace(" ","")
        proteinChange2type[mutant]=mutantion_type
        proteinChange2line[mutant]=line
        
        if proteinChange2num.has_key(mutant):
            proteinChange2num[mutant]+=1
        else:
            proteinChange2num[mutant]=1

print len(proteinChange2num.keys())
keys=list (proteinChange2num.keys())
keys.sort()

for k in keys:
    geneName=k.split(":")[0]
    proteinChange=k.split(":")[1]
    num=proteinChange2num[k]

    mutantion_type=proteinChange2type[k]
    start_pos=0
    end_pos=0
    #print mutantion_type
    if mutantion_type == "Unknown":
        continue
    if mutantion_type.startswith("Substitution"):
        change=proteinChange.split(".")[1]
        start_pos=change[1:-1]
        end_pos=start_pos
        newline="%s\t%s\t%s\t%s\t%s\t%s\n"%(geneName,start_pos,end_pos,num,proteinChange,mutantion_type)
        fw.write(newline)
        #print "%s\t%s\t%s\t%s\t%s"%(geneName,start_pos,end_pos,proteinChange,mutantion_type)
        
    elif mutantion_type.startswith("Deletion"):
        start_pos=0

        if mutantion_type.endswith("Frameshift"):
            change=proteinChange.split(".")[1].split("fs")[0]
            start_pos=change[1:]
            end_pos=start_pos
            newline="%s\t%s\t%s\t%s\t%s\t%s\n"%(geneName,start_pos,end_pos,num,proteinChange,mutantion_type)
            fw.write(newline)
            #print "%s\t%s\t%s\t%s\t%s\t%s"%(geneName,start_pos,end_pos,num,proteinChange,mutantion_type)
        elif mutantion_type.endswith("frame"):
            #print proteinChange
            if proteinChange[1]==".":
                change=proteinChange.split(".")[1].split("del")[0]
            else:
                change=proteinChange[1:]
            
            if change.rfind("_") !=-1:
                print change
                start_pos=change[1:]
                end_pos=start_pos
                newline="%s\t%s\t%s\t%s\t%s\t%s\n"%(geneName,start_pos,end_pos,num,proteinChange,mutantion_type)
                fw.write(newline)
                print "%s\t%s\t%s\t%s\t%s"%(geneName,start_pos,end_pos,proteinChange,mutantion_type)
            else:
                
                start_pos=change[1:]
                end_pos=start_pos
                newline="%s\t%s\t%s\t%s\t%s\t%s\n"%(geneName,start_pos,end_pos,num,proteinChange,mutantion_type)
                fw.write(newline)
                #print "%s\t%s\t%s\t%s\t%s"%(geneName,start_pos,end_pos,proteinChange,mutantion_type)

    elif mutantion_type.startswith("Insertion"):
        start_pos=0
    elif mutantion_type.startswith("Complex"):
        #print mutantion_type
        #print change
        #print "%s\t%s\t%s\t%s\t%s"%(geneName,start_pos,end_pos,proteinChange,mutantion_type)
        start_pos=0
        
    elif mutantion_type.startswith("Nonstop"):
        continue
        start_pos=0
        print geneName,proteinChange,mutantion_type
        print proteinChange2line[k]
