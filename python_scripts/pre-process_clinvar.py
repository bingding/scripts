#!/usr/bin/python
# not used any more
import gzip,sys
    

fn=sys.argv[1]

prefix=fn.split(".gz")[0]
snp="%s.snp.txt"%prefix
fw_snp=open(snp,"w")

#indel="%s.indel.txt"%prefix
#fw_indel=open(indel,"w")

num_snp=0
num_indel=0
CLNSIG_Table={}

CLNSIG_Table={"0":  "Uncertain significance", "1": "not provided", "2": "Benign", "3": "Likely benign", "4": "Likely pathogenic", "5": "Pathogenic",\
              "6": "drug response", "7": "histocompatibility", "255":  "other"}


##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Variant Clinical Significance, 0 - Uncertain significance, 1 - not provided, 2 - Benign, 3 - Likely benign, 4 - Likely pathogenic, 5 - Pathogenic, 6 - drug response, 7 - histocompatibility, 255 - other">

for line in gzip.open(fn):
    if not line.startswith("#"):
        line=line.rstrip()
        tags=line.split("\t")
        chrom=tags[0]
        pos=tags[1]
        ID=tags[2]
        ref=tags[3]
        alt=tags[4]
        info=tags[7]
        CLNSIG=info.split("CLNSIG=")[1].split(";")[0]
        CLNACC=info.split("CLNACC=")[1].split(";")[0]
        if chrom =="MT":
            chrom="M"
        # only consider SNV
        if len(ref) ==1:
            if len(alt) ==1:
                print "single SNV ",chrom,pos,ref,alt,CLNSIG
                newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,alt,CLNSIG,CLNACC)
                fw_snp.write(newline)
                
            else: 
                if alt.rfind(",") != -1:
                    #clns=CLNSIG.split(",")
                    #clns_str=""
                    print "multi SNV ",chrom,pos,ref,alt,CLNSIG,CLNACC
                    snvs=alt.split(",")
                    sigs=CLNSIG.split(",")
                    accs=CLNACC.split(",")
                    for i in range(len(snvs)):
                        s=snvs[i]
                        if len(s) ==1:
                            if i <len(sigs):
                                sig=sigs[i]
                            else:
                                sig=sigs[0]
                            if i <len(accs):
                                acc=accs[i]
                            else:
                                acc=accs[0]
                                

                            #print "BBB ",chrom,pos,ref,s,sig,CLNACC,CLNSIG,alt,acc 
                            #newline
                            newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,s,sig,acc)
                            fw_snp.write(newline)
                        

fw_snp.close()
                


                   
        
