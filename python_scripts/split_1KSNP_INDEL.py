#!/usr/bin/python
import gzip,sys
"""
seperate the SNV and INDEL of 1K project
not used any more
"""


fn=sys.argv[1]
prefix=fn.split(".gz")[0]
snp="%s.snp.txt"%prefix
fw_snp=open(snp,"w")
indel="%s.indel.txt"%prefix
fw_indel=open(indel,"w")

num_snp=0
num_indel=0

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
        #SNP and INDEL in the same position
        if info.rfind("SNP") !=-1 and info.rfind("INDEL") !=-1 :
            #print line
            #
            if len(ref)==1:
                alts=alt.split(",")
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                afs=af.split(",")
                ids=ID.split(";")
                
                for i in range(len(alts)):
                    snp=alts[i]
                    af=afs[i]
                    if i <len(ids):
                        rsID=ids[i]
                    else:
                        rsID=ids[0]
                        #print line
                    # SNP
                    if len(snp)==1:
                        newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,snp,rsID,af)
                        #print line 
                        #print newline
                        fw_snp.write(newline)
                        num_snp +=1
                    # indel/insertion
                    else:
                        newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,snp,rsID,af)
                        #print line 
                        #print newline
                        fw_indel.write(newline)
                        num_indel +=1
                        
                        
            #1	11256034	rs17036508;rs533913726	TGTGA	CGTGA,T	100	PASS	AC=768,5;AF=0.153355,0.000998403;AN=5008;NS=2504;DP=18341;EAS_AF=0.126,0;AMR_AF=0.1225,0;AFR_AF=0.2194,0.0008;EUR_AF=0.0249,0.004;SAS_AF=0.2464,0;VT=SNP,INDEL;MULTI_ALLELIC
            # the SNP is the first position 
            else:
                alts=alt.split(",")
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                afs=af.split(",")
                ids=ID.split(";")
                
                for i in range(len(alts)):
                    snp=alts[i]
                    af=afs[i]
                    if i <len(ids):
                        rsID=ids[i]
                    else:
                        #take the first rsID
                        rsID=ids[0]
                        #print line," dddd"
                    
                    if len(snp)==len(ref):
                        real_snp=snp[0]
                        real_ref=ref[0]
                        newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,real_ref,real_snp,rsID,af)
                        #print line 
                        #print newline
                        fw_snp.write(newline)
                        num_snp +=1
                    # deletions
                    else:
                        newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,snp,rsID,af)
                        #print line 
                        #print newline
                        fw_indel.write(newline)
                        num_indel +=1

        #only SNP, it can be multi-allelic SNPS in the same position with different alt alle
        if info.rfind("SNP") !=-1 and info.rfind("INDEL") == -1:
            #print line 
            if len(alt)==1:
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,alt,ID,af)
                fw_snp.write(newline)
                num_snp +=1
            # mutli-allelic SNPs    
            else:
                alts=alt.split(",")
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                afs=af.split(",")
                
                for i in range(len(alts)):
                    snp=alts[i]
                    af=afs[i]
                    newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,snp,ID,af)
                    #print line 
                    #print newline
                    fw_snp.write(newline)
                    num_snp +=1
        #only INDEL,
        if info.rfind("SNP") ==-1 and info.rfind("INDEL") != -1:
            #print line
            #multiple indels
            if alt.rfind(",") != -1:
                alts=alt.split(",")
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                afs=af.split(",")
                
                for i in range(len(alts)):
                    indel=alts[i]
                    af=afs[i]
                    newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,indel,ID,af)
                    #print line 
                    #print newline
                    fw_indel.write(newline)
                    num_indel +=1
            #only one indel
            else:
                af_tag=info.split(";AF=")[1]
                af=af_tag.split(";")[0]
                newline="chr%s:%s\t%s\t%s\t%s\t%s\n"%(chrom,pos,ref,alt,ID,af)
                fw_indel.write(newline)
                num_indel +=1
                
fw_snp.close()
fw_indel.close()
print "total snp: ",num_snp
print "total indel: ",num_indel


