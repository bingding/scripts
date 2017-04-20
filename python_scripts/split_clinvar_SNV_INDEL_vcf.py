#!/usr/bin/python
"""
split the vcf file into SNV and INDEL vcf file 
then bgzip it and use tabix to index it
ensure bgzip and tabix are in your PATH
"""
import gzip,sys,os 
    

fn=sys.argv[1]

prefix=fn.split(".vcf.gz")[0]
snp="%s.SNV.vcf"%prefix
indel="%s.INDEL.vcf"%prefix
fw_snp=open(snp,"w")
fw_indel=open(indel,"w")

#indel="%s.indel.txt"%prefix
#fw_indel=open(indel,"w")

num_snp=0
num_indel=0
num=0

for line in gzip.open(fn):
    if line.startswith("#"):
        fw_snp.write(line)
        fw_indel.write(line)
    else:
        #line=line.rstrip()
        num +=1 
        tags=line.split("\t")
        chrom=tags[0]
        pos=tags[1]
        ID=tags[2]
        ref=tags[3]
        alt=tags[4]
        info=tags[7]
        # check whether it is SNV or INDEL
        if len(ref) ==1:
            if len(alt) ==1:
                fw_snp.write(line)
                num_snp +=1
            else:
                # it is multiple SNVs or multiple INDELs
                if alt.rfind(",") != -1:
                    is_SNV=True
                    alts=alt.split(",")
                    for a in alts:
                        if len(a)>1:
                            is_SNV=False
                            break
                    if is_SNV == True:
                        fw_snp.write(line)
                        num_snp +=1
                    else:
                        fw_indel.write(line)
                        num_indel +=1
                else:
                    # it is indel
                    fw_indel.write(line)
                    num_indel +=1
        # it is indel 
        else:
            fw_indel.write(line)
            num_indel +=1

fw_snp.close()
fw_indel.close()
print "total number of variants ",num
print "total number of snv ",num_snp
print "total number of indel ",num_indel
print num_snp+num_indel

# now gzip it and use tabix to  index it
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/bgzip %s"%snp)
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/tabix -p vcf  %s.gz"%snp)
                
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/bgzip %s"%indel)
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/tabix -p vcf  %s.gz"%indel)


                   
        
