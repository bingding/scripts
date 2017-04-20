#!/usr/bin/python
import gzip,sys,os 
"""
seperate the SNV and INDEL of cosmic codint mutatation vcf file 
and then use bgzip to gzip it, tabix to index it 
not used any more
"""


fn=sys.argv[1]
prefix=fn.split(".vcf.gz")[0]
snp="%s.SNV.vcf"%prefix
fw_snp=open(snp,"w")
indel="%s.INDEL.vcf"%prefix
fw_indel=open(indel,"w")

num_snp=0
num_indel=0
num=0

for line in gzip.open(fn):
    if  line.startswith("#"):
        fw_snp.write(line)
        fw_indel.write(line)
    else:
        tags=line.split("\t")
        chrom=tags[0]
        pos=tags[1]
        ID=tags[2]
        ref=tags[3]
        alt=tags[4]
        info=tags[7]
        num +=1
        # SNP and INDEL in the same position
        if len(ref)==1 and len(alt)==1:
            fw_snp.write(line)
            num_snp +=1
        else:
            # newline="chr%s:%s\t%s\t%s\t%s\n"%(chrom,pos,ref,alt,ID)
            fw_indel.write(line)
            num_indel +=1
            
            
                
fw_snp.close()
fw_indel.close()
print "total mutations: ",num
print "total snp: ",num_snp
print "total indel: ",num_indel

print num_snp+num_indel
# now gzip it and use tabix to  index it
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/bgzip %s"%snp)
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/tabix -p vcf  %s.gz"%snp)
                
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/bgzip %s"%indel)
os.system("/mnt/pipeline-programs/tabix/tabix-0.2.6/tabix -p vcf  %s.gz"%indel)


