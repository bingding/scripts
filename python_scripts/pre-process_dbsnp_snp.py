#!/usr/bin/python
"""
this python script read the dbsnp SNV vcf file and generate snp.txt per chromosome
the snp.txt looks like 
rsID chrom Allele_A Allele_B

can be used to replace get_DBSNP.sh 
"""
import gzip, sys

fn = sys.argv[1]

prefix = fn.split(".gz")[0]
chroms = []
for i in range(1, 23):
    chroms.append("%s" % i)
chroms.append("X")
chroms.append("Y")
chroms2out = {}
reverse = {"C": "G", "G": "C", "A": "T", "T": "A"}

for line in gzip.open(fn):
    if not line.startswith("#"):
        line = line.rstrip()
        tags = line.split("\t")
        chrom = tags[0]
        pos = tags[1]
        ID = tags[2]
        ref = tags[3]
        alt = tags[4]
        info = tags[7]
        allNuc = []
        RV = False
        # INFO=<ID=RV,Number=0,Type=Flag,Description="RS orientation is reversed">
        if info.find("RV") != -1:
            RV = True;
        if len(ref) == 1:
            if RV:
                ref = reverse[ref]
            allNuc.append(ref)
            if len(alt) == 1:
                if RV:
                    alt = reverse[alt]
                allNuc.append(alt)
            else:
                ss = alt.split(",")
                for s in ss:
                    if len(s) == 1:
                        if RV:
                            s = reverse[s]
                        allNuc.append(s)

            # sort the nuclitides in A,C,G,T
            allNuc.sort()
            if len(allNuc) >= 2:
                n1 = allNuc[0]
                n2 = allNuc[1]
                output = "%s\tchr%s\t%s\t%s\t%s" % (ID, chrom, pos, n1, n2)
                if chroms2out.has_key(chrom):
                    chroms2out[chrom].append(output)
                else:
                    print chrom
                    chroms2out[chrom] = []
                    chroms2out[chrom].append(output)
for chrom in chroms:
    print chrom
    # write to snp.txt per chromosome 
    fw = open("hg19_snp_chr%s.txt" % chrom, "w")
    for s in chroms2out[chrom]:
        fw.write("%s\n" % s)
    fw.close()
