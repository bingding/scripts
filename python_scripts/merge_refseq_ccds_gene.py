#!/usr/bin/python
"""
merge refseq_gene and ccdsgene downloaded from NCBI
"""
import os,sys,gzip

# it is refGene.txt.gz
refseq_fn=sys.argv[1]
# it is ccdsGene.txt.gz 
#ccdsgene_fn=sys.argv[2]
# mapping ccds id to NM Id
# it is ccdsInfo.txt.gz
ccdsInfo_fn=sys.argv[2]
line=""
NMID_2_ccds={}
NMID_no_version_2_ccds={}
# get the mapping of NM ID to ccdsID using ccdsInfo.txt.gz downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/ccdsInfo.txt.gz
#for line in gzip.open(ccdsInfo_fn):
#    line=line.rstrip()
#    if not line.startswith("#"):
#        tags=line.split("\t")
#        ccdsID=tags[0]
#        source=tags[1]
#        refID=tags[2]
#        if source == "N":
#            if refID.startswith("NM_"):
#                NMID=refID.split(".")[0]
#                NMID_2_ccds[NMID]=ccdsID

#if another version e.g. downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/ is used the script has to be adapted!
# now also the version of the transcript NM_ ID is taken into account for the mapping
for line in open(ccdsInfo_fn):
    line=line.rstrip()
    if not line.startswith("#"):
        tags=line.split("\t")
        ccdsID=tags[0]
        source=tags[3]
        refID=tags[4]
        status=tags[6]
        #print tags
        if source == "NCBI":
            if refID.startswith("NM_"):
                #print refID
                #print NMID_no_version
                NMID=refID
                NMID_2_ccds[NMID]=ccdsID
                NMID_no_version=refID.split(".")[0]
                NMID_no_version_2_ccds[NMID_no_version]=ccdsID

#print "length of NMID_2_ccds",len(NMID_2_ccds.keys())
fw=open("hg19_refSeq_genes.txt","w")
#for c in ccds2refseq.keys():
#    print c,ccds2refseq[c]

#now fo through refSeq file. new format without bin field and with version number of transcript
num_transcript=0
for line in gzip.open(refseq_fn):
    line=line.rstrip()
    if not line.startswith("#"):
        tags=line.split("\t")
        NMID=tags[0]+"."+tags[16]
        #print NMID
        NMID_no_version=tags[0]
        chrom=tags[1]
        strand=tags[2]
        UTR_start=tags[3]
        UTR_end=tags[4]
        first_exon_start=tags[5]
        last_exon_end=tags[6]
        numExons_str=tags[7]
        exonStart_positions=tags[8]
        exonEnd_positions=tags[9]
        # use following algorithm exchange first_exon_start position with the following of the exonStart_positions:
        # if corresponding entry in exonEnd_positions is the first entry larger than first_exon_start. 
        # Remove the coordinate pairs (start/end),for which entry in exonEnd_positions is smaller than first_exon_start
        # analogous for last_exon_end...  
        numExons=0
        if first_exon_start!=last_exon_end:    # in this case no UTR assigned and replacement, else: pseudogene etc.
            starts=exonStart_positions.split(",")
            ends=exonEnd_positions.split(",")
            exonStart_positions_NEW=[]
            exonEnd_positions_NEW=[]
            for index in range(len(starts)):
                #print x
                if starts[index]<last_exon_end and ends[index]>first_exon_start:
                    numExons+=1
                    #replace first/last exons:
                    if starts[index]<first_exon_start and ends[index]>first_exon_start and ends[index]<=last_exon_end:
                        exonStart_positions_NEW.append(first_exon_start)
                        exonEnd_positions_NEW.append(ends[index])
                    elif starts[index]<last_exon_end and ends[index]>last_exon_end and starts[index]>=first_exon_start:
                        exonStart_positions_NEW.append(starts[index])
                        exonEnd_positions_NEW.append(last_exon_end)
                    elif ends[index]>last_exon_end and starts[index]<first_exon_start:
                        exonStart_positions_NEW.append(first_exon_start)
                        exonEnd_positions_NEW.append(last_exon_end)
                        #print "special case"
                        #print NMID 
                    else:
                        exonStart_positions_NEW.append(starts[index])
                        exonEnd_positions_NEW.append(ends[index])
                        
            exonStart_positions=",".join(exonStart_positions_NEW)
            exonEnd_positions=",".join(exonEnd_positions_NEW)
            exonStart_positions+= ','
            exonEnd_positions+= ','
            numExons_str=str(numExons)
        #print exonStart_positions,exonEnd_positions
        #print
        geneName=tags[11]
        # map to ccds ID
        # but also check whether some mapping cannot be performed only due to version number of transcript
        ccdsID="."
        if NMID in NMID_2_ccds.keys():
            ccdsID=NMID_2_ccds[NMID]
        #else:
        #    if (len(chrom)<=5 and chrom !="chrM" and NMID.startswith("NM_")):
        #        if NMID_no_version in NMID_no_version_2_ccds.keys():
        #            print "no match due to version:"
        #            print NMID
        #            print NMID_no_version

        output="\t".join([chrom,strand,UTR_start,UTR_end,first_exon_start,last_exon_end,numExons_str,exonStart_positions,exonEnd_positions,geneName,NMID,ccdsID])
        # only consider chr1,chr2...chr22,chrX,chrY
        if (len(chrom)<=5 and chrom !="chrM" and NMID.startswith("NM_")):
            fw.write("%s\n"%output)
            num_transcript +=1
            #print line
            #print output

fw.close()
#print "total transcript ",num_transcript
