#!/usr/bin/python
"""
generate all the exonic coordinates for P30 transcript
Given a pre-defined  transcript ID, all exonic start and end position is retrieved
from refseq gene model file used in our pipeline-reference data
"""
import os,sys
from sets import  Set

def get_transcript_exons():
    # the refgene annotation file used in current production pipeline
    refseq_fn = "/mnt/pipeline-reference-data/cage_tables/gene_models_tables/gene_models_table_1/hg19_refSeq_genes.txt"
    transcript2exons={}
    transcript2Strand={}
    for line in open(refseq_fn):
        line = line.rstrip()
        tags = line.split("\t")
        gene = tags[9]
        chrom = tags[0]
        transcriptID=tags[10]
        strand=tags[1]
        transcript2Strand[transcriptID]=strand
        num_exons = int(tags[6])
        # delete the last ,
        exon_startPos = tags[7].split(",")[0:-1]
        # delete the last ,
        exon_endPos = tags[8].split(",")[0:-1]


        exon_regions = []
        # add the start-end position into a set
        for i in range(len(exon_startPos)):
            region = "%s-%s-%s" % (chrom, exon_startPos[i], exon_endPos[i])
            exon_regions.append(region)
        transcript2exons[transcriptID]=exon_regions
    return transcript2exons,transcript2Strand

# this file define the exons in which transcript covered by P30 panel
gene_file="P30_GOI.txt"
if not os.path.exists(gene_file):
    print "the P30 GOI file is not there"
    print "please make sure you have the GOI txt file in current folder"

    exit()
transcript2exons,transcript2Strand =get_transcript_exons()
# write out all the exonic coordinates
all_exon_file="P30_GOI_all_exons.bed"
fw=open(all_exon_file,"w")
added_transcripts=Set()
for line in open(gene_file).readlines()[1:]:
    line=line.rstrip()
    tags=line.split("\t")
    transcriptID=tags[6]
    gene_name=tags[0]

    exons=transcript2exons[transcriptID]
    strand=transcript2Strand[transcriptID]
    if not transcriptID in added_transcripts:
        added_transcripts.add(transcriptID)
        # exon num start from 1 for + strand
        exon_num=1
        # exon num start from the number of exons for - strand
        if strand=="-":
            exon_num=len(exons)
        # one line for one exon
        # the output format is chr, start_pos, end_pos, gene_name, transcriptID, strand, exon_num
        for e in exons:
            ex=e.replace("-","\t")
            output_line="%s\t%s\t%s\t%s\t%s\n"%(ex,gene_name,transcriptID,strand,exon_num)
            #print output_line
            fw.write(output_line)
            if strand=="+":
                exon_num +=1
            else:
                exon_num -=1


fw.close()