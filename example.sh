#!/bin/bash

export PYTHONPATH=/usr/local/lib/python2.7/site-packages/

sample=GHOUL_p_TCGASNP_b85and51R_N_GenomeWideSNP_6_F04_735226
tcgaMeta=broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt
snpAnno=GenomeWideSNP_6.na35.annot.csv

python ./birdseed2vcf.py \
   --birdseed $sample.birdseed.data.txt \
   --snp_annotation_file $snpAnno \
   --array_sample $sample \
   --fasta GRCh37-lite.fa \
   --mapfile $tcgaMeta \
   --dont_convert_MT_to_M \
   --show_dropped_sites \
   --drop_tri_allelic
