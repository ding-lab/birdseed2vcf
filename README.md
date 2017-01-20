# birdseed2vcf

Convert birdseed variant file with genotypes to a VCF file.

This enhanced version allows users to provide a mapfile/table containing 
birdseed filenames and corresponding sample identifiers (e.g., TCGA barcodes) for automatically
naming the output VCF files and sample ids therein, either of which may be superseded by
the appropriate options.

### Requirements

o Python
o python module pyfasta
o additional files (see, e.g., http://archive.broadinstitute.org/mpg/birdsuite/birdseed.html)


### Summary

birdseed2vcf.py [-h] --birdseed BIRDSEED [--output_vcf OUTPUT_VCF]
                     --array_sample ARRAY_SAMPLE [--vcf_sample VCF_SAMPLE]
                     [--mapfile MAPFILE] --snp_annotation_file
                     SNP_ANNOTATION_FILE --fasta FASTA [--add_chr]
                     [--dont_convert_MT_to_M] [--show_dropped_sites]
                     [--sort_on_the_fly] [--contig_ordering CONTIG_ORDERING]
                     [--drop_tri_allelic]


### Example conversion script (example.sh):

```sh
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
```

Refer to the global variables section in the python script to modify which strings are used for searching column headers in the mapfile.
