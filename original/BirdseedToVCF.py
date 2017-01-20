#!/usr/bin/env python
"""
Given a birdseed data file, produce a VCF file for the target sample.

@author aaron <aaron@broadinstitute.org>
@date 4/27/2012
@version 1.2

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the "Software"), to deal in the Software without restriction, including without limitation 
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE 
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""

import os
import sys
import subprocess
import argparse
import cPickle as pickle
import datetime
import shlex
from pyfasta import Fasta

'''
#####################################################
Utility functions -- see the bottom for the main loop
#####################################################
'''


'''
load up a SNP annotation file, and produce a mapping of header fields in the annotation file to their values for
each SNP probe. Store this association as a pickle file to speed up further runs
'''
def load_annotation_file(file,pickle_file):
    annote = open(file,"r")
    htz = shlex.shlex(annote.readline().strip("\n"))
    htz.whitespace += ","
    header = [k.strip().strip("\"") for k in htz]
    # print "header: " + "-",join(header)
    header_fields = {"allele_a": "Allele A", "allele_b": "Allele B", \
                         "chromosome" :"Chromosome", "pos" : "Physical Position",  "ma":  "Minor Allele", "sd": "Strand Versus dbSNP", \
                         "dbsnp": "dbSNP RS ID", "strand": "Strand"}
    if not "Probe Set ID" in header:
            raise NameError("Unable to find the required header field: Probe Set ID in the header of your annotations file.  Please check that this exists")
    probe_index = header.index("Probe Set ID")

    header_positions = {}
    for hf in header_fields.values():
        if not hf in header:
            raise NameError("Unable to find the required header field: " + hf + " in the header of your annotations file.  Please check that this exists")
        header_positions[header.index(hf)] = hf
    
    annotations = {}
    spi = shlex.shlex()
    
    # this part is slow -- but shlex allows for delimiters inside of quoted strings  
    for line in annote:
        spi = shlex.shlex(line.strip("\n"))
        spi.whitespace += ","
        dct = {}
        tk_index = 0
        name = ""
        for token in spi:
            if tk_index in header_positions:
                dct[header_positions[tk_index]] = token.strip("\"")
            elif tk_index == probe_index:
                name = token.strip("\"")
            tk_index += 1
        if name != "":
            annotations[name] = dct

    pickle.dump( annotations, open( pickle_file, "wb" ) )
    return annotations

'''
add some more direct output when a file doesn't exist
'''
def check_file_parameter(file):

    if not os.path.exists(file):
        print "-=-=-=-=-=-=-=-=-=- ERROR -=-=-=-=-=-=-=-=-=-"
        print file + " doesn't seem to exists, please correct this parameter"
        sys.exit(1)

'''
add a basic header to the VCF file we create
'''
def vcf_header(file_hd,sample):
    now = datetime.datetime.now()
    file_hd.write("##fileformat=VCFv4.1\n##CreatedBy=BirdseedToVCF-v1.1\n##CreatedOn=")
    file_hd.write(now.strftime("%Y-%m-%d:%H:%M"))
    file_hd.write("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")

'''
figure out the call string
'''
def create_call(allele_a,allele_b,reference,call):
    if reference != allele_a and reference != allele_b:
        if call == "0":
            return "1/1"
        elif call == "1":
            return "1/2"
        else:
            return "2/2"
    else:
        if call == "1":
            return "0/1"
        elif call == "2":
            if reference.upper() == allele_a.upper():
                return "1/1"
            else:
                return "0/0"
        else:
            if reference.upper() == allele_a.upper():
                return "0/0"
            else:
                return "1/1"
'''
change up the call, based on if we're reversed stranded or not
'''
def reverse_call(allele_a,allele_b,reverse):
    if reverse == "-":
        allele_a = reverse_base(allele_a)
        allele_b = reverse_base(allele_b)
    return (allele_a,allele_b)

def reverse_base(base):
    if base.upper() == "A":
        return "T"
    if base.upper() == "T":
        return "A"
    if base.upper() == "C":
        return "G"
    if base.upper() == "G":
        return "C"
    return base # if N

'''
output a single VCF line with the details from the probe annotations
'''
def vcf_line(file_hd,chromosome,pos,allele_a,allele_b,ref,call,name,probe_name,reverse):
    prev_a = allele_a
    prev_b = allele_b
    ret = reverse_call(allele_a,allele_b,reverse)
    allele_a = ret[0]
    allele_b = ret[1]
    # if the ref isn't the a or b allele, this is a special case 
    if ref != allele_a and ref != allele_b:
        line = chromosome + "\t" + pos + "\t" + name + "\t" + ref + "\t" + allele_a + "," + allele_b + "\t60\t.\t"
        line += "PID=" + probe_name + "\tGT:DP:GQ\t" + create_call(allele_a,allele_b,ref,call) + ":0:60\n"
    elif ref == allele_a:
        line = chromosome + "\t" + pos + "\t" + name + "\t" + ref + "\t" + allele_b + "\t60\t.\t"
        line += "PID=" + probe_name + "\tGT:DP:GQ\t" + create_call(allele_a,allele_b,ref,call) + ":0:60\n"
    else:
        line = chromosome + "\t" + pos + "\t" + name + "\t" + ref + "\t" + allele_a + "\t60\t.\t"
        line += "PID=" + probe_name + "\tGT:DP:GQ\t" + create_call(allele_a,allele_b,ref,call) + ":0:60\n"
    file_hd.write(line)

''' 
find the minor allele from the probe information
'''
def get_minor_allele(pop_field,pop):
    pop_split = [tk.strip().strip() for tk in pop_field.split("///")]
    for tk in pop_split:
        allele = [pp.strip() for pp in tk.split("//")]
        if len(allele) == 2 and allele[1] == pop:
                return allele[0]
    print pop_field
    return "N"
      
''' utility to check if a string is a number '''
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
      
''' add a vcf record information to the map '''
def add_record(mp,chromosome,position,allele_a,allele_b,ref,call,dbsnp,probe_name,reverse):
    if not mp.has_key(chromosome):
        mp[chromosome] = {}
    mp[chromosome][position] = (allele_a,allele_b,ref,call,dbsnp,probe_name,reverse)

''' output the stored records to the set output file '''
def save_records(output_file_handle,mp,contig_ordering):
    contig_names = mp.keys()
    if contig_ordering:
        print "loading contig file " + contig_ordering
        contig_names = []
        order = open(contig_ordering,"r")
        for line in order:
            contig_names.append(line.strip())
    for key in contig_names:
        if not mp.has_key(key):
            continue
        value = mp[key]
        keys = value.keys()
        keys.sort()
        for ckey in keys:
            cval = value[ckey]
            vcf_line(output_file_handle,key,str(ckey),cval[0],cval[1],cval[2],cval[3],cval[4],cval[5],cval[6])

'''
#########################################################################
the main guts of the program -- load up the annotation file, and then
go probe by probe, converting each into a VCF entry
#########################################################################
'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert a birdseed call file to a VCF for contamination input')
    parser.add_argument('--birdseed', help='the birdseed file, containing SNP positions and their genotype call column',required=True)
    parser.add_argument('--output_vcf', help='the output vcf file',required=True)
    parser.add_argument('--array_sample', help='the sample name to look for in the array (probe) file, this is converted to the output sample name for the vcf',required=True)
    parser.add_argument('--vcf_sample', help='the sample name used for the VCF file; this should be the same sample used in the BAM file in the contamination estimation',required=True)
    parser.add_argument('--snp_annotation_file', help='the annotation file, describing the SNP positions and their alleles from the array manufacturer',required=True)
    parser.add_argument('--fasta', help='the reference fasta file',required=True)
    parser.add_argument('--add_chr', help='sometime the reference sequence uses chr, where the arrays dont, in this case you want to add a \'chr\' to the array contig name',required=False,action="store_true", default=False)
    parser.add_argument('--dont_convert_MT_to_M', help='convert the mitochondrial name of MT to M',required=False,action="store_true", default=False)
    parser.add_argument('--show_dropped_sites', help='if we drop a site due to missing information, output a line',required=False,action="store_true", default=False)
    parser.add_argument('--sort_on_the_fly', help='sometimes probes can come out of order, depending on the reference.  If this flag is set we store every record in memory, and sort when outputting (large memory req!)',required=False,action="store_true", default=False)
    parser.add_argument('--contig_ordering', help='if the contigs have a canonical ordering, we\'ll order the contigs by the ordering in this file. Each contig should be on a seperate line, in the order you\'d like',required=False)
    parser.add_argument('--drop_tri_allelic', help='drop all the tri-allelic sites.  Why arrays have these, who\'s to know. Drop them.',required=False,action="store_true", default=False)

    args = parser.parse_args()

    # check that our files are appropriate
    check_file_parameter(args.snp_annotation_file)
    check_file_parameter(args.birdseed)
    
    # get the fasta file loaded
    print "[Status] loading the FASTA file..."
    fasta = Fasta(args.fasta, key_fn=lambda key: key.split(" ")[0])
    
    # now load up the snp annotation pickle if it exists, otherwise create it from the annotation file
    annotations = {}  
    print "[Status] loading the annotation file..."  
    if args.snp_annotation_file.endswith(".pickle"):
        pickle_file = args.snp_annotation_file
    else:
        pickle_file = args.snp_annotation_file + ".pickle"

    if (os.path.exists(pickle_file)):
        print "[Status] loading from pickled annotation file..."
        annotations = pickle.load( open( pickle_file, "rb" ))
    else:
        print "[Status] loading the raw annotation file (and pickling), this may take a while..."
        annotations = load_annotation_file(args.snp_annotation_file,pickle_file)
    
    # now open up the birdseed file and the output,and convert the probes one by one
    birdseed = open(args.birdseed,"r")
    output = None
    
    header = birdseed.readline().strip("\n").split("\t")

    # ignore any comment/header lines that begin with '#'
    while header[0].startswith('#'):
        header = birdseed.readline().strip("\n").split("\t")

    # check that their probeset file has the right parameters
    if len(header) < 2 or header[0] != "probeset_id":
        raise NameError("The genotypes (birdseed file) seems to have an unexpected header; we expect probeset_id followed by tab seperated sample names")
    samples_raw = header[1:len(header)]
    sample_index = 0

    # find the sample name
    for sample in samples_raw:
        if sample == args.array_sample:
            if output != None:
                raise NameError("Duplicate names in the probe header; please fix this")
            output = open(args.output_vcf,"w")
            vcf_header(output,args.vcf_sample)
            sample_index = samples_raw.index(sample) + 1
    if output == None:
        raise NameError("Unable to find the array sample name " + args.array_sample + " in the birdseed file header")

    print "[Status] processing birdseed lines..."
    total_lines = 0
    vcf_lines = {} # if we're sorting on the fly, use this to fill up with records
 
    for line in birdseed:
        sp = line.strip("\n").split("\t")
        if not annotations.has_key(sp[0]):
            if args.show_dropped_sites:
                print "Dropping site " + sp[0] + " because it's not in the annotations file"
            continue
        annot = annotations[sp[0]]
        if not is_number(annot["Physical Position"]):
            if args.show_dropped_sites:
                print "Dropping site " + sp[0] + " because the annotation seems to be wonky (position isn't a number)"
            continue

        position = int(annot["Physical Position"])
        chromosome = str(annot["Chromosome"])

        if args.add_chr:
            chromosome = "chr" + chromosome 
        
        if not args.dont_convert_MT_to_M and chromosome == "chrMT":
            chromosome = "chrM"

        ref_fasta = fasta[chromosome][position-1].upper()
        
        # output to the VCF file
        if ref_fasta == annot["Allele A"] or ref_fasta == annot["Allele B"] or not args.drop_tri_allelic:
            if not args.sort_on_the_fly:
                vcf_line(output,chromosome,str(position),annot["Allele A"],annot["Allele B"],ref_fasta,sp[sample_index],annot["dbSNP RS ID"],sp[0],annot["Strand"])
            else:
                add_record(vcf_lines,chromosome,position,annot["Allele A"],annot["Allele B"],ref_fasta,sp[sample_index],annot["dbSNP RS ID"],sp[0],annot["Strand"])
        
        total_lines += 1
        if total_lines % 10000 == 0:
            print "[Status] Processed " + str(total_lines) + " lines from the birdseed file"
    if args.sort_on_the_fly:
        save_records(output,vcf_lines,args.contig_ordering)
