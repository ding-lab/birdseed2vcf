#!/usr/bin/env python
"""
Given a birdseed data file, produce a VCF file for the target sample.

@author R Jay Mashl <rmashl@wustl.edu>
@date 2017-01-19
@version 1.2-DL1
@comment: expanded on header parsing; autodetermine sample ids and filenames from file metadata information

@author aaron <aaron@broadinstitute.org>
@date 4/27/2012
@version 1.2
@URL archive.broadinstitute.org/cancer/cga/sites/default/files/data/tools/contest/BirdseedToVCF.py

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
import commands

# --Global variables--
# set substrings for searching column headers
sStr = ['Aliquot Barcode', 'associated_entities__entity_submitter_id', 'Hybridization Name']
#sStr = ['Aliquot UUID', 'associated_entities__entity_submitter_id']

version = "1.2-DL1"

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
    # Ignore any comment/header lines that begin with hash
    # htz = shlex.shlex(annote.readline().strip("\n"))
    header = annote.readline()
    while header.startswith('#'):
        header = annote.readline()
    htz = shlex.shlex(header.strip("\n"))
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
    annote.close()

    pickle.dump( annotations, open( pickle_file, "wb" ) )
    return annotations

'''
get column number of string
'''
def get_column_num(line, s):
    cols = line.split("\t")
    if s in cols:
        return cols.index(s)
    else:
        for i in range(len(cols)):
            if s in cols[i]:
                return i
    return None

'''
set columns numbers based n sample mapping file
'''
def set_columns(mapfile, birdseedfile):
    global sStr
    myCols = {'file': None, 'id': None}

    for s in sStr:
        if myCols['id'] is None:
            cmd = "grep '%s' %s"  %  (s, mapfile)
            status, output = commands.getstatusoutput(cmd)
            lines = output.split("\n")
            if output != "":
                myCols['id'] = get_column_num( lines[0], s)
    if myCols['id'] is None:
        print "Error: did not recognize sample id from among the known choices %s"  %  (sStr)
        os._exit(1)

    # -- Data file --
    cmd = "grep %s %s"  %  (birdseedfile, mapfile)
    status, output = commands.getstatusoutput(cmd)
    lines = output.split("\n")
    myCols['file'] = get_column_num( lines[0], birdseedfile )
    if myCols['file'] is None:
        print "Error: did not locate birdseed file name in mapfile"
        os._exit(1)
    return myCols

'''
build mapping from file to sample id
'''
def get_mapping(mapfile, myCols):
    myMap = {}
    cmd = "/usr/bin/env perl -lane '@a=split/\t/;print join(\"\t\",@a[%d,%d])' %s"  %  (myCols['file'], myCols['id'], mapfile)
    status, output = commands.getstatusoutput(cmd)
    lines = output.split("\n")
    k=0
    for i in lines:
        c = i.split("\t")
        myMap[ c[0] ] = c[1]
        k += 1
        if k % 1000 == 0:
            print "[Status] Processed %d lines from the mapfile"  %  (k)
    return myMap

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
    global version
    now = datetime.datetime.now()
    file_hd.write("##fileformat=VCFv4.1\n##CreatedBy=BirdseedToVCF-v" + version + "\n##CreatedOn=")
    file_hd.write(now.strftime("%Y-%m-%d:%H:%M"))
    file_hd.write("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")

'''
figure out the call string
'''
def create_call(allele_a,allele_b,reference,call):
    if reference != allele_a and reference != allele_b:  #Chen SHH: for triallele?
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
        elif call == "0": #else:   #Chen SHH: modified on 20201/2/4 
            if reference.upper() == allele_a.upper():
                return "0/0"
            else:
                return "1/1"
        elif call == "-1":  #Chen SHH: modified on 20201/2/4
            return "./."
        else:
            pass
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
        order.close()
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
    parser.add_argument('--output_vcf', help='name of output VCF file created (supersedes mapfile)',required=False)
    parser.add_argument('--array_sample', help='the birdseed sample name in the birdseed file to process',required=True)
    parser.add_argument('--vcf_sample', help='the sample ID to place in the CHROM header line of produced VCF file (supersedes mapfile)',required=False)
    parser.add_argument('--mapfile', help='TCGA metadata file, mapping birdseed sample name/filename to sample/barcode; if given, --output_vcf and --vcf_sample are determined automatically (likely based on TCGA sample name) unless overridden',required=False)
    parser.add_argument('--snp_annotation_file', help='the annotation file from the array manufacturer describing the SNP positions and their alleles',required=True)
    parser.add_argument('--fasta', help='the reference fasta file',required=True)
    parser.add_argument('--add_chr', help='add prefix \'chr\' to the array contig name',required=False,action="store_true", default=False)
    parser.add_argument('--dont_convert_MT_to_M', help='convert the mitochondrial name of MT to M',required=False,action="store_true", default=False)
    parser.add_argument('--show_dropped_sites', help='if we drop a site due to missing information, output a line',required=False,action="store_true", default=False)
    parser.add_argument('--sort_on_the_fly', help='sometimes probes can come out of order, depending on the reference.  If this flag is set we store every record in memory, and sort when outputting (large memory required!)',required=False,action="store_true", default=False)
    parser.add_argument('--contig_ordering', help='if the contigs have a canonical ordering, we\'ll order the contigs by the ordering in this file. Each contig should be on a seperate line, in the order you\'d like',required=False)
    parser.add_argument('--drop_tri_allelic', help='drop all the tri-allelic sites.  Why arrays have these, who\'s to know. Drop them.',required=False,action="store_true", default=False)

    args = parser.parse_args()

    # check that our files are appropriate
    check_file_parameter(args.snp_annotation_file)
    check_file_parameter(args.birdseed)
    check_file_parameter(args.fasta)
    if args.mapfile:
        check_file_parameter(args.mapfile)
    else:
        if not args.output_vcf or not args.vcf_sample:
            raise NameError("Please specify both --output_vcf and --vcf_sample when --mapfile is not specified")

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

    # load filename mapfile
    mapfile_annotations = {}
    if args.mapfile:
        print "[Info] Search strings set as %s"  %  (sStr)
        print "[Status] analyzing files..."
        myCols = set_columns(args.mapfile, args.birdseed)
        print "[Info] Initially basing sample id and filename on mapfile %s using columns %d,%d (0-based), respectively"  %  (args.mapfile, myCols['id'], myCols['file'])
        print "[Status] loading the mapfile..."
        mapfile_annotations = get_mapping(args.mapfile, myCols)

    # now open up the birdseed file and the output,and convert the probes one by one
    birdseed = open(args.birdseed,"r")
    output = None

    # Parse header
    checking_header = 1
    samples_raw = []
    while checking_header:
        header = birdseed.readline().strip("\n").split("\t")
        if header[0].startswith('#') or header[0].startswith('Composite'):
            continue
        # check that their probeset file has the right parameters
        if header[0] == "probeset_id":
            if len(header) < 2:
                raise NameError("The genotypes (birdseed file) seems to have an unexpected header; we expect probeset_id followed by tab separated sample names")
            else:
                # Allow duplicate names for this header keyword to be handled as before below
                samples_raw = header[1:len(header)]
                checking_header = 0
        elif header[0].startswith('Hybridization'):
            if len(header) < 2:
                raise NameError("The genotypes (birdseed file) seems to have an unexpected header; we expect Hybridization REF followed by tab separated sample names")
            else:
                # currently assume single sample; have yet to see multisample
                samples_raw.append( header[1] )
                checking_header = 0

    if checking_header != 0:
        raise NameError("No sample names were found in the genotypes (birdseed file) header; no expected header keywords were found")

    sample_index = 0

    # find the sample name
    for sample in samples_raw:
        if sample == args.array_sample:
            if output != None:
                raise NameError("Duplicate names in the probe header; please fix this")
            # determine output filenames
            if args.mapfile:
                if not args.birdseed in mapfile_annotations:
                    print "Error: the birdseed filename %s was not found in the mapfile %s" % (args.birdseed, args.mapfile)
                    continue
                myOutputVcf = myVcfSample = mapfile_annotations[ args.birdseed ]
                myOutputVcf += ".vcf"
            if args.output_vcf:
                myOutputVcf = args.output_vcf
                if args.mapfile:
                    print "[Info] Overriding output VCF filename based by mapfile..."
            if args.vcf_sample:
                myVcfSample = args.vcf_sample
                if args.mapfile:
                    print "[Info] Overriding sample id indicated by mapfile..."

            print "[Info] Writing calls for sample %s to file %s ..."  %  (myVcfSample, myOutputVcf)
            output = open(myOutputVcf,"w")
            vcf_header(output,myVcfSample)
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

        if not args.dont_convert_MT_to_M:
            if chromosome == "chrMT":
                chromosome = "chrM"
            if chromosome == "MT":
                chromosome = "M"

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

    output.close()
    birdseed.close()
