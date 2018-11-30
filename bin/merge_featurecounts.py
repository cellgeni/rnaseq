#!/usr/bin/env python

import argparse
import os
import sys
import re
import logging
from collections import defaultdict

def merge_featureCounts(dest_dir,out_file,input_files):

   logger = logging.getLogger(__name__)
   logger.addHandler(logging.StreamHandler())
   logger.setLevel(logging.INFO)
   
   table_dict=defaultdict(dict)
   sample_names=[]
   genes=set()
   num=1

   for input_file in input_files:
       logger.info("{:>5} Reading from {}".format(num, input_file))
       num += 1
       sample_name=os.path.basename(input_file)
       if (len(thesuffix)):
           sample_name = sample_name.replace(thesuffix,"")
           # featureCounts uses      Aligned.sortedByCoord.out_gene.featureCounts.txt
           # salmon uses             _1.quant.genes.sf
       sample_names.append(sample_name)
       table_dict[sample_name]=dict()
       must_read_header = expect_header
     
       with open(input_file, 'r') as f:
           while must_read_header:
               line = f.readline()
               if skip_comments and re.search("^\s*#", line):
                   continue
               else:
                   break        # only a single-line header supported currently, disregarding comments.

           for line in f:
               if skip_comments and re.search("^\s*#", line):
                   continue
               #save the genes to a list for the first file
               line_info=line.split('\t')
               gene=line_info[0]
               genes.add(gene)
               try:
                   gene_count = line_info[thecolumn].rstrip()
                   table_dict[sample_name][gene] = gene_count
               except TypeError:
                       logger.warning("Detected discrepancy in {}  line {}".format(input_file, line))
                
               table_dict[sample_name][gene]=gene_count

   #write Output
   logger.info("Writing to file {}".format(out_file))
   with open(out_file, 'w') as f:
       #Generate header
       line_to_write="ENSEMBL_ID"
       sample_names.sort()
       for sample_name in sample_names:
           line_to_write+=('\t{}'.format(sample_name))
       line_to_write += "\n"
       f.write(line_to_write)    
       #Write the rest of the lines
       gene_list=list(genes)
       gene_list.sort()
       for gene in gene_list:
           line_to_write=gene
           for sample_name in sample_names:
               try:
                   line_to_write += ('\t{}'.format(table_dict[sample_name][gene]))
               except KeyError:
                    #Missing gene in one of the filess.
                    line_to_write += ('\tNa')
           line_to_write += "\n"
           f.write(line_to_write)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Merges the counts for all the samples in a project
    """)
    parser.add_argument("-d", "--dest_dir", dest='dest_dir', default='.',
                                   help="Path to output.")
    parser.add_argument("-c", dest='colidx', default=-1, type=int,
                                   help="Index of column containing counts")
    parser.add_argument("--rm-suffix", default='',
                                   help="Remove suffix from file name to obtain sample name")
    parser.add_argument("--skip-comments", action="store_true",
                                   help="skip lines starting with #")
    parser.add_argument("--header", action="store_true",
                                   help="first non-commment line is a header line")
    parser.add_argument("-o", dest='out_file', default='all_counts.txt',
                                   help= "Name of the output file that will be created")
    parser.add_argument("-i", dest='input_files', metavar='<input_files>', nargs='+',
                                   help="Path to the outputfiles from FeatureCounts. ")
    parser.add_argument("-I", dest='metafile', default="", type=str,
                                   help="File containing input file names, one per line")
    args = parser.parse_args()
    skip_comments = args.skip_comments
    expect_header = args.header
    thesuffix = args.rm_suffix
    thecolumn = args.colidx
    metafile = args.metafile

    if len(metafile) > 0:
      if args.input_files:
        sys.exit("Do not mix -I and -i flags!")
      print(metafile)
      args.input_files = [line.rstrip('\n') for line in open(metafile)]

    merge_featureCounts(args.dest_dir, args.out_file, args.input_files)

