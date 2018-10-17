
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
                  # to prevent numpy changing type size warnings.
                  # https://github.com/numpy/numpy/issues/11788

import collections
import pandas
import re
import argparse

parser = argparse.ArgumentParser(description="Summarise samtools idxstats output even further")
parser.add_argument("-m", dest = 'name_mito', required=True, default='MT', help="Name of mitochondrial chromosome")
parser.add_argument("-t", dest = 'name_table', required=True, default='', help="Name of file with samtools idxstats output")
args = parser.parse_args()

name_mito  = args.name_mito
name_table = args.name_table

                  # samtools idxstats output:
                  # seq-name(0), seq-length(1) #mapped(2) #unmapped(3)
                  # last row indexed *: only #unmapped

df = pandas.read_table(name_table, sep='\t', header=None)

re_ercc = re.compile("^ERCC-")

ids_ercc = [ i for i, word in enumerate(df[0]) if re_ercc.match(word) ]
id_unmapped = -1
id_mito = -1
for i, word in enumerate(df[0]):
  if word == '*':
    id_unmapped = i
  elif word == name_mito: 
    id_mito = i

if id_mito < 0:
  raise ValueError("Chromosome %s not found" % name_mito)
if id_unmapped < 0:
  raise ValueError("Unmapped tag '*' not found")


sum_mapped_ercc   = sum(df[2][ids_ercc])
sum_mapped_mito   = df[2][id_mito]
sum_mapped_other  = sum(df[2]) - sum_mapped_ercc - sum_mapped_mito
sum_unmapped      = df[3][id_unmapped]


thedict = collections.OrderedDict(
[("mapped_other"  ,   sum_mapped_other)
,("unmapped"      ,   sum_unmapped)
,("mapped_ercc"   ,   sum_mapped_ercc)
,("mapped_mito"   ,   sum_mapped_mito)
])

for item in thedict:
  num = thedict[item]
  print "%s\t%d"     %   (item, num)




