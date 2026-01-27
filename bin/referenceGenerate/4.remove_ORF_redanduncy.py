####bin/referenceGenerate/4.remove_ORF_redanduncy.py
import argparse

import sys

def getOptions(args=sys.argv[1:]):

    parser = argparse.ArgumentParser(description="Parses command.")

    parser.add_argument("-i", "--input", help="Your input fasta file.")

    parser.add_argument("-o", "--outpath", help="Your output path.")

    parser.add_argument("-p", "--prefix",default='ORF_no_redundancy', help="Your output name prefix, default will be ORF_no_redundancy.")

    options = parser.parse_args(args)

    return options


options = getOptions(sys.argv[1:])

import pyfastx

import sys

import pysam

import re

import pandas as pd

import numpy as np


def insert(x):
  x.insert(2,'orf_merge_id',x['orf_id'].str.cat(sep=','));
  return x;

def remove_redundancy(input=options.input,outpath=options.outpath,prefix=options.prefix):
  fa=list(pyfastx.Fasta(input,build_index=False));
  fa=pd.DataFrame(fa,columns=('orf_id','seq'));
  fa=fa.groupby('seq').apply(insert);
  fa=fa.drop_duplicates(subset='orf_merge_id');
  with open(outpath+"/"+prefix+".fa","w") as orf_out:
    fa.apply(lambda x:orf_out.write(">"+x['orf_id']+"\n"+x['seq']+"\n"),axis=1);
  orf_out.close();
  with open(outpath+"/"+"redundancy_id.txt","w") as id_out:
    fa.apply(lambda x:id_out.write(x['orf_merge_id']+"\n"),axis=1);
  id_out.close();

if __name__ == '__main__':

  remove_redundancy();
