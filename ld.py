#!/usr/bin/env python

import argparse, sys
# import math, time, re
import gzip
import numpy as np
import pandas as pd
import scipy.sparse as sparse
# from scipy import stats
# from collections import Counter
from argparse import RawTextHelpFormatter

__author__ = "Colby Chiang (colby.chiang@childrens.harvard.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2021-08-12 13:47 $"

# --------------------------------------
# define functions

def get_args():
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="\
ld.py\n\
author: " + __author__ + "\n\
version: " + __version__ + "\n\
description: calculate LD \n\
             from UKBB")
    parser.add_argument('--vars',
                        metavar='FILE', dest='vars_path',
                        type=str, default=None,
                        help='''
Tab-delimited file of variant pairs [stdin].')
Columns:
    1. chrom_a
    2. pos_a
    3. ref_a
    4. alt_a
    5. chrom_b
    6. pos_b
    7. ref_b
    8. alt_b
''')
    parser.add_argument('--npz',
                        metavar='FILE', dest='npz_file',
                        type=str, required=True,
                        help='LD npz file, from https://alkesgroup.broadinstitute.org/UKBB_LD')

    # parser.add_argument('-c', '--flagC',
    #                     required=False, action='store_true',
    #                     help='sets flagC to true')
    # parser.add_argument('input', nargs='?', type=argparse.FileType('r'),
    #                     default=None,
    #                     help='file to read. If \'-\' or absent then defaults to stdin.')


    # parse the arguments
    args = parser.parse_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.vars_path == None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)

    # send back the user input
    return args

def load_ld_npz(npz_file):
    ld_prefix = npz_file.removesuffix('.npz')
    
    #load the SNPs metadata
    gz_file = '%s.af.gz'%(ld_prefix)
    df_ld_snps = pd.read_table(gz_file, sep='\s+')
    df_ld_snps.rename(columns={'rsid':'SNP', 'chromosome':'CHR', 'position':'BP', 'allele1':'A1', 'allele2':'A2', 'aaf':'AAF'}, inplace=True, errors='ignore')
    assert 'SNP' in df_ld_snps.columns
    assert 'CHR' in df_ld_snps.columns
    assert 'BP' in df_ld_snps.columns
    assert 'A1' in df_ld_snps.columns
    assert 'A2' in df_ld_snps.columns
    assert 'AAF' in df_ld_snps.columns
    df_ld_snps.index = df_ld_snps['CHR'].astype(str) + '.' + df_ld_snps['BP'].astype(str) + '.' + df_ld_snps['A1'] + '.' + df_ld_snps['A2']

    # #load the SNP allele frequencies
    # df_freq_snps = pd.read_table(freq_file, sep='\s+')
    # df_freq_snps.rename(columns={}, inplace=True, errors='ignore')

    #load the LD matrix
    try:
        R = sparse.load_npz(npz_file).toarray()
        R += R.T
    except ValueError:
        raise IOError('Corrupt file: %s'%(npz_file))

    #create df_R and return it
    df_R = pd.DataFrame(R, index=df_ld_snps.index, columns=df_ld_snps.index)
    return df_R, df_ld_snps

# open file (either plaintext or zip)
def get_file(filename):
    if filename.endswith('.gz'):
        data = gzip.open(filename, 'rb')
    else:
        data = open(filename, 'r')
    return data    

# --------------------------------------
# main function

def main():
    # parse the command line args
    args = get_args()

    # if no input file, check if part of pipe and if so, read stdin.
    if args.vars_path == None:
        vars_file = sys.stdin
    else:
        vars_file = get_file(args.vars_path)

    print("loading npz...", file=sys.stderr)
    ld_mat = load_ld_npz(args.npz_file)
    print("done", file=sys.stderr)

    # print header
    print('\t'.join(['#chrA', 'posA', 'refA', 'altA', 'aafA', 'chrB', 'posB', 'refB', 'altB', 'aafB', 'eGene', 'slope', 'caddA', 'caddB', 'r', 'D']))
    
    # print variant pairs with LD line-by-line
    for line in vars_file:
        if line.startswith("#"):
            continue

        # map the columns to variables
        [chrA, posA, refA, altA, chrB, posB, refB, altB, eGene, slope, caddA, caddB] = line.rstrip().split('\t')[:12]

        try:
            # retrieve LD (as r) from npz file
            r = ld_mat[0]['.'.join((chrA, posA, refA, altA))]['.'.join((chrB, posB, refB, altB))]

            # retrive the allele frequencies for the variant pair
            aafA = float(ld_mat[1].loc['.'.join((chrA, posA, refA, altA)), 'AAF'])
            aafB = float(ld_mat[1].loc['.'.join((chrB, posB, refB, altB)), 'AAF'])

            # convert r to signed D, using allele frequencies
            # r^2 = D^2 / ( pA*(1-pA)*pB*(1-pB) )
            D = r * (aafA * (1 - aafA) * aafB * (1 - aafB)) ** 0.5

            # print the output
            print('\t'.join(map(str, (chrA, posA, refA, altA, aafA, chrB, posB, refB, altB, aafB, eGene, slope, caddA, caddB, r, "{:.8e}".format(D)))))

        except KeyError:
            print("Could not find variant in .npz file: %s" %('.'.join((chrA, posA, refA, altA)) + ', ' + '.'.join((chrB, posB, refB, altB))), file=sys.stderr)
            # raise IOError('Could not find variant in .npz file: %s'%('.'.join((chrA, posA, refA, altA)) + ', ' + '.'.join((chrB, posB, refB, altB))))
            
# initialize the script
if __name__ == '__main__':
    try:
        sys.exit(main())
    except IOError as e:
        if e.errno != 32:  # ignore SIGPIPE
            raise 
