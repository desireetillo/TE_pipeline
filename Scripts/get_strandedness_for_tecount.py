#!/usr/bin/env python3
"""
:Program Name: *get_strandedness_for_tecount.py*
:Developer:  Des Tillo
:email:  desiree.tillo@nih.gov
:Date: 2022-08-15

This script takes output from RNAseQC's infer_experiment.py and determines the 
strandedness parameter needed for TE_count (forward/reverse/no).  Modified 
from the RNAseq workflow from nf-core (https://github.com/nf-core/rnaseq)

"""

import os
import sys
import re
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", help="Required data input file from infer_experiment.py . ", required=True,dest="input")
    return parser.parse_args()


def inferstrandedness(infile):
    # initialize types
    cutoff=30
    sense = 0
    antisense = 0
    undetermined = 0
    undetermined_matcher = re.compile(r'Fraction of reads failed to determine:\s(\d.\d+)')
    se_sense_matcher     = re.compile(r'Fraction of reads explained by "\++,--":\s(\d.\d+)')
    se_antisense_matcher = re.compile(r'Fraction of reads explained by "\+-,-\+":\s(\d.\d+)')
    pe_sense_matcher     = re.compile(r'Fraction of reads explained by "1\++,1--,2\+-,2-\+":\s(\d.\d+)')
    pe_antisense_matcher = re.compile(r'Fraction of reads explained by "1\+-,1-\+,2\+\+,2--":\s(\d.\d+)')
    with open(infile) as file:
        for line in file:
            if undetermined_matcher.match(line):
                m=undetermined_matcher.match(line)
                undetermined=float(m.group(1)) * 100
            if se_sense_matcher.match(line):
                m=se_sense_matcher.match(line)
                sense=float(m.group(1)) * 100
            if se_antisense_matcher.match(line):
                m=se_antisense_matcher.match(line)
                antisense=float(m.group(1))*100
            if pe_sense_matcher.match(line):
                m=pe_sense_matcher.match(line)
                sense=float(m.group(1)) * 100
            if pe_antisense_matcher.match(line):
                m=pe_antisense_matcher.match(line)
                antisense=float(m.group(1))*100
    strandedness='no'
    if (sense >= 100-cutoff):
        strandedness = 'forward'
    elif (antisense >= 100-cutoff):
        strandedness = 'reverse'
    return strandedness;


def main():
    args=get_args()
    infile = args.input.rstrip("")
    strandedness=inferstrandedness(infile)
    print(strandedness)
    
if __name__ == "__main__":
    main()


