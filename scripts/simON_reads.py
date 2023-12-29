

# -*- coding: utf-8 -*-
"""

Python 3.10 or higher

@author: Astra Bertelli

Generate MinION-like long reads to store in artificial fastq files

"""
from argparse import ArgumentParser
from datetime import datetime
import sys
import os
import matplotlib.pyplot as plt
from simON_reads_module import *

argparse = ArgumentParser()
argparse.add_argument("-v", "--version", help="Print the version of the code", required=False, default=False, action='store_true')
argparse.add_argument("-i", "--infile", help="Path to the input fasta file that contains the reference sequence(s)", required=False)
argparse.add_argument("-snp", "--single_nucleotide_polymorphism", help="Insert a single nucleotide variant; Insert a single nucleotide variant; the syntax of this option should be SAMPLE:POS:REF>ALT,SAMPLE:POS:REF>ALT:1/0,...,SAMPLE:POS:REF>ALT:1/0 (it should be separated by commas without blank spaces) where SAMPLE is the header of the sequence (withouth \">\") in the original fasta file, POS is an integer that indicates the position (0-based) of the polymorphic site, REF is the reference allele, ALT is the alternative allele you want to be put and 1/0 (where you should report either 1 or 0, not both of them) is the haplotype phasing information: all the SNPs referred to 1 will end up on the same sequences, separate from the ones attributed to 0: this will generate a diploid-like distribution of variants. (Default is \"NO_SNP\")", required=False, default="NO_SNP")
argparse.add_argument("-n", "--nreads", help="Insert the number of reads you want to generate for each provided reference sequence: deafult is 2000", required=False, default=2000, type=int)
argparse.add_argument("-ehp", "--enable_homopolymer_error", help="This will set a 30%% chance of getting an extra nucleotide around homopolymeric regions", required=False, default=False, action='store_true')
argparse.add_argument("-ese", "--enable_sequencing_error", help="This will set a 5%% chance of getting a random single nucleotide variant or insertion, while it retains also a 5%% chance of skipping a base (single delition)", required=False, default=False, action='store_true')

args = argparse.parse_args()

inf = args.infile
snpstring = args.single_nucleotide_polymorphism
nr = args.nreads
enhompol= args.enable_homopolymer_error
enseqerr= args.enable_sequencing_error
vers = args.version

__version__="2.0.1"

if __name__=="__main__":
    if vers:
        print(f"Version: {__version__}")
        sys.exit()
    else:
        start = datetime.now()
        dic=load_data(inf) ## Load reference sequences from infile
        avgs=seqs_to_file(dic,snpstring,nr,enhompol,enseqerr) ## Generate reads and print them
        end=datetime.now()
        plt.style.use('ggplot') ## Plot average read quality as histogram, show it and save it as avg_quality_distribution.png in current directory
        fig,ax=plt.subplots(figsize=(10,5))
        ax.hist(avgs, bins=40, color='blue', edgecolor='black')
        ax.set_xlabel('Quality')
        ax.set_ylabel('Frequency')
        ax.set_title('Average read quality distribution')
        fig.show()
        fig.savefig(os.path.join(os.getcwd(),'avg_quality_distribution.png'))
        print('Duration: {}'.format(end - start), file=sys.stderr)
