from argparse import ArgumentParser
from datetime import datetime
import sys
import os
import matplotlib.pyplot as plt
from genfastqmodule import *

argparse = ArgumentParser()
argparse.add_argument("-i", "--infile", help="Path to the input fasta file that contains the reference sequence(s)", required=True)
argparse.add_argument("-snp", "--single_nucleotide_polymorphism", help="Insert a single nucleotide variant; the syntax of this option should be SAMPLE:POS:REF>ALT,SAMPLE:POS:REF>ALT,...,SAMPLE:POS:REF>ALT (it should be separated by commas without blank spaces) where SAMPLE is the header of the sequence (withouth \">\") in the original fasta file, POS is an integer that indicates the position (0-based) of the polymorphic site, REF is the reference allele, ALT is the alternative allele you want to be put. This will generate a diploid-like distribution of variants", required=False, default="NO_SNP")

args = argparse.parse_args()

inf = args.infile
snpstring = args.single_nucleotide_polymorphism


if __name__=="__main__":
    start = datetime.now()
    dic=load_data(inf)
    avgs=seqs_to_file(dic,snpstring)
    end=datetime.now()
    plt.style.use('ggplot')
    fig,ax=plt.subplots(figsize=(10,5))
    ax.hist(avgs, bins=40, color='blue', edgecolor='black')
    ax.set_xlabel('Quality')
    ax.set_ylabel('Frequency')
    ax.set_title('Average read quality distribution')
    fig.show()
    fig.savefig(os.path.join(os.getcwd(),'avg_quality_distribution.png'))
    print('Duration: {}'.format(end - start), file=sys.stderr)
