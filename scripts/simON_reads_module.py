"""

Python 3.10 or higher

@author: Astra Bertelli

Functions to generate MinION-like long reads to store in artificial fastq files

"""

import gzip
import random as r
from math import ceil
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re



def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())



def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if infile.endswith(".fasta.gz") or infile.endswith(".fna.gz") or infile.endswith(".fas.gz") or infile.endswith(".fa.gz"):
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
        else:
            y.close()
            raise ValueError("File is the wrong format")
    # Read file directly as fasta if it is a not zipped fasta: handle also more uncommon extensions :-)
    elif infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fas") or infile.endswith(".fa"):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")

def ascii_conv_and_mean(line):
    phred_quality_dict = {
    '!' : 0, '"' : 1, '#' : 2, '$' : 3, '%' : 4, '&' : 5, "'" : 6, '(' : 7, ')' : 8, '*' : 9,
    '+' : 10, ',' : 11, '-' : 12, '.' : 13, '/' : 14, '0' : 15, '1' : 16, '2' : 17, '3' : 18, '4' : 19,
    '5' : 20, '6' : 21, '7' : 22, '8' : 23, '9' : 24, ':' : 25, ';' : 26, '<' : 27, '=' : 28, '>' : 29,
    '?' : 30, '@' : 31, 'A' : 32, 'B' : 33, 'C' : 34, 'D' : 35, 'E' : 36, 'F' : 37, 'G' : 38, 'H' : 39, 'I' : 40
    }
    mean_list = []
    for i in line:
        if i != '\n':
            mean_list.append(phred_quality_dict[i])
    return round(sum(mean_list)/len(mean_list), 3)


def quality_string(seq):
    """Produce random quality string with ASCII phred score code (from 1 to 30)"""
    keys =['!', '"', '#', '$', '%', '&', "'", '(', ')', '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I']
    qlt = []
    for i in seq:
        qlt.append(keys[int(round(r.random()*(len(keys)-1),0))])
    a=""
    quality = a.join(qlt)
    return quality

def generate_snp(snp_string, seq, header):
    """Insert SNP variation in the provided seq, starting from header and information provided with --single_nucleotide_polymorphism"""
    if snp_string == "NO_SNP":
        pass
    else:
        SNPs=snp_string.split(",")
        seq_snps = {}
        for i in SNPs:
            a=i.split(":")
            if a[0] == header: ## There is the possibility to insert every SNP referred to the sequence we are examining
                seq_snps.update({int(a[1]):[a[2].split(">")[0], a[2].split(">")[1], int(a[3])]})
        n = int(round(r.random(),0))
        for key in list(seq_snps.keys()):
            if seq[key] == seq_snps[key][0] and seq_snps[key][2]==n: ## insert SNP more or less 50% of the times
                seq=seq[:key]+seq_snps[key][1]+seq[key+1:]
            elif seq[key] != seq_snps[key][0] and seq_snps[key][2]==n: ## Warn that the information provided in the SNP string is not correct, so the outcome might not be the same as expected
                seq=seq[:key]+seq_snps[key][1]+seq[key+1:]
                print("WARNING! SNP in position " + str(key) + " for sequence " + str(header) + " do not match given genomic information; ignoring REF allele information...", file=sys.stderr)
            else:
                pass
    return seq




def homopolymers(s): 
    """Find homopolymeric regions, with 4 or more repeated nucleotides, such as AAAA, and return a dictionary that contains the starting site of each homopolymeric strain"""
    a = re.compile(r'A{4,}')
    amatches = a.finditer(s)
    t = re.compile(r'T{4,}')
    tmatches = t.finditer(s)
    c = re.compile(r'C{4,}')
    cmatches = c.finditer(s)
    g = re.compile(r'G{4,}')
    gmatches = g.finditer(s)
    return {'A':[match.start() for match in amatches],'T':[match.start() for match in tmatches],'C':[match.start() for match in cmatches],'G':[match.start() for match in gmatches]}

def perform_random_mutation(seq):
    """Perform a single nucleotide variant sequencing error or an indel error with 5% of the chance, perform a delition with 5% of the chance"""
    combinations = [
    'A', 'C', 'G', 'T',
    'AA', 'AC', 'AG', 'AT',
    'CA', 'CC', 'CG', 'CT',
    'GA', 'GC', 'GG', 'GT',
    'TA', 'TC', 'TG', 'TT',
    'AAA', 'AAC', 'AAG', 'AAT',
    'ACA', 'ACC', 'ACG', 'ACT',
    'AGA', 'AGC', 'AGG', 'AGT',
    'ATA', 'ATC', 'ATG', 'ATT',
    'CAA', 'CAC', 'CAG', 'CAT',
    'CCA', 'CCC', 'CCG', 'CCT',
    'CGA', 'CGC', 'CGG', 'CGT',
    'CTA', 'CTC', 'CTG', 'CTT',
    'GAA', 'GAC', 'GAG', 'GAT',
    'GCA', 'GCC', 'GCG', 'GCT',
    'GGA', 'GGC', 'GGG', 'GGT',
    'GTA', 'GTC', 'GTG', 'GTT',
    'TAA', 'TAC', 'TAG', 'TAT',
    'TCA', 'TCC', 'TCG', 'TCT',
    'TGA', 'TGC', 'TGG', 'TGT',
    'TTA', 'TTC', 'TTG', 'TTT']
    n=r.random()
    if n<=0.05: # SNVs and INDELs error
        g = ceil(r.random()*(len(seq)-1))
        seq=seq[:g]+combinations[ceil(r.random()*(len(combinations)-1))]+seq[g+1:]
        return seq
    elif 0.05<n<=0.1: # Delition error
        g = ceil(r.random()*(len(seq)-1))
        seq=seq[:g]+seq[g+1:]
        return seq
    else:
        return seq

def perform_homopolimer_mutation(seq, hp):
    """Insert am extra nucleotide to a homopolimeric site with 30% of the chances""" 
    ks=[]
    for i in list(hp.keys()): 
        if len(hp[i])>0: ## Check if there are homopolimers
            ks.append(i)
        else:
            pass
    ind=ceil(r.random()*(len(ks)-1)) ## Choose a random homopolimeric nucleotide 
    ind2=ceil(r.random()*(len(hp[ks[ind]])-1)) ## Choose a random homopolimeric site from that nucleotide
    n = r.random()
    if n <= 0.3:
        seq=seq[:ind2]+ks[ind]+seq[ind2+1:] ## Insert the extra nucleotide (from AAAAA > AAAAAA)
        return seq
    return seq


def seqs_to_file(genomes_dict, snp_string,nreads, enabled_hompolymer, enabled_sequencing_error):
    """Write a specified number of reads for each of the provided reference sequences (with a 5% of them being reverse complemented and a 0.5% being chimeric), inserting SNPs and random sequencing errors in them""" 
    genomes_list = [value for value in list(genomes_dict.values())]
    headers = [key for key in list(genomes_dict.keys())]
    seqs=[]
    for i in range(len(headers)):
        hp = homopolymers(genomes_list[i])
        revcomp=[int(nreads*r.random()) for k in range(int(nreads*0.05))] ## Add 5% of reverse complemented reads
        chimers=[int(nreads*r.random()) for k in range(int(nreads*0.005)) if k not in revcomp] ## Add 0.5% of chimeric reads 
        for j in range(nreads):
            if j in revcomp: ## Create reverse complemented reads
                seq=generate_snp(snp_string, genomes_list[i], headers[i])
                if enabled_sequencing_error:
                    seq=perform_random_mutation(seq)
                if enabled_hompolymer:
                    seq=perform_homopolimer_mutation(seq,hp)
                seqs.append(reverse_complement(seq))
            elif j in chimers: ## Create chimeric reads by merging two sequences
                n=ceil(r.random()*(len(genomes_list)-1))
                seq=genomes_list[n]+genomes_list[n-1]
                seqs.append(seq)
            else: ## Create a normal read
                seq=generate_snp(snp_string, genomes_list[i], headers[i])
                if enabled_sequencing_error:
                    seq=perform_random_mutation(seq)
                if enabled_hompolymer:
                    seq=perform_homopolimer_mutation(seq,hp)
                seqs.append(seq)
    means=[] ## Calculate the average base quality
    for i in range(len(seqs)):
        print("@seq"+str(i+1))
        print(seqs[i])
        print("+seq"+str(i+1))
        qt=quality_string(seqs[i])
        print(qt)
        means.append(ascii_conv_and_mean(qt))
    return means
    
