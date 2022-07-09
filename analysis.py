!/bin/python
-*- coding: utf-8 -*-
"""
This program takes a FASTA file and returns various analysis of the sequences.
"""

''' modules I going to use
import sys
import urllib
import contextlib
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
'''

from Bio import SeqIO

fasta_file = input("Enter the name of the FASTA file: ") # input file using terminal
def parse_fasta(fasta_file):
    '''
    This function takes a FASTA file and returns the ID, sequence, and length of the sequences.
    :param fasta_file:
    :return: seq_id, sequences, lenghts
    '''
    seq_id = []
    sequences = []
    lengths = []
    # created lists for the IDs, sequences, and lengthsq
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id.append(record.id)
        sequences.append(record)
        lengths.append(len(record)) # added the ids sequences and lengths of the sequence to the list

    return seq_id, sequences, lengths # now works fine with one sequence in the file. Input file will be genomic sequence.

# print(parse_fasta(fasta_file))

def mRNA(fasta_file):
    '''
    This function makes mRNA from the DNA sequence.
    '''
    seq_id, sequences, lenghts = parse_fasta(fasta_file) # calling the function to get the IDs, sequences, and lengths
    mRNA_seq = [] # created a list for the mRNA sequences
    for record in sequences:
        mRNA_seq.append(str(record.seq).replace("T", "U")) # replacing the T with U in the DNA sequence
    return mRNA_seq

# print(mRNA(fasta_file))