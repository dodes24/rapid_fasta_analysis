#!/bin/python
# -*- coding: utf-8 -*-
"""
This program takes a FASTA file and returns various analysis of the sequences.
"""


import sys
import urllib
import contextlib
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import transcribe


# fasta_file = input("Enter the name of the FASTA file: ")  # input file using terminal
fasta_file = "genomic.fna"


def parse_fasta(fasta_file):
    """
    This function takes a FASTA file and returns the ID, sequence, and length of the sequences.
    :param: fasta_file
    :return: seq_id, sequences, lengths
    """
    seq_id = []
    sequences = []
    lengths = []
    # created lists for the IDs, sequences, and lengths
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id.append(record.id)
        sequences.append(record)
        lengths.append(len(record))  # added the ids sequences and lengths of the sequence to the list

    return seq_id, sequences, lengths  # now works fine with one sequence in the file. Input file will be genomic sequence.


print(parse_fasta(fasta_file))


def gc_content(fasta_file):
    """
    This function calculates the GC content of the DNA sequence.
    :param: fasta_file
    :return: gc_cont_percent
    """
    seq_id, sequences, lengths = parse_fasta(fasta_file)
    for record in sequences:
        gc_cont_percent =  round(GC(str(record.seq)), 2)
        return gc_cont_percent


print(gc_content(fasta_file))


def transcribe_to_mrna(fasta_file):
    """
    This function transcribes the DNA sequence to mRNA.
    :param: fasta_file
    :return: mrna
    """
    seq_id, sequences, lengths = parse_fasta(fasta_file)
    for record in sequences:
        mrna = transcribe(str(record.seq))
        return mrna


print(transcribe_to_mrna(fasta_file))



# next step is to get the amino acid sequence, plots












"""
def gc_content(fasta_file):
    '''
    This function returns the GC content of the sequences.
    '''
    global all_nucleotides, gc_percent
    seq_id, sequences, lengths = parse_fasta(fasta_file)  # calling the function to get the IDs, sequences, and lengths
    gc_cont = []  # created a list for the GC content
    for record in sequences:
        gc_cont.append(str(record.seq).count("G") + str(record.seq).count("C"))  # counting the number of G and C
        # in the DNA sequence
        all_nucleotides = str(record.seq).count("A") + str(record.seq).count("T") + str(record.seq).count("G") + str(
            record.seq).count("C")  # counting the number of A, T, G, and C in the DNA sequence
        gc_percent = round((gc_cont[0] / all_nucleotides) * 100, 2)  # calculating the GC content in percent

    return gc_cont, all_nucleotides, gc_percent # returning the GC content, all nucleotides, and GC percent


print(gc_content(fasta_file))



def mRNA(fasta_file):
    '''
    This function makes mRNA from the DNA sequence.
    '''
    seq_id, sequences, lengths = parse_fasta(fasta_file)  # calling the function to get the IDs, sequences, and lengths
    mRNA_seq = []  # created a list for the mRNA sequences
    for record in sequences:
        mRNA_seq.append(str(record.seq).replace("T", "U"))  # replacing the T with U in the DNA sequence
    return mRNA_seq

print(mRNA(fasta_file))
"""