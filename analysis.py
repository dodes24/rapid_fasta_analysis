#!/bin/python
# -*- coding: utf-8 -*-
"""
This program is used to analyze the genomic sequence of a bacterial genome.
Bacteria is a prokaryotic organism and therefore it does not have introns.
This program takes a FASTA file and returns various analysis of the sequence.

"""
# from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import transcribe
from orffinder import orffinder
import numpy as np
from matplotlib import pyplot as plt
from collections import Counter

# fasta_file = input("Enter the name of the FASTA file: ")  # input file using terminal
fasta_file = "tuberculosis_genomic.fna"


def parse_fasta(input_file):
    """
    This function takes a FASTA file and returns the ID, sequence, and length of the sequences.
    :param: fasta_file
    :return: seq_id, sequences, lengths
    """
    seq_id = []
    sequences = []
    lengths = []
    # created lists for the IDs, sequences, and lengths
    for record in SeqIO.parse(input_file, "fasta"):
        seq_id.append(record.id)
        sequences.append(record)
        lengths.append(len(record))  # added the ids sequences and lengths of the sequence to the list

    return seq_id, sequences, lengths
    # now works fine with one sequence in the file. Input file will be genomic sequence.


def get_sequence():
    sequence_to_orf = SeqIO.read(fasta_file, "fasta")
    return sequence_to_orf


def gc_content(input_file):
    """
    This function calculates the GC content of the DNA sequence.
    :param: fasta_file
    :return: gc_cont_percent
    """
    seq_id, sequences, lengths = parse_fasta(input_file)
    for record in sequences:
        gc_cont_percent = round(GC(str(record.seq)), 2)
        return gc_cont_percent


def orf_finder():
    """
    This function finds the ORFs in the mRNA sequence.
    :param: mrna
    :return: open_reading_frames
    """
    open_reading_frames = orffinder.getORFs(get_sequence(), minimum_length=75, remove_nested=False)
    return open_reading_frames


def orf_nucleotides():
    """
    This function returns the nucleotide sequence of the ORFs.
    """
    orf_nucleotide_seq = orffinder.getORFNucleotides(get_sequence(), return_loci=False)
    return orf_nucleotide_seq


def list_seq_nucleotides():
    """
    This function returns list of the nucleotides of the ORFs.
    """
    seq_nucleotides = []
    for orf in orf_nucleotides():
        seq_nucleotides.append(str(orf))
    return seq_nucleotides


def translate_sequence(input_seq):
    """
    Wrapper for Biopython translate function.  Bio.Seq.translate will complain if input sequence is
    not a multiple of 3.  This wrapper function passes an acceptable input to Bio.Seq.translate in order to
    avoid this warning.
    """

    trailing_bases = len(str(input_seq)) % 3

    if trailing_bases:
        input_seq = ''.join([str(input_seq), 'NN']) if trailing_bases == 1 else ''.join([str(input_seq), 'N'])

    output_seq = Seq.translate(input_seq, to_stop=True, table="Bacterial")

    if trailing_bases:
        # remove last residue if input needed to be extended because of trailing bases
        output_seq = output_seq[:-1]

    return output_seq


def genes_to_protein():
    """
    This function translates the ORFs to proteins.
    returns proteins in a list of Seq objects.
    """
    proteins = []
    for orf in orf_nucleotides():
        translated = translate_sequence(orf)
        proteins.append(translated)
    return proteins


def gc_all_sequences():
    """
    This function returns the GC content of all sequences.
    :return:
    """
    all_gc_cont = []
    for gc in orf_nucleotides():
        all_gc_cont.append(str(round(GC(gc), 2)))
    return all_gc_cont


def create_fasta_protein(d=orf_finder(), e=None, f=None):
    """
    This function creates a FASTA file of protein sequences with header info abut nucleotide sequence
    from which the protein was translated. (start, stop, frame, sense, length, GC content...)
    """

    if f is None:  # if no file name is provided, create a default name
        f = gc_all_sequences()  # create a list of GC content for each ORF
    if e is None:
        e = genes_to_protein()  # create a list of protein sequences
    e_list_str = list(map(str, e))  # convert list of Seq objects to list of strings
    d_list_str = list(map(str, d))  # convert list of Seq objects to list of strings

    prot_number_list = np.arange(1, len(e_list_str) + 1)

    with open("proteins.fasta", 'w') as fp:
        for j in range(len(e_list_str)):
            print("> PROTEIN NO." + str(prot_number_list[j]) + " |NUCLEOTIDE SEQ INFO| GC content: " +
                  f[j] + " %|" + d_list_str[j] + "\n" + e_list_str[j], file=fp)

    return


create_fasta_protein()  # later on in the program


def transcription_orf():
    """
    This function transcribes the ORFs to mRNA.
    :return:
    """
    mrna = []
    for orf in orf_nucleotides():
        mrna.append(transcribe(str(orf)))
    return mrna


# print(transcription_orf())


def create_fasta_nucleotide_mrna(d=orf_finder(), e=None, f=None):
    """
    This function creates a FASTA file of mRNA sequences of found genes with header info abut nucleotide sequence.
    (start, stop, frame, sense, length, GC content...)
    """

    if f is None:
        f = gc_all_sequences()
    if e is None:
        e = transcription_orf()
    e_list_str = list(map(str, e))
    d_list_str = list(map(str, d))

    gene_number_list = np.arange(1, len(e_list_str) + 1)

    with open("mrna_sequences.fasta", 'w') as fp:
        for j in range(len(e_list_str)):
            print("> GENE NO." + str(gene_number_list[j]) + " |NUC SEQ INFO| GC content: " +
                  f[j] + " %|" + d_list_str[j] + "\n" + e_list_str[j], file=fp)

    return


create_fasta_nucleotide_mrna()


def plot_nucleotide_frequency_histogram():
    """
    This function plots the nucleotide frequency of the genome.
    """
    nuc_freq = {}
    for nuc in get_sequence():  # get the nucleotide sequence
        if nuc not in nuc_freq:  # if the nucleotide is not in the dictionary, add it
            nuc_freq[nuc] = 1  # add the nucleotide to the dictionary with a value of 1
        else:
            nuc_freq[nuc] += 1  # if the nucleotide is already in the dictionary, add 1 to its value
    nuc_freq_list = list(nuc_freq.values())  # convert the dictionary to a list
    nuc_freq_keys = list(nuc_freq.keys())  # convert the dictionary to a list
    plt.bar(nuc_freq_keys, nuc_freq_list)  # plot the list of values against the list of keys
    plt.title("Nucleotide Frequency", fontsize=14)  # set the title
    plt.xlabel("Nucleotide", fontsize=14)  # set the x-axis label
    plt.ylabel("Frequency", fontsize=14)  # set the y-axis label
    plt.tight_layout()  # set the layout of the plot
    plt.savefig("plot_nucleotide_frequency_histogram.jpg", dpi=250)  # save the plot as a .jpg file
    plt.show()  # show the plot

    return


plot_nucleotide_frequency_histogram()


def plot_nucleotide_frequency_pie():
    """
    This function plots the nucleotide frequency of the genome.
    """
    nuc_freq = {}
    for nuc in get_sequence():
        if nuc not in nuc_freq:
            nuc_freq[nuc] = 1
        else:
            nuc_freq[nuc] += 1
    nuc_freq_list = list(nuc_freq.values())
    nuc_freq_keys = list(nuc_freq.keys())
    plt.figure(figsize=(10, 6))
    plt.pie(nuc_freq_list, labels=nuc_freq_keys, autopct='%1.1f%%')
    plt.title("Nucleotide Frequency", fontsize=14)
    plt.savefig("plot_nucleotide_frequency_pie.jpg", dpi=250)
    plt.show()
    return


plot_nucleotide_frequency_pie()


def overall_sequence_info():
    """
    This function returns the overall sequence information of the genome.
    """
    nuc_freq = {}
    for nuc in get_sequence():
        if nuc not in nuc_freq:
            nuc_freq[nuc] = 1
        else:
            nuc_freq[nuc] += 1
    nuc_freq_list = list(nuc_freq.values())
    nuc_freq_keys = list(nuc_freq.keys())
    nuc_freq_sum = sum(nuc_freq_list)
    nuc_freq_percent = list(map(lambda x: round(x / nuc_freq_sum * 100, 2), nuc_freq_list))
    nuc_freq_percent_keys = list(map(lambda x: str(x) + "% ", nuc_freq_percent))
    final_list_freq = str([item for sublist in zip(nuc_freq_keys, nuc_freq_percent_keys) for item in sublist])
    nuc_freq_percent_keys_str = ''.join(final_list_freq)
    print(f"Sequence length: {nuc_freq_sum} bp")
    print(f"Nucleotide Percentage: {nuc_freq_percent_keys_str}")
    print(f"GC Content: {round((nuc_freq_list[1] + nuc_freq_list[3]) / nuc_freq_sum * 100, 2)}%")
    print(f"AT Content: {round((nuc_freq_list[2] + nuc_freq_list[0]) / nuc_freq_sum * 100, 2)}%")
    return


print(overall_sequence_info())


def dna_to_protein():
    """
    This function translates the DNA sequence to protein.
    """
    protein = get_sequence().translate()
    return protein


def protein_frequency():
    """
    This function returns the frequency of each amino acid in the protein sequence.
    :return:
    """
    prt_freq = Counter(dna_to_protein())
    return prt_freq


def protein_frequency_histogram():
    """
    This function plots the protein frequency of the genome.
    """
    freq = protein_frequency()
    plt.figure(figsize=(10, 6))
    plt.bar(freq.keys(), freq.values())
    plt.title("Amino acid Frquency Distribution", fontsize=14)
    plt.xlabel("Amino acid", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.tight_layout()
    plt.savefig("protein_frequency.jpg", dpi=250)
    plt.show()


protein_frequency_histogram()


def protein_frequency_pie():
    """
    This function plots the protein frequency of the genome.
    """
    freq = protein_frequency()
    plt.figure(figsize=(8, 6))
    plt.pie(freq.values(), labels=freq.keys(), autopct='%1.1f%%')
    plt.title("Amino acid Frequency", fontsize=14)
    plt.savefig("protein_frequency_pie.jpg", dpi=250)
    plt.show()


protein_frequency_pie()


# Todo: add stat output in right format, find more stats to add about hypothetic proteins I found, graphs, etc.