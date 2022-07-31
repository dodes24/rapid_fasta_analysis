#!/bin/python
# -*- coding: utf-8 -*-
"""
Title: rapid_fasta_analysis
Date: 2022-31-07
Author: Jozef Fülöp
Version: 1.0

Description:
This program is used to analyze the genomic sequence of a bacterial genome.
Bacteria is a prokaryotic organism, and therefore it does not have introns.
This program takes a FASTA file and returns various analysis of the sequence.
Output: two fasta files, one with mRNA sequences and some stats about the nucleotide sequence. The other with the
protein sequences and some stats abut the protein sequence.
Stats:
    nucleotide seq: GC content, ORFs, start-end, frame, sense, length, trailing
    protein seq: isoelectric point, molecular weight, instability index, secondary structure
    secondary structures: helix, strand, coil = number predicts probability of being a helix, strand, or coil
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis as Pa
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import transcribe
from orffinder import orffinder
from matplotlib import pyplot as plt
from collections import Counter
import textwrap
import numpy as np

fasta_file = input("Enter your FASTA file: ")  # input file using terminal
#  fasta_file = "tuberculosis_genomic.fna"


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


def isoelectric_point():  # prot analysis
    """
    This function returns the charge at pH 7.
    """
    iso_op_list = []
    for protein in genes_to_protein():
        protein_str = Pa(str(protein))
        iso_op_list.append(round(protein_str.isoelectric_point(), 2))
    return iso_op_list


def molecular_weight():
    """
    This function returns the molecular weight of the protein sequences.
    :return:
    """
    molecular_weight_l = []
    for protein in genes_to_protein():
        protein_str = Pa(str(protein))
        molecular_weight_l.append(round(protein_str.molecular_weight(), 2))
    return molecular_weight_l


def instability_index():
    """
    Thi function returns the instability index of the protein sequence.
    If the index is < 40, the protein is stable.
    """
    instability_index_list = []
    for protein in genes_to_protein():
        protein_str = Pa(str(protein))
        instability_index_list.append(round(protein_str.instability_index(), 2))
    return instability_index_list


def secondary_structure():
    """
    This function returns the prediction of secondary structure of the proteins.
    :return:
    """
    secondary_structure_list = []
    for protein in genes_to_protein():
        protein_str = Pa(str(protein))
        sec_struc = protein_str.secondary_structure_fraction()
        secondary_structure_list.append(sec_struc)

    return secondary_structure_list


def create_fasta_protein(e=None, f=None, g=None, h=None, i=None):
    """
    This function creates a FASTA file of protein sequences with header info abut nucleotide sequence
    from which the protein was translated. (start, stop, frame, sense, length, GC content...)
    """

    if f is None:  # if no file name is provided, create a default name
        f = isoelectric_point()  # create a list of GC content for each ORF
    if e is None:
        e = genes_to_protein()  # create a list of protein sequences
    if g is None:
        g = molecular_weight()  # create a list of molecular weights for each protein
    if h is None:
        h = instability_index()
    if i is None:
        i = secondary_structure()

    e_list_str = list(map(str, e))  # convert list of Seq objects to list of strings
    f_list_str = list(map(str, f))  # convert list of float objects to list of strings
    g_list_str = list(map(str, g))  # convert list of float objects to list of strings
    h_list_str = list(map(str, h))  # convert list of float objects to list of strings
    i_list_str = list(map(str, i))

    prot_number_list = np.arange(1, len(e_list_str) + 1)

    wrapper = textwrap.TextWrapper(width=80)

    with open("output/proteins.fasta", 'w') as fp:
        for j in range(len(e_list_str)):
            string = wrapper.fill(text=e_list_str[j])
            print("> PROTEIN NO." + str(prot_number_list[j]) + " |PROT SEQ INFO: IEP = " +
                  f_list_str[j] + ", MW = " + g_list_str[j] + ", Instability idx = " + h_list_str[j] +
                  ", 2nd. structure = " + i_list_str[j] +
                  "\n" + string + "\n", file=fp)
    return


def transcription_orf():
    """
    This function transcribes the ORFs to mRNA.
    :return:
    """
    mrna = []
    for orf in orf_nucleotides():
        mrna.append(transcribe(str(orf)))
    return mrna


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

    wrapper = textwrap.TextWrapper(width=80)

    with open("output/mrna_sequences.fasta", 'w') as fp:
        for j in range(len(e_list_str)):
            string = wrapper.fill(text=e_list_str[j])
            print("> GENE NO." + str(gene_number_list[j]) + " |NUC SEQ INFO| GC content: " +
                  f[j] + " %|" + d_list_str[j] + "\n" + string + "\n", file=fp)

    return


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
    plt.savefig("nucleotide_frequency_histogram.jpg", dpi=250)  # save the plot as a .jpg file
    # plt.show()  # show the plot

    return


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
    plt.savefig("nucleotide_frequency_pie.jpg", dpi=250)
    #  plt.show()

    return


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
    This function plots the amino acid frequency of the genome.
    """
    freq = protein_frequency()
    plt.figure(figsize=(10, 6))
    plt.bar(freq.keys(), freq.values())
    plt.title("Amino acid Frequency Distribution", fontsize=14)
    plt.xlabel("Amino acid", fontsize=14)
    plt.ylabel("Frequency", fontsize=14)
    plt.tight_layout()
    plt.savefig("protein_frequency.jpg", dpi=250)
    # plt.show()

    return


def protein_frequency_pie():
    """
    This function plots the  amino acid frequency of the genome.
    """
    freq = protein_frequency()
    plt.figure(figsize=(8, 6))
    plt.pie(freq.values(), labels=freq.keys(), autopct='%1.1f%%')
    plt.title("Amino acid Frequency", fontsize=14)
    plt.savefig("protein_frequency_pie.jpg", dpi=250)
    # plt.show()


protein_frequency_pie()
protein_frequency_histogram()
plot_nucleotide_frequency_pie()
plot_nucleotide_frequency_histogram()

print(overall_sequence_info())

create_fasta_nucleotide_mrna()
create_fasta_protein()
