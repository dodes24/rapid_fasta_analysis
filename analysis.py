#!/bin/python
# -*- coding: utf-8 -*-
"""
This program takes a FASTA file and returns various analysis of the sequences.
"""

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Seq import transcribe
from orffinder import orffinder


# fasta_file = input("Enter the name of the FASTA file: ")  # input file using terminal
fasta_file = "tuberculosis_genomic.fna"


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


#print(parse_fasta(fasta_file))


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


def transcribe_to_mrna():
    """
    This function transcribes the DNA sequence to mRNA.
    :return: mrna
    """
    seq_id, sequences, lengths = parse_fasta(fasta_file)
    for record in sequences:
        mrna = transcribe(str(record.seq))
        return mrna


# print(transcribe_to_mrna(fasta_file))

def get_sequence():
    sequence_to_orf = SeqIO.read(fasta_file, "fasta")
    return sequence_to_orf

def orf_finder():
    '''
    This function finds the ORFs in the mRNA sequence.
    :param: mrna
    :return: orfs
    '''
    orfs = orffinder.getORFs(get_sequence(), minimum_length=75, remove_nested=False)  # remove_nested=True)
    return orfs


print(orf_finder()[-1])


def orf_nucleotides():
    '''
    This function returns the nucleotides of the ORFs.
    '''
    orf_nucleotides = orffinder.getORFNucleotides(get_sequence(), return_loci=False)
    return orf_nucleotides
'''
    for x in orf_nucleotides:
        if (len(str(x))) % 3 != 0:
           printed = print(len(x))
        #printed = print(len(str((x))))
    return printed
'''
#print(orf_nucleotides()[0])

def list_seq_nucleotides():
    '''
    This function returns list of the nucleotides of the ORFs.
    '''
    seq_nucleotides = []
    for orfs in orf_nucleotides():
        seq_nucleotides.append(str(orfs))
    return seq_nucleotides

#print(list_seq_nucleotides())



def translate_sequence(input_seq):
    """Wrapper for Biopython translate function.  Bio.Seq.translate will complain if input sequence is
    not a mulitple of 3.  This wrapper function passes an acceptable input to Bio.Seq.translate in order to
    avoid this warning."""

    trailing_bases = len(str(input_seq)) % 3

    if trailing_bases:
        input_seq = ''.join([str(input_seq), 'NN']) if trailing_bases == 1 else ''.join([str(input_seq), 'N'])

    output_seq = Seq.translate(input_seq, to_stop=True)

    if trailing_bases:
        # remove last residue if input needed to be extended because of trailing bases
        output_seq = output_seq[:-1]

    return output_seq


def genes_to_protein():
    '''
    This function translates the ORFs to proteins.
    returns proteins in a list of Seq objects.
    '''
    proteins = []
    for orfs in orf_nucleotides():
        translated = translate_sequence(orfs)
        proteins.append(translated)
    return proteins


#print(genes_to_protein()[-1])


def gc_all_sequences():
    all_gc_cont = []
    for gc in orf_nucleotides():
        all_gc_cont.append(str(round(GC(gc), 2)))
    return all_gc_cont


def create_fasta_protein(d=orf_finder(), e=genes_to_protein(), f=gc_all_sequences()):
    """
    This file creates a FASTA file of protein sequences with header info abut nucleotide sequence
    from which the protein was translated. (start, stop, frame, sense, lenght, GC conent...)
    """

    e_list_str = list(map(str, e))
    d_list_str = list(map(str, d))

    with open("proteins.fasta", 'w') as fp:
        for j in range(len(e_list_str)):
            print("> |NUCLEOTIDE SEQ INFO| GC content: " + f[j] + " %|" + d_list_str[j] + "\n" + e_list_str[j], file=fp)

    return

#create_fasta_protein() #later on in the program.











"""
def export_found_proteins(infile=genes_to_protein(), outfile):
    '''
    This function exports the proteins to a file.
    '''
    for line in infile:
        if line.startswith('>'):
            header = line.lststrip('>')
            outfile.write(header)
        else:

"""

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