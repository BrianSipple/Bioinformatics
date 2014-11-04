from data_structures.dictionaries import FrequencyDict
from utils import _CODON_PEPTIDE_SYMBOL, _PEPTIDE_SYMBOL_RNA_ORIGINS, reverse_complement, dna_to_rna, rna_to_dna, \
    peptide_symbol_to_rna
import numpy as np


def translate_to_peptides(RNA):
    """
    Find the translation of an RNA string into an amino acid

    :param RNA: String - the pattern of nulcleotides
    :param genetic_code: Array - the 64-element array of codons that translate
    into amino acids

    :return: peptide: String - the final animo acids that RNA has been coded to
    """
    assert len(RNA) >= 3
    assert len(RNA) % 3 == 0

    peptides = ""

    for i in range(0, len(RNA) - 3 + 1, 3):
        codon = RNA[i: i + 3]
        if codon in _CODON_PEPTIDE_SYMBOL:
            peptides += _CODON_PEPTIDE_SYMBOL[codon]

    return peptides


def compute_possible_rna_origins(peptide_list):
    """
    Compute the possible DNA origins for a list of amino acid peptides

    Ideally, the list would represent a larger organism made up of multiple
    peptides.

    :param peptide_list: Array - a list of peptides

    :return: count - The number of possible DNA origins
    :return: possible_origins: Array - the list of possible origins... if desired
    """
    possible_origins = []

    for peptide in peptide_list:
        if peptide in _PEPTIDE_SYMBOL_RNA_ORIGINS:
            for origin in _PEPTIDE_SYMBOL_RNA_ORIGINS[peptide]:
                possible_origins.append(origin)

    return possible_origins


def recursive_find_rna_encoders(encoders, amino_seq):
    """
    Recursively generate all RNA sequences that encode a given amino acid sequence
    """
    if not encoders:
        [encoders.add(t) for t in peptide_symbol_to_rna(amino_seq[0])]
        encoders = recursive_find_rna_encoders(encoders, amino_seq[1:])

    elif amino_seq:
        aux = set()
        for t in peptide_symbol_to_rna(amino_seq[0]):
            for e in encoders:
                e += t
                aux.add(e)
        encoders = recursive_find_rna_encoders(aux, amino_seq[1:])

    return encoders



def compute_possible_dna_origins(DNA, final_peptides):
    """
    Given a string of DNA and a string of peptides, find the
    subsets of DNA withing the String that could have encoded the peptide

    (Encoding process: DNA -> RNA -> peptide)

    NOTE: for each String of DNA, we have to get the reverse complement that
    we know will be attached to it during transcription into RNA

    :param DNA: String - The strand of DNA
    :param peptides: String - The peptide produced after transcription and translation

    :return: origins: Array - A list of possible DNA origins for the peptides
    """
    encoders = recursive_find_rna_encoders(set(), final_peptides)
    freq_dict = FrequencyDict(DNA, len(final_peptides) * 3)

    res = []

    for codon in encoders:

        enc_dna = rna_to_dna(codon)
        enc_rev = reverse_complement(enc_dna, as_string=True)

        freq = freq_dict.get(enc_dna, 0)
        res.extend([enc_dna] * freq)

        freq = freq_dict.get(enc_rev, 0)
        res.extend([enc_rev] * freq)

    return res















