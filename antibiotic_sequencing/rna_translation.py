from utils import _GENETIC_CODE, _PEPTIDE_SYMBOL_ORIGINS, reverse_complement, dna_to_rna, rna_to_dna
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
        if codon in _GENETIC_CODE:
            peptides += _GENETIC_CODE[codon]

    return peptides


def compute_possible_rna_origins(peptide_list, return_list=False):
    """
    Compute the possible DNA origins for a list of amino acid peptides

    Ideally, the list would represent a larger organism made up of multiple
    peptides.

    :param peptide_list: Array - a list of peptides

    :return: count - The number of possible DNA origins
    :return: possible_origins: Array - the list of possible origins... if desired
    """
    count = 0
    possible_origins = []

    for peptide in peptide_list:
        if peptide in _PEPTIDE_SYMBOL_ORIGINS:
            count += len(_PEPTIDE_SYMBOL_ORIGINS[peptide])
            if return_list:
                for origin in _PEPTIDE_SYMBOL_ORIGINS[peptide]:
                    possible_origins.append(origin)

    if return_list:
        return count, possible_origins

    return count


def compute_possible_dna_origins(DNA, final_peptide):
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
    possible_origins = []
    reverse_comp = reverse_complement(DNA, as_string=True, unreversed_order=True)
    rna = dna_to_rna(DNA)
    rna_complement = dna_to_rna(reverse_comp)
    final_peptide_length = len(final_peptide)

    for i in range(0, len(rna) - (3 * final_peptide_length) + 1, 3):

        corresponding_codons = rna[ i : i + (3 * final_peptide_length) ]
        codon_peptides = translate_to_peptides(corresponding_codons)

        if codon_peptides == final_peptide:
            # We have a match!
            possible_origins.append(rna_to_dna(corresponding_codons))

    for i in range(0, len(rna_complement) - (3 * final_peptide_length) + 1, 3):

        corresponding_codons = rna_complement[ i : i + (3 * final_peptide_length)]
        codon_peptides = translate_to_peptides(corresponding_codons)

        if codon_peptides == reversed(final_peptide):
            possible_origins.append(rna_to_dna(corresponding_codons))

    return possible_origins










