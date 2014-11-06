from itertools import cycle, dropwhile, islice
from data_structures.dictionaries import FrequencyDict
from utils import _CODON_PEPTIDE_SYMBOL, _PEPTIDE_SYMBOL_RNA_ORIGINS, reverse_complement, dna_to_rna, rna_to_dna, \
    peptide_symbol_to_rna, _PEPTIDE_INTEGER_MASS
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


def compute_cyclopeptides(peptide):
    """
    Models a mass spectrometer that breaks copies of a cyclic peptide
    at every possible two bonds, so that the resulting experimental specturm
    contains the masses of all possible linear fragments of the peptide, called "subpeptides"

    Example: the cyclic peptide NQEL has 12 subpeptides: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, and LNQ

    :return: the list of subpeptides
    """
    res = []

    amino_acids_to_grab = 1
    cycles_remaining = len(peptide) - 1

    while cycles_remaining > 0:

        cycled = cycle([amino_acid for amino_acid in peptide])

        for i in range(len(peptide)):

            skipped = dropwhile(lambda x: x != peptide[i], cycled)

            # Grab the amino slice from the iterable of amino acids, starting from the first not
            # skipped and extending to the length of the slice we need to take
            amino_slice = islice(skipped, amino_acids_to_grab)

            res.append("".join(amino_slice))
        amino_acids_to_grab += 1
        cycles_remaining -= 1

    return res


def compute_mass_spectrum(cyclopeptide, sorted=True):
    """
    Given a cyclopeptides, we can compute its THEORETICAL SPECTRUM by
    assigning a total integer mass to each of its subpeptides

    NOTE: We also start with 0, and end with the mass of the full
    peptide itself.

    :param cyclopeptide: String  - the peptide we aim to smash into subpeptides
    and compute the mass for

    :return: list - integers representing the mass for each cyclopeptide, sorted ascending by default
    """
    masses = [0]

    subpeptides = compute_cyclopeptides(cyclopeptide)
    print(subpeptides)

    for subpep in subpeptides:

        mass = 0
        for amino_acid in subpep:
            mass += _PEPTIDE_INTEGER_MASS[amino_acid]

        masses.append(mass)

    # One final mass computation for the entire cyclopeptide
    mass = 0
    for amino_acid in cyclopeptide:
        mass += _PEPTIDE_INTEGER_MASS[amino_acid]
    masses.append(mass)

    if sorted:
        masses = list(np.sort(masses))

    return masses


















