from itertools import cycle, dropwhile, islice
from data_structures.dictionaries import FrequencyDict
from utils import _CODON_PEPTIDE_SYMBOL, _PEPTIDE_SYMBOL_RNA_ORIGINS, reverse_complement, dna_to_rna, rna_to_dna, \
    peptide_symbol_to_rna, _PEPTIDE_INTEGER_MASS, _AMINO_ACID_SYMBOLS
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


def count_subpeptides(n):
    """
    Return the number of sub-peptides of a cyclic peptide of length n.

    Being cyclic, a peptide can always form n groups of length 1, 2, 3, ..., n-1. I.e
    the peptide NQEL can form: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, LNQ == 12

    Thus we have our answer... n*n-1
    """
    return n*(n-1)



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


def compute_mass_spectrum(peptide, cyclic=False):
    """
    Given a peptide, we can compute its THEORETICAL SPECTRUM by
    assigning a total integer mass to each of its subpeptides.

    If Peptide represents a cyclic peptide instead, then the
    masses in its theoretical spectrum can be divided into
    those found by LinearSpectrum and those corresponding to
    subpeptides wrapping around the end of Peptide

    NOTE: We also start with 0, and end with the mass of the full
    peptide itself.

    :param peptide: String  - the peptide we aim to smash into subpeptides
    and compute the mass for

    :return: list - integers representing the mass for each cyclopeptide,
    computed in linear form by default
    """
    prefix_masses = [0]

    # We'll start by calculating the mass prefixes
    for index, amino_acid in enumerate(peptide):
        prefix_masses.append(_PEPTIDE_INTEGER_MASS[amino_acid] + prefix_masses[index])

    peptide_mass = prefix_masses[-1]

    # Now, we'll calculate the linear spectrum using the prefix masses
    spectrum = [0]

    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            spectrum.append(prefix_masses[j] - prefix_masses[i])

            #####
            # Cyclic scenario...
            # Each such subpeptide has mass equal to the difference between
            # Mass(Peptide) and a subpeptide mass identified by LinearSpectrum
            ####
            if cyclic and (i > 0 and j < len(peptide)):
                spectrum.append(peptide_mass - (prefix_masses[j] - prefix_masses[i]))

    return sorted(spectrum)



def find_peptide_from_spectrum(spectrum):
    """
    Finds the amino acid peptide that a given spectrum
    must have originated from

    :param spectrum: list - integers corresponding to the masses of a peptide's subpeptides

    :return: the amino acid peptide that originated the spectrum
    """
    pass




def _peptides_with_mass_recursive(mass_table, mass):
    """
    Recursive computation of peptides with a given mass
    """
    if (len(mass_table)) == 0 and mass > 0:
        return 0

    elif mass < 0:
        return 0

    elif mass == 0:
        return 1

    else:

        return _peptides_with_mass_recursive(mass_table, mass - mass_table[0]) + \
               _peptides_with_mass_recursive(mass_table[1:], mass)




def count_peptides_with_mass(mass, recursive=False):
    """
    Given a mass, we compute which peptides produced it
    :param mass - integer - the total integer mass of a peptide

    :return: count - integer - the number of linear peptides having a the mass
    """

    if recursive:
        mass_table = sorted([m for m in list(np.unique(list(_PEPTIDE_INTEGER_MASS.values())))])

        print(mass_table)

        return _peptides_with_mass_recursive(mass_table, mass)

    else:
        raise NotImplementedError("No non-recursive method implemented")



def compute_peptide_total_mass(peptide):
    """
    Computes the total mass of a peptide (the sum of all single subpeptide masses)
    """
    total_mass = 0

    for p in peptide:
        if p in _PEPTIDE_INTEGER_MASS:
            total_mass += _PEPTIDE_INTEGER_MASS[p]

    return total_mass















