_NUCLEOTIDES = ['A', 'C', 'G', 'T']


# DNA transcribes to RNA, which can then be broken
# into the following 3-mer codons.
#
# Each codon (exepct for 3 stop codons) is converted into
# its corresponding amino acid peptide
#
_GENETIC_CODE = {
    "AAA": "K",
    "AAC": "N",
    "AAG": "K",
    "AAU": "N",
    "ACA": "T",
    "ACC": "T",
    "ACG": "T",
    "ACU": "T",
    "AGA": "R",
    "AGC": "S",
    "AGG": "R",
    "AGU": "S",
    "AUA": "I",
    "AUC": "I",
    "AUG": "M",
    "AUU": "I",
    "CAA": "Q",
    "CAC": "H",
    "CAG": "Q",
    "CAU": "H",
    "CCA": "P",
    "CCC": "P",
    "CCG": "P",
    "CCU": "P",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "CGU": "R",
    "CUA": "L",
    "CUC": "L",
    "CUG": "L",
    "CUU": "L",
    "GAA": "E",
    "GAC": "D",
    "GAG": "E",
    "GAU": "D",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCU": "A",
    "GGA": "G",
    "GGC": "G",
    "GGG": "G",
    "GGU": "G",
    "GUA": "V",
    "GUC": "V",
    "GUG": "V",
    "GUU": "V",
    "UAA": "",
    "UAC": "Y",
    "UAG": "",
    "UAU": "Y",
    "UCA": "S",
    "UCC": "S",
    "UCG": "S",
    "UCU": "S",
    "UGA": "",
    "UGC": "C",
    "UGG": "W",
    "UGU": "C",
    "UUA": "L",
    "UUC": "F",
    "UUG": "L",
    "UUU": "F",
}


# We can map each peptide to the RNA codons that convert to it
_PEPTIDE_SYMBOL_ORIGINS = {
    "V": ["GUA", "GUC", "GUG", "GUU"],
    "": ["UGA", "UAG", "UAA"],
    "Y": ["UAC", "UAU"],
    "T": ["ACU", "ACG", "ACC", "ACA"],
    "D": ["GAU", "GAC"],
    "G": ["GGU", "GGA", "GGC", "GGG"],
    "Q": ["CAG", "CAA"],
    "A": ["GCG", "GCA", "GCU", "GCC"],
    "P": ["CCC", "CCU", "CCG", "CCA"],
    "I": ["AUU", "AUA", "AUC"],
    "W": ["UGG"],
    "S": ["UCU", "AGC", "UCA", "UCC", "UCG", "AGU"],
    "E": ["GAA", "GAG"],
    "R": ["AGG", "AGA", "CGU", "CGA", "CGC", "CGG"],
    "M": ["AUG"],
    "C": ["UGC", "UGU"],
    "L": ["CUU", "CUA", "CUC", "UUA", "UUG", "CUG"],
    "F": ["UUC", "UUU"],
    "N": ["AAU", "AAC"],
    "K": ["AAG", "AAA"],
    "H": ["CAU", "CAC"]
}

_ABBREVIATED_PEPTIDE_ORIGINS = {
    "Val": ["GUA", "GUC", "GUG", "GUU"],
    "": ["UGA", "UAG", "UAA"],
    "Tyr": ["UAC", "UAU"],
    "Thr": ["ACU", "ACG", "ACC", "ACA"],
    "Asp": ["GAU", "GAC"],
    "Gly": ["GGU", "GGA", "GGC", "GGG"],
    "Gln": ["CAG", "CAA"],
    "Ala": ["GCG", "GCA", "GCU", "GCC"],
    "Pro": ["CCC", "CCU", "CCG", "CCA"],
    "Ile": ["AUU", "AUA", "AUC"],
    "Trp": ["UGG"],
    "Ser": ["UCU", "AGC", "UCA", "UCC", "UCG", "AGU"],
    "Glu": ["GAA", "GAG"],
    "Arg": ["AGG", "AGA", "CGU", "CGA", "CGC", "CGG"],
    "Met": ["AUG"],
    "Cys": ["UGC", "UGU"],
    "Leu": ["CUU", "CUA", "CUC", "UUA", "UUG", "CUG"],
    "Phe": ["UUC", "UUU"],
    "Asn": ["AAU", "AAC"],
    "Lys": ["AAG", "AAA"],
    "His": ["CAU", "CAC"]
}


def make_spaced_string_from_comma_separated_array(array):
    return " ".join([str(x) for x in array])


def reverse_complement(DNA, as_string=False, unreversed_order=False):
    """
    Returns the reverse complement of a strand of DNA
    """
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    if unreversed_order:
        complements = [comp_dict[nuc] for nuc in DNA]
    else:
        complements = [comp_dict[nuc] for nuc in DNA[::-1]]

    if as_string:
        complements = "".join(complements)

    return complements


def dna_to_rna(DNA):
    """
    Transcribes a String of DNA into RNA
    """
    return "".join([nuc if nuc != 'T' else 'U' for nuc in DNA])


def rna_to_dna(RNA):
    """
    Reverse-transcribes a String of RNA back into DNA
    """
    return "".join([nuc if nuc != 'U' else 'T' for nuc in RNA])

def peptide_symbol_to_rna(peptide_symbol):
    return _PEPTIDE_SYMBOL_ORIGINS.get(peptide_symbol.upper())


if __name__ == "__main__":
    pass


    ##################################################################################
    # Quick helper for generating the _PEPTIDE_ORIGINS dict
    ##################################################################################
    # peptides = {}
    #
    # for codon, peptide in _GENETIC_CODE.items():
    #     if peptide not in peptides:
    #         peptides[peptide] = [codon]
    #     else:
    #         peptides[peptide].append(codon)
    #
    # print(peptides)
    #
    #
    # with open('./translations', 'w') as outfile:
    #
    #     for pep, mapping in peptides.items():
    #         outfile.write('"{}": [{}]\n'.format(
    #             pep,
    #             ", ".join(['"' + codon + '"' for codon in mapping])
    #         ))



