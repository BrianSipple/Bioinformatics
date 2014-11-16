_NUCLEOTIDES = ['A', 'C', 'G', 'T']


# DNA transcribes to RNA, which can then be broken
# into the following 3-mer codons.
#
# Each codon (exepct for 3 stop codons) is converted into
# its corresponding amino acid peptide
#
_CODON_PEPTIDE_SYMBOL = {
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
_PEPTIDE_SYMBOL_RNA_ORIGINS = {
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

_ABBREVIATED_PEPTIDE_RNA_ORIGINS = {
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


# For mass spectrometry, the exact measurement is Daltons,
# but having Integer mass conversions get us super-close
_PEPTIDE_TO_MASS = {
    "G": 57,
    "A": 71,
    "S": 87,
    "P": 97,
    "V": 99,
    "T": 101,
    "C": 103,
    "I": 113,
    "L": 113,
    "N": 114,
    "D": 115,
    "K": 128,
    "Q": 128,
    "E": 129,
    "M": 131,
    "H": 137,
    "F": 147,
    "R": 156,
    "Y": 163,
    "W": 186
}

## In some cases, we're only looking for unique masses,
## so we can treat the pairs ("I", "L") and ("K", "Q")
## as interchangeable keys to the same value
_PEPTIDE_TO_UNIQUE_MASS = {
    "G": 57,
    "A": 71,
    "S": 87,
    "P": 97,
    "V": 99,
    "T": 101,
    "C": 103,
    "I": 113,
    #"L": 113,
    "N": 114,
    "D": 115,
    "K": 128,
    #"Q": 128,
    "E": 129,
    "M": 131,
    "H": 137,
    "F": 147,
    "R": 156,
    "Y": 163,
    "W": 186
}

_UNIQUE_MASS_TO_PEPTIDE = {
    57: "G",
    71: "A",
    87: "S",
    97: "P",
    99: "V",
    101: "T",
    103: "C",
    113: "I",
    114: "N",
    115: "D",
    128: "K",
    129: "E",
    131: "M",
    137: "H",
    147: "F",
    156: "R",
    163: "Y",
    186: "W"
}


_AMINO_ACID_SYMBOLS = [
    "G",
    "A",
    "S",
    "P",
    "V",
    "T",
    "C",
    "I",
    "L",
    "N",
    "D",
    "K",
    "Q",
    "E",
    "M",
    "H",
    "F",
    "R",
    "Y",
    "W"
]


def spectrum_to_peptides(spectrum):
    """
    Converts a mass spectrum to its peptide representation
    """
    return [_UNIQUE_MASS_TO_PEPTIDE[mass] for mass in spectrum]



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
    return _PEPTIDE_SYMBOL_RNA_ORIGINS.get(peptide_symbol.upper())


def factorial(n):
    """
    Classic recursive factorial helper

    Why? Why not?
    """
    if n == 0 or n == 1:
        return n
    return n * factorial(n-1)


def nth_triangle(n):
    """
    Like a factorial, but with addition!  (http://en.wikipedia.org/wiki/Triangular_number)
    """
    return (n * (n + 1)) / 2



if __name__ == "__main__":
    pass




