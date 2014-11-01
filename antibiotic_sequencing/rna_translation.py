from utils import _GENETIC_CODE

def find_rna_translation(RNA):
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



