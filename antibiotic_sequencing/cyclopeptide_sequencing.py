#!/usr/bin/env python
from itertools import cycle, dropwhile, islice
from utils import nth_triangle, _PEPTIDE_TO_MASS, _PEPTIDE_TO_UNIQUE_MASS, spectrum_to_peptides


def main():
    with open("resources/convolution_cyclopeptide_sequencing.txt", 'r') as fi:
        m = int(fi.readline())
        n = int(fi.readline())
        spectrum = [int(i) for i in fi.readline().strip().split(' ')]

    top_elements = convolution(spectrum, m)
    print(top_elements)
    leadercyclopeptide = leaderboard_cyclopeptide_sequencing(spectrum, n, top_elements)

    if leadercyclopeptide:
        print('-'.join([str(i) for i in leadercyclopeptide]))





def count_subpeptides_in_cyclopeptides(pep_length):
    """
    Return the number of sub-peptides of a cyclic peptide of length n.

    Being cyclic, a peptide can always form n groups of length 1, 2, 3, ..., n-1. I.e
    the peptide NQEL can form: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, LNQ == 12

    Essentially... n * n-1... plus 2 for factoring in the first empty string and last whole number
    """
    return pep_length * (pep_length - 1) + 2


def count_subpeptides_in_linear_peptide(pep_length):
    """
    The number of linear subpeptides in a peptide is the
    nth_triangle number of the length  (http://en.wikipedia.org/wiki/Triangular_number)

    ...plus 2 for factoring in the first empty string and last whole number
    """
    return nth_triangle(pep_length) + 2


def compute_cyclopeptides(peptide):
    """
    Models a mass spectrometer that breaks copies of a cyclic peptide
    at every possible two bonds, so that the resulting experimental specturm
    contains the masses of all possible linear fragments of the peptide, called "subpeptides"

    Example: the cyclic peptide NQEL has 12 subpeptides: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, and LNQ

    :return: the list of subpeptides
    """
    if len(peptide) == 1:
        return peptide

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
        prefix_masses.append(_PEPTIDE_TO_MASS[amino_acid] + prefix_masses[index])

    peptide_mass = prefix_masses[-1]

    # Now, we'll calculate the linear spectrum using the prefix masses
    spectrum = [0]

    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            spectrum.append(prefix_masses[j] - prefix_masses[i])

            # ####
            # Cyclic scenario...
            # Each such subpeptide has mass equal to the difference between
            # Mass(Peptide) and a subpeptide mass identified by LinearSpectrum
            # ###
            if cyclic and (i > 0 and j < len(peptide)):
                spectrum.append(peptide_mass - (prefix_masses[j] - prefix_masses[i]))

    return sorted(spectrum)




def compute_peptide_total_mass(peptide):
    """
    Computes the total mass of a peptide (the sum of all single subpeptide masses)
    """
    total_mass = 0

    for p in peptide:
        if p in _PEPTIDE_TO_MASS:
            total_mass += _PEPTIDE_TO_MASS[p]

    return total_mass


def find_cyclopeptide_in_mass_spectrum(exp_spectrum):
    """
    Iteratively build a list of candidate peptides for a total mass, where the final
    candidates will be those whose theoretical spectra are “consistent” with the experimental
    spectrum.

    We start with an empty 0-mer string ("") for a candidate, and a mass of 0, and at each
    iteration we'll add 18 new k+1-mer candidates to the list.

    In the same loop, we'll trim the list, only keeping the
    linear peptides that remain consistent with the experimental spectrum.

    :param mass: the total mass of all codons in a peptide
    :return: the cyclopeptides
    """
    results = []
    candidate_spectrums = [[]]
    parent_mass = int(exp_spectrum[-1])  # treat the largest mass of the spectrum as the parent mass

    while len(candidate_spectrums) > 0:

        candidate_spectrums = expand(candidate_spectrums)

        # Make a copy so we can delete from the original
        spectrums = candidate_spectrums[:]

        #print("Candidate spectrums: {}".format(candidate_spectrums))

        for linear_spec in spectrums:

            linear_candidate_mass = sum(linear_spec)

            #print("Linear candidate mass: {}".format(linear_candidate_mass))

            if linear_candidate_mass == parent_mass:  # potential winner
                #print("Mass of {} is equal to parent mass! Checking for theoretical cyclospectrum match".format(
                #    linear_spec))

                candidate_peptides = spectrum_to_peptides(linear_spec)
                theor_spec = compute_mass_spectrum(candidate_peptides, cyclic=True)

                if theor_spec == exp_spectrum:
                    results.append(linear_spec)
                else:
                    # Remove -- it's not consistent
                    candidate_spectrums.remove(linear_spec)

            elif linear_candidate_mass > parent_mass:
                # We've crossed the limit... something went wrong
                candidate_spectrums.remove(linear_spec)


            # If we're here, we still need to remove any "inconsistent" new candidates to prep
            # for the next iteration
            elif not is_linear_spectrum_consistent(linear_spec, exp_spectrum):
                candidate_spectrums.remove(linear_spec)

    return results


def expand(spectrums):
    """
    For every spectrum... in a list of spectrums... we create 18
    new k+1-mer spectrum representations

    :return: a list of lists containing spectrums
    """

    expanded = []

    for spectrum in spectrums:
        for peptide in _PEPTIDE_TO_UNIQUE_MASS:
            new_spectrum = spectrum + [_PEPTIDE_TO_UNIQUE_MASS[peptide]]
            expanded.append(new_spectrum)

    return expanded


def is_linear_spectrum_consistent(linear_spec, exp_spec):
    """
    checks whether the linear spectrum of a candidate peptide is "consistent" with
    the cyclic experimental spectrum we're working with.

    We define consistency as

    :return: consistent - boolean
    """

    for subpeptide in linear_spec:
        if subpeptide not in exp_spec:
            return False

    return True


def spectrum_score(peptide, exp_spec):
    """
	Returns the number of matching masses from a peptide's linear spectrum
	when compared with a given experimental_spectrum
    """
    linear_spectrum = compute_mass_spectrum(peptide, cyclic=False)
    # Return -1 is the linear spectrum has more mass than the experimental spectrum
    if linear_spectrum[-1] > exp_spec[-1]:
        return -1

    return sum([min(linear_spectrum.count(mass), exp_spec.count(mass)) for mass in set(linear_spectrum)])




def convolution(spectrum, m):
    spectrum = sorted(spectrum)

    elements = {}
    for i in spectrum:
        for j in spectrum:
            if i == j:
                break
            else:
                difference = abs(i - j)
                if 57 <= difference <= 200:
                    if difference in elements:
                        elements[difference] += 1
                    else:
                        elements.update({difference: 1})

    n = 0
    last_score = None
    top_elements = []
    for i in sorted(elements, key=lambda x: elements.get(x), reverse=True):
        n += 1
        score = elements.get(i)
        if n > m and last_score and score == last_score:
            top_elements.append(i)
        elif n > m:
            break
        else:
            top_elements.append(i)
        last_score = score

    return top_elements


def expand(spectrum=None, limit_spec=None):

    if limit_spec:
        new_mass = []
        for i in _PEPTIDE_TO_UNIQUE_MASS:
            if i in limit_spec:
                new_mass.append(i)
        mass = new_mass
    if spectrum:
        spectrum_ext = []
        for i in spectrum:
            for j in _PEPTIDE_TO_UNIQUE_MASS:
                spectrum_ext.append(i + (j,))
        return spectrum_ext
    else:
        return [(i,) for i in _PEPTIDE_TO_UNIQUE_MASS]


def cyclospectrum(spectrum):
    """Return the spectra of the cyclopeptide"""
    spec = spectrum * 2

    for i in range(len(spectrum)):
        for j in range(i, i + len(spectrum)):
            yield(spec[i:j])


def cut(leaderboard, spectrum, n):
    scores = {}
    for peptide in leaderboard:
        if peptide:
            scores.update({peptide: spectrum_score(peptide, spectrum)})
    if not scores:
        return []
    rank = 1
    new_peptides = []
    ranked_peptides = sorted(scores, key=lambda x: scores.get(x), reverse=True)
    prepep = ranked_peptides[0]
    new_peptides.append(prepep)
    for i in ranked_peptides[1:]:
        if rank > n and scores.get(i) == scores.get(prepep):
            new_peptides.append(i)
        elif rank > n:
            break
        else:
            new_peptides.append(i)
            rank += 1
            prepep = i

    return new_peptides


def leaderboard_cyclopeptide_sequencing(spectrum, n, limit_spec=None):
    mass = [
        57,
        71,
        87,
        97,
        99,
        101,
        103,
        113,
        114,
        115,
        128,
        129,
        131,
        137,
        147,
        156,
        163,
        186,
    ]
    leaderboard = [(i,) for i in mass]
    if limit_spec:
        new_leaderboard = []
        for i in leaderboard:
            if i[0] in limit_spec:
                new_leaderboard.append(i)
        leaderboard = new_leaderboard
    leaderpeptide = None
    max_mass = max(spectrum)
    while leaderboard:
        leaderboard = expand(leaderboard, limit_spec)
        for i, j in enumerate(leaderboard):
            if sum(j) == max_mass:
                if leaderpeptide:
                    if spectrum_score(j, spectrum) > spectrum_score(leaderpeptide, spectrum):
                        leaderpeptide = j
                else:
                    leaderpeptide = j
            elif sum(j) > max_mass:
                leaderboard[i] = 0
        leaderboard = cut(leaderboard, spectrum, n)

    return leaderpeptide








if __name__ == '__main__':
    main()