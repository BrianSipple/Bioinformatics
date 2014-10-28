# -*- coding: utf-8 -*-
"""
How do we find the MOST FREQUENT word, or pattern in a genome?

We will use the term k-mer to refer to a string of 
length k.

We say that Pattern is a most frequent k-mer in Text if it maximizes 
Count(Text, Pattern) among all k-mers.

We can say that that ACTAT is a most 
frequent 5-mer of ACAACTATGCATACTATCGGGAACTATCCT, and ATA is 
a most frequent 3-mer of CGATATATCCATAG.

The problem is now well-defined:


Frequent Words Problem: Find the most frequent k-mers in a string.
     Input: A string Text and an integer k.
     Output: All most frequent k-mers in Text.


We aspire to slide a window down Text only once. As we slide this window, 
we will keep track of the number of times that each k-mer Pattern has already appeared in Text, 
updating these numbers as we proceed.

To achieve this goal, we will first order all 4^k k-mers 
lexicographically (i.e., according to how they would appear in the dictionary) and 
then convert them into the 4^k different integers between 0 and 4^k − 1. Given an 
integer k, we define the frequency array of a string Text as an array of 
length 4^k, where the i-th element of the array holds the number of times 
that the i-th k-mer (in the lexicographic order) appears in Text
"""

import numpy as np

try:
    import matplotlib.pyplot as plt
except ImportError:
    pass

import sys
from utils import numberToPattern, patternToNumber
from utils import hammingDistance
from utils import _NUCLEOTIDES
import queue


def pattern_count(text, pattern, start=0, end=0):
    """
     Counts number of overlapping occurrences of `pattern` in `text`.

    :param text: String - Text (DNA) to examine
    :param pattern: String - Pattern (K-mer) to find in text
    :param start: Integer - Start in position <start> in the DNA
    :param end: Integer - End in position <end> in the DNA
    :returns: Integer - n number of occurrences of patter in text
    :raises: ValueError - If start < 0 or >= len(t)
    """
    if start < 0 or start >= len(text):
        raise ValueError("The starting position should be between 0 and the size " + \
                         "of the DNA")

    k = len(pattern)
    count = 0
    end = len(text) - k + 1 if end == 0 else end

    for i in range(0, end):
        if text[i: i + k] == pattern:
            count += 1

    return count


def most_frequent_kmers(DNA, k):
    """
    Returns a list of most frequent k-mers in DNA

    We'll use a Priorty Queue to track each pattern
    along with its frequency in the DNA.

    :param DNA: String - DNA
    :param k: Integer - Length of the K-mer
    :return: Set - Set of most frequent K-mers
    """
    kmers = set()
    most_frequent = queue.PriorityQueue()

    # ###
    # Go through all the K-mers in the DNA string.
    # If the K-mer is found, check whether it has been found previously.
    # If not, record the find and count its occurances.
    ####
    for i in range(len(DNA) - k + 1):

        kmer = ''.join(DNA[i: i + k])

        if kmer in kmers:
            continue
        else:
            kmers.add(kmer)
            ####
            # A priority queue will return the lowest-ordered item first,
            # so we can negate the frequency count it if we want return to the highest-ordered
            # (i.e. highest-frequency k-mer) first
            ####
            most_frequent.put((-pattern_count(DNA, kmer), kmer))

    # Extract ze most-frequent k-mers
    result = set()
    first = True
    freq = -1

    if not most_frequent.empty():

        while not most_frequent.empty():
            kmer = most_frequent.get()

            # We'll need to grab items out of queue until the freq number changes
            if not first and kmer[0] != freq:
                break
            else:
                result.add(kmer[1])
                freq = kmer[0]
                first = False

    return result


def reverse_complement(DNA, as_string=False):
    """
    Returns the reverse complement of a strand of DNA
    """
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complements = [comp_dict[nuc] for nuc in DNA[::-1]]

    if as_string:
        complements = "".join(complements)

    return complements


def find_pattern_positions(pattern, DNA):
    """
    Given a pattern and a larger genome, we want to be
    able to find all starting points of said pattern
    within the genome
    """
    positions = []
    for i in range(len(DNA) - (len(pattern) + 1)):

        current_pattern = DNA[i: (i + len(pattern))]
        if pattern == current_pattern:
            positions.append(i)

    return positions


def find_clumps(DNA, k, L, t):
    """
    Find k-mers forming clumps in DNA

    For a given K-mer, we say it forms an (L, t)-Clump, if K-mer appears at
    least t times in a region of length L in DNA

    :param DNA: String - DNA String
    :param k: Integer - Length of the K-mer
    :param L: Integer - Size of the clump
    :param t: Integer - Number of times kmer must appear in the DNA region

    :return: List - K-mers forming (L, t)-Clumps in DNA
    """
    assert len(DNA) >= L

    clumps = set()

    # Go through all regions of length L in the DNA
    for i in range(0, len(DNA) - L + 1):

        # For each k-mer, find the number of occurrances in the L-length region
        region = DNA[i: i + L]
        kmers = set()
        for j in range(i, i + L - 1 - k):
            kmer = DNA[j: j + k]
            if not kmer in kmers:
                kmers.add(kmer)
                _t = pattern_count(region, kmer)
                if _t >= t:
                    clumps.add(kmer)

        return clumps


def find_highest_frequencies_with_mismatches_by_sorting(text, k, d):
    """
    By knowing that we can have "close" patterns be valid for our search,
    we can find said "close" patterns using the "neighbors" functions, and
    only apply "find_approximate_pattern_matches()" to those as opposed to all k-mers
    """
    num_possibilities = 4 ** k
    most_frequent_patterns = []
    neighborhoods = []
    index = []
    counts = []

    close_patterns = [
                         0] * num_possibilities  # These will be set to 1 whenever k-mer Pattern = NumberToPattern(i, k) is close to some k-mer in Text
    frequency_array = [0] * num_possibilities

    for i in range(len(text) - k + 1):
        current_pattern = text[i: i + k]
        neighborhoods.append(
            neighbors(current_pattern, d))  # Forms an array of arrays holding the neighborhoods of each k-mer

    for neighborhood in neighborhoods:
        index.append(patternToNumber(neighborhood))
        counts.append(1)

    sortedIndex = np.sort(index)

    for i in range(len(neighborhoods)):
        if sortedIndex[i] == sortedIndex[i + 1]:  # repeating values at the index positions indicate a repeating pattern
            counts[i + 1] = counts[i] + 1

    max_count = np.max(counts)

    for i in range(len(neighborhoods)):
        if counts[i] == max_count:
            pattern = numberToPattern(sortedIndex[i], k)
            most_frequent_patterns.append(pattern)

    return most_frequent_patterns


def neighbors(pattern, d):
    """
    Generates the d-neighborhood for a PATTERN.

    In other words: the set of all k-mers
    whose Hamming distance from PATTERN does not exceed `d`

    :return a list of k-mer combinations that comprise the d-neighborhood
    """
    if d == 0:
        return pattern

    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []

    # ##########
    # We can use recursion to successively compute neighbors(suffix(pattern), d),
    # where suffix(pattern) = pattern[1:]
    #
    # The reason being: if we have neighbors(suffix(pattern, d)), then we know
    # that the Hamming Distance between `pattern` and `suffix(pattern)` is either equal
    # to d or less than d.
    #
    # In the first case, we can add `pattern[0]` to the beginning of
    # `suffix(pattern)` in order to obtain a k-mer belonging to
    # Neighbors(Pattern, d). In the second case, we can add any symbol
    # to the beginning of `suffix(pattern)` and obtain a k-mer belonging
    # to Neighbors(Pattern, d).
    # ##########

    suffix_pattern = pattern[1:]
    suffix_neighbors = neighbors(suffix_pattern, d)
    print(suffix_pattern)
    print(suffix_neighbors)

    for i in range(len(suffix_neighbors)):

        neighboring_pattern_text = suffix_neighbors[i]

        if hammingDistance(suffix_pattern, neighboring_pattern_text) < d:
            for n in _NUCLEOTIDES:
                neighborhood.append(n + neighboring_pattern_text)

        else:
            neighborhood.append(pattern[0] + neighboring_pattern_text)

    return neighborhood


def find_approximate_pattern_positions(pattern, genome, d):
    """
    Applying the Hamming Distance Theorem, we'll often want to find
    patterns in a genome having d or fewer mistmatches with a given `pattern`.

    Approximate Pattern Matching Problem:
        - Find all approximate occurrences of a pattern in a string.

     Input:
         - Strings Pattern and Text along with an integer d.

     Output:
         - All starting positions where Pattern appears as a substring
           of Text with at most d mismatches.
    """
    positions = []
    pattern_length = len(pattern)

    found_patterns = find_highest_frequencies_with_mismatches_by_sorting(genome, pattern_length, d)

    for pat in found_patterns:
        positions.append(patternToNumber(pat))

    return positions


"""
Although the frequencies of A and T are practically identical
on the two half-strands, C is more frequent on the reverse half-strand
than on the forward half-strand.

It turns out that we observe these discrepancies because cytosine (C)
has a tendency to mutate into thymine (T) through a process called DEANIMATION.
Deamination rates rise 100-fold when DNA is single-stranded, which leads
to a decrease in cytosine on the forward half-strand. Also,
since C-G base pairs eventually change into T-A base pairs,
deamination results in the observed decrease in guanine (G) on the reverse half-strand.

This means that we can compute the G-C skew at every point in a genome...
"""


def compute_gc_skew(DNA, chart=False):
    """
    Give all values of Skew(DNA)

    We define skew(DNA) as the difference between the total number of occurrences
    of G and the total number of occurrences of C in DNA.

    :param DNA: String - DNA to calculate skew
    :param chart: Boolean - If True, will save a skew.png chart in the current directory.

    :return: Tuple - With skew array and list of positions where the skew minimizes.
    """
    res = [0]
    G_C = 0
    _min = 0
    indexes = []

    for i in range(len(DNA)):

        if DNA[i] == "G":
            G_C += 1

        elif DNA[i] == "C":
            G_C -= 1

        # Compute the min. NOTE: We have to sum one to the indexes because we already
        # start with an extra element in the res (a 0)
        if G_C < _min:
            indexes = [i + 1]
        elif G_C == min:
            indexes.append(i + 1)
        res.append(G_C)

    if chart:
        if sys.modules.get('matplotlib', None):
            plt.plot(res)
            plt.ylabel('G - C diff')
            plt.title('Skew diagram')
            plt.savefig('skew.png')
        else:
            print("No matplotlib module found -- no skew diagram for you :-(")

    return (res, indexes)


"""
Having computed the skews, we can find oriC!...

We know that the skew is decreasing along the reverse half-strand
and increasing along the forward half-strand. Thus, the skew should
achieve a minimum at the position where the reverse half-strand ends
and the forward half-strand begins... which is exactly the location of oriC!

This insight can now be used as its own oriC detection algorithm: finding the location
of the minimum G-C skew
"""

def find_min_gc_skew_locations(DNA):
    """
    Input: A DNA string Genome.
    Output: All integer(s) i minimizing Skew[i] among all values of i (from 0 to |Genome|)
     """
    min_skew_locations = []

    skew_values = compute_gc_skew(DNA)

    min_skew_number = np.min(skew_values)

    for i in range(len(skew_values)):
        if skew_values[i] == min_skew_number:
            min_skew_locations.append(i)

    return min_skew_locations


if __name__ == "__main__":
    pass
    # print("Example Pattern : {}".format(pattern))
    # print ("Count: {}".format(pattern_count(genome, pattern)))


    # with open('resources/oric_vibrio_cholerae.txt', 'rb') as file:
    # oric_vibrio_cholerae = file.read()
    #
    # with open('resources/oric_thermotoga_petrophila.txt', 'rb') as file:
    # oric_thermotoga_petrophila = file.read()
    #
    #     # k = 8
    #
    # #print("Most frequent 1-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 1)))
    # # print("Most frequent 2-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 2)))
    # #print("Most frequent 3-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 3)))
    # # print("Most frequent 4-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 4)))
    # # print("Most frequent 5-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 5)))
    # # print("Most frequent 6-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 6)))
    # # print("Most frequent 7-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 7)))
    # # print("Most frequent 8-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 8)))
    # #print("Most frequent 9-mers: {}".format(computeFrequenciesBySorting(oric_thermotoga_petrophila, 9)))
    # # print("Most frequent 10-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 10)))
    # # print("Most frequent 11-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 11)))
    # # print("Most frequent 12-mers: {}".format(computeFrequenciesBySorting(oric_vibrio_cholerae, 12)))
