import numpy as np


def sequence_composition(k, sequence):
    """
    Given a string `seqeunce`, its k-mer composition Composition(k, sequence)
    is the collection of all k-mer substrings of sequence (including repeated k-mers)

    :return: the collection of all k-mers in LEXICOGRAPHIC order
    """
    res = []

    for i in range(len(sequence) - k + 1):
        res.append(sequence[i:i+k])

    return list(np.sort(res))


def sequence_reconstruction(kmers, string_result=True):
    """
    The inverse of sequence composition... GIVEN a composition of k-mers, we
    need to reconstruct the k-mer composition back into the sequence string
    from which they came.

    To do this, we analyze the suffix of the first
    k-mer (kmner[1:] and look to see if it matches the prefix of the
    second k-mer (kmer[:-1])

    If so, it's added to our genome string.

    :param kmers: a list of k-mers representing a DNA string's k-mer composition
    :return: the DNA string that would have generated the composition
    """


    k = len(kmers[0])
    current_string_suffix = kmers[0][0]  # start with the first character of the first k-mer


    # for k-mers 1 to n-1, (1 being the second), overwrite the current string at i+1 with the first k-1 characters
    i = 1
    for kmer in kmers[1:]:
        next_string_prefix = kmer[0:-1]
        current_string_suffix = current_string_suffix[:i] + next_string_prefix
        i+=1

    current_string_suffix += kmers[-1][-1]  # end with the last character of the last k-mer

    result = current_string_suffix  # To be clear: the final suffix is our result

    return "".join([char for char in result])





def dna_path(dna, return_list=True):

    # Lambda functions to check for overlapping suffixes and prefixes and print them
    check_overlap = lambda pair: pair[0][1:] == pair[1][:-1]
    print_overlap = lambda pair: ' -> '.join(pair)

    # Get all pairs, filter out non-overlapping pairs, then print remaining overlapping pairs
    pairs = ([dna1, dna2] for i, dna1 in enumerate(dna) for j, dna2 in enumerate(dna) if i != j)
    overlaps = map(print_overlap, filter(check_overlap, pairs))

    # Print and save the answers.
    print('\n'.join(overlaps))
    with open('../output_data/genome_seq__reconstruction.txt', 'w') as output_data:
        output_data.write('\n'.join(overlaps))

    # By default, return list to any callers
    if return_list:
        return [overlap for overlap in overlaps]

