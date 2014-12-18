from collections import OrderedDict
from itertools import product
from collections import defaultdict
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



def string_from_path(kmers):

    result = kmers[0]

    for kmer in kmers[1:]:
        result += kmer[-1]

    return result





def overlap_graph(kmers):
    """
    Constructs an overlap graph (In the form of an adjacency list)
    for a list of k-mers.

    I.E: If we treat each k-mer in a set as a node, we can map edges to
    each of the other k-mers where the k-1-mer prefix of kmer and k-1k-mer
    suffix of kmer' match.
    """
    prefix = {}
    for i in kmers:
        p = i[:-1]
        if p in prefix:
            prefix.get(p).append(i)
        else:
            prefix.update({p: [i]})

    for i in kmers:
        s = i[1:]
        if s in prefix:
            for k in prefix.get(s):
                yield '%s -> %s' % (i, k)



def universal_string(k):
    """
    Contructs a k-univeral binary string_result
    """
    kmers = [''.join(x) for x in list(product('01', repeat=k - 1))]

    edges = []
    for i in overlap_graph(kmers):
        data = i.split(' -> ')
        for i in data[1].split(','):
            edges.append((data[0], i))

    small, large = classify_edges(edges)
    origin = '0' * (k - 1)
    start = large.get(origin)
    cycle_edges = [(origin, small.get(origin)), (origin, start)]
    while 1:
        if start in large:
            cycle_edges.append((start, large.get(start)))
            end = large.get(start)
            large.pop(start)
            start = end
        elif start in small:
            cycle_edges.append((start, small.get(start)))
            end = small.get(start)
            small.pop(start)
            start = end
        else:
            break

    return glue(cycle_edges)



def glue(edges):
    k = len(edges[0][0])
    string = [edges[0][0] + edges[0][1][-1]]

    for i, j in edges[1:]:
        string.append(j[-1])

    return ''.join(string[:-(k + 1)])


def classify_edges(edges):
    small_edges = {}
    large_edges = {}

    for i, j in edges:
        if i in small_edges:
            large_edges.update({i: j})
        else:
            small_edges.update({i: j})

    return small_edges, large_edges



def debruijn_graph(k, sequence):
    """
    Constructs a graph of length |text| - k + 1 in the form
    of an adjacency list
    """
    kmers = set(gen_kmers(k, sequence))
    return overlap_graph(kmers)


def gen_kmers(k, sequence):
    for i in range(len(sequence)):
        yield sequence[i: i + k]



