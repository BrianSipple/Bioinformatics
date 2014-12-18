import unittest
import numpy as np
from genome_sequencing.genome_sequencing import (
    sequence_composition,
    universal_string, debruijn_graph, string_from_path, overlap_graph)
from genome_sequencing.sequence_reconstruction import sequence_reconstruction, eulerian_cycle


class GenomeSequencingTest(unittest.TestCase):

    def setUp(self):
        super(GenomeSequencingTest, self).setUp()


    def tearDown(self):
        super(GenomeSequencingTest, self).tearDown()


    def test_sequence_composition(self):

        k = 5
        sequence = "CAATCCAAC"
        composition = sequence_composition(k, sequence)

        #print(composition)
        self.assertListEqual(composition, [
            "AATCC",
            "ATCCA",
            "CAATC",
            "CCAAC",
            "TCCAA"
        ])

        # k = 100
        # with open('../input_data/genome_sequencing/sequence_composition.txt', 'r') as f:
        #     sequence = f.read()
        #
        # with open('../output_data/genome_sequencing/sequence_composition.txt', 'w') as f:
        #     f.write("\n".join(comp for comp in sequence_composition(k, sequence)))


    def test_string_from_path(self):

        kmers = [
            "AAAT",
            "AATG",
            "ACCC",
            "ACGC",
            "ATAC",
            "ATCA",
            "ATGC",
            "CAAA",
            "CACC",
            "CATA",
            "CATC",
            "CCAG",
            "CCCA",
            "CGCT",
            "CTCA",
            "GCAT",
            "GCTC",
            "TACG",
            "TCAC",
            "TCAT",
            "TGCA"
        ]

        #reconstructed_dna = string_from_path(kmers)
        #print(reconstructed_dna)

        kmers = [
            "ACCGA",
            "CCGAA",
            "CGAAG",
            "GAAGC",
            "AAGCT"
        ]

        reconstructed_dna = string_from_path(kmers)
        #print(reconstructed_dna)
        self.assertEqual("ACCGAAGCT", reconstructed_dna)



    def test_overlap_graph(self):

        kmers = [
            "ATGCG",
            "GCATG",
            "CATGC",
            "AGGCA",
            "GGCAT"
        ]

        overlaps = []

        for overlap in overlap_graph(kmers):
            overlaps.append(overlap)

        overlaps = list(np.sort(overlaps))

        expected_results = [
            "AGGCA -> GGCAT",
            "CATGC -> ATGCG",
            "GCATG -> CATGC",
            "GGCAT -> GCATG"
        ]

        self.assertListEqual(overlaps, expected_results)


    def est_eulerian_path(self):

        with open('../input_data/genome_sequencing/eulerian_path.txt') as input_data:
            edges = {}
            for edge in [line.strip().split(' -> ') for line in input_data.readlines()]:

                if ',' in edge[1]:
                    strings = []
                    for string in edge[1].split(','):
                        strings.append(int(string))

                    edges[int(edge[0])] = strings
                else:
                    edges[int(edge[0])] = [int(edge[1])]

        path = eulerian_caycle(edges)
        #print(path)


    def test_universal_string(self):

        k = 4
        univ_string = universal_string(k)
        self.assertEqual("0000111101100101", univ_string)


    def test_debruijn_graph(self):

        k = 3
        text = "TAATGCCATGGGATGTT"

        # k = 12
        #
        # with open("../input_data/genome_sequencing/debruijn_graph.txt", "r") as f:
        #     text = f.read().strip()

        graph = debruijn_graph(k, text)
        #print(graph)

        with open("../output_data/genome_sequencing/debruijn_graph.txt", 'w') as f:
            f.write("\n".join([item for item in graph]))


    def test_string_reconstruction(self):

        with open('../input_data/genome_sequencing/sequence_reconstruction.txt', 'r') as f:

            read_pairs = []
            for line in f:
                data = line.strip().split(' -> ')
                for i in data[1].split(','):
                    read_pairs.append((data[0], i))

        string = sequence_reconstruction(read_pairs)





if __name__ == "__main__":
    unittest.main()
