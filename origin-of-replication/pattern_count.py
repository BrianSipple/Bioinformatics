# -*- coding: utf-8 -*-
"""
We will use the term k-mer to refer to a string of 
length k and define Count(Text, Pattern) as the number of times 
that a k-mer Pattern appears as a substring of Text. 

Following the above example,
Count(ACAACTATGCATACTATCGGGAACTATCCT, ACTAT) = 3.

We note that Count(CGATATATCCATAG, ATA) is equal to 3 (not 2) since we should
account for overlapping occurrences of Pattern in Text.

To compute Count(Text, Pattern), our plan is to “slide a window” 
down Text, checking whether each k-mer substring of Text matches Pattern
"""

import numpy as np

genome = "CAATAAAGGAATCAACGGCTCGCGCTGAAACCAATAAACGCTGAAACCGCTGAAACCGCTGAAACCAATAAACAATAAAACGGCTCGCGCTGAAACCGCTGAAACGGAATCAGGAATCACGCTGAAACGGAATCAGGAATCACAATAAAGGAATCAGGAATCAGGAATCACAATAAACCCATAGTGGGAATCACCCATAGTGCAATAAACAATAAACCCATAGTGGGAATCACAATAAACAATAAACGCTGAAACCGCTGAAACCCCATAGTGCAATAAACGCTGAAACCGCTGAAACCCCATAGTGCCCATAGTGCGCTGAAACACGGCTCGGGAATCAACGGCTCGACGGCTCGACGGCTCGCAATAAACAATAAACGCTGAAACCGCTGAAACCAATAAACAATAAACGCTGAAACACGGCTCGGGAATCACGCTGAAACGGAATCACGCTGAAACACGGCTCGCCCATAGTGACGGCTCGCCCATAGTGACGGCTCGGGAATCAACGGCTCGACGGCTCGCAATAAACGCTGAAACACGGCTCGCCCATAGTGCAATAAAGGAATCAGGAATCACAATAAACCCATAGTGCAATAAACGCTGAAACCCCATAGTGCGCTGAAACCCCATAGTGGGAATCAACGGCTCGCCCATAGTGCGCTGAAACCCCATAGTGCCCATAGTGCGCTGAAACCCCATAGTGACGGCTCGGGAATCACAATAAAACGGCTCGCAATAAAGGAATCACGCTGAAACCCCATAGTGCAATAAACCCATAGTGCCCATAGTGCGCTGAAAC"

pattern = "TAATAGATA"



def pattern_count(text, pattern):
	count = 0;
	for i in range(0, (len(text) - len(pattern))):
		if text[i : i + len(pattern)] == pattern:
			count += 1
	return count


"""
But now, how do we find the MOST FREQUENT word, or pattern in a genome?

We say that Pattern is a most frequent k-mer in Text if it maximizes 
Count(Text, Pattern) among all k-mers.

We can say that that ACTAT is a most 
frequent 5-mer of ACAACTATGCATACTATCGGGAACTATCCT, and ATA is 
a most frequent 3-mer of CGATATATCCATAG.

The problem is now well-defined:

Frequent Words Problem: Find the most frequent k-mers in a string.
     Input: A string Text and an integer k.
     Output: All most frequent k-mers in Text.

"""

def most_frequent_patterns(text, word_length):
	"""
	Checks all k-mers (words of a given-length) that appear in 
	a text, then computes HOW MANY TIMES that k-mer appears in the text
	"""
	pattern_frequencies = [] # stores counts for each pattern in the text
	most_frequent_patterns = []

	for i in range(0, len(text) - word_length):
		pattern_to_count = text[i : i + word_length]
		pattern_frequencies.append(pattern_count(text, pattern_to_count))

	max_pattern_freq = np.max(pattern_frequencies)
	
	for i in range(len(pattern_frequencies)):
		if pattern_frequencies[i] == max_pattern_freq:
			pattern = text[i: (i + word_length)]
			if pattern not in most_frequent_patterns:
				most_frequent_patterns.append(pattern)

	return most_frequent_patterns










if __name__ == "__main__":
	print("Example Pattern : {}".format(pattern))
	print ("Count: {}".format(pattern_count(genome, pattern)))

	print("Most frequent 1-mers: {}".format(most_frequent_patterns(genome, 1)))
	print("Most frequent 2-mers: {}".format(most_frequent_patterns(genome, 2)))
	print("Most frequent 3-mers: {}".format(most_frequent_patterns(genome, 3)))
	print("Most frequent 4-mers: {}".format(most_frequent_patterns(genome, 4)))
	print("Most frequent 5-mers: {}".format(most_frequent_patterns(genome, 5)))
	print("Most frequent 6-mers: {}".format(most_frequent_patterns(genome, 6)))
	print("Most frequent 7-mers: {}".format(most_frequent_patterns(genome, 7)))
	print("Most frequent 8-mers: {}".format(most_frequent_patterns(genome, 8)))
	print("Most frequent 9-mers: {}".format(most_frequent_patterns(genome, 9)))
	print("Most frequent 10-mers: {}".format(most_frequent_patterns(genome, 10)))
	print("Most frequent 11-mers: {}".format(most_frequent_patterns(genome, 11)))
	print("Most frequent 12-mers: {}".format(most_frequent_patterns(genome, 12)))








