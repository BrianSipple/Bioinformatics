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


_SYMBOLS = ['A', 'C', 'G', 'T']


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



def frequentPatterns(text, k):
	"""
	Fast algorithm for finding the most frequenct words or patterns
	of k-length (referred to as "k-mers") within a string of text.

	The steps are as follows: (Also, see additional helper-method documentation)
	
	1) Generate a frequency array for every k-mer in the text.
	
	2) Find all most frequent k-mers by finding all k-mers corresponding
	   to the maximum element(s) in the frequency array.

   	3) Return the k-mers as a list.
	"""

	frequency_array = computeFrequencies(text, k)

	max_frequency = np.max(frequency_array)

	indicies_of_most_freqrent_words = [i for i, freq in enumerate(frequency_array)
		if freq == max_frequency]

	most_frequent_patterns = [numberToPattern(i, k) for i in indicies_of_most_freqrent_words]

	return most_frequent_patterns






def computeFrequenciesBySorting(text, k):
	"""
	This performs the setup of our nucleotide sequence by 
	initializing every element of the frequency array to 0 (k^4 operations)
	and then making a single pass down Text (approximately |Text| · k operations). 
	
	For each k-mer Pattern that we encounter, we add 1 to the value of the array 
	corresponding to Pattern. As before, we refer to the k-mer beginning at 
	position i of Text as Text(i, k).
	"""
	
	pattern_size = (4**k)

	unique_patterns = []  # Keep track of the unique patterns we find
	most_frequent_patterns = []
	frequency_counts = []  # computes the number of times that the integer at position i in the array SortedIndex appears in the first i elements of this array
	indicies = []

	for i in range( (len(text) - k) + 1):
		
		current_pattern = text[i:(i+k)]
		if current_pattern not in unique_patterns:
			unique_patterns.append(current_pattern) 

		
		# if current_pattern not in patterns:
		# 	patterns.append(current_pattern)		

		indicies.append(patternToNumber(current_pattern))  # This gives us the location within the array of counts!
		
		frequency_counts.append(1)

	sortedIndicies = list(np.sort(indicies))

	for i in range(1, len(sortedIndicies)):
		if sortedIndicies[i] == sortedIndicies[i-1]:  # Repeating patterns will have created repeating index numbers
			frequency_counts[i] +=1

	max_count = np.max(frequency_counts)

	for i in range(len(sortedIndicies)):
		if frequency_counts[i] == max_count:
			most_frequent_patterns.append(numberToPattern(sortedIndicies[i], k))

	return most_frequent_patterns



def patternToNumber(pattern):
	"""
	Transforms a k-mer Pattern into an integer to mark its index location.
	
	Our set of locations will range from 0 to ((4^k)-1)
	"""
	if len(pattern) == 0:
		return 0
	symbol = pattern[-1]
	pattern = pattern[:-1]
	return (4 * patternToNumber(pattern)) + symbolToNumber(symbol)



def symbolToNumber(symbol):
	"""
	"""
	if symbol in _SYMBOLS:
		return _SYMBOLS.index(symbol)


def numberToSymbol(number):
	return _SYMBOLS[int(number % 4)]



def numberToPattern(index, k):
	"""
	Reverses patternToNumber, transforming an integer between 0 and 4^k − 1 into a k-mer

	Steps:
		1) Divide the index # by 4 to obtain a QUOTIENT and a REMAINDER
		2) The REMAINDER, r, represents the final nucleotide of the pattern. 
			- e.g: NumberToSymbol(r)
		3) Recurse...
			- divide each subsequent quotient by 4, until we obtain a quotient of 0.
			- PUSH each symbol onto the FRONT of the string
		4) Return the final nucleotide pattern
	"""
	if k == 1:
		return numberToSymbol(index)

	prefixIndex = quotient(index, 4)

	r = remainder(index, 4)

	prefix_pattern = numberToPattern(prefixIndex, k-1)
	symbol = numberToSymbol(r)

	return prefix_pattern + symbol



def quotient(n, m):
	return int(n / m)

def remainder(n, m):
	return n % m


if __name__ == "__main__":
	# print("Example Pattern : {}".format(pattern))
	# print ("Count: {}".format(pattern_count(genome, pattern)))

	genome = "atcaatgatcaacgtaagcttctaagcATGATCAAGgtgctcacacagtttatccacaacctgagtggatgacatcaagataggtcgttgtatctccttcctctcgtactctcatgaccacggaaagATGATCAAGagaggatgatttcttggccatatcgcaatgaatacttgtgacttgtgcttccaattgacatcttcagcgccatattgcgctggccaaggtgacggagcgggattacgaaagcatgatcatggctgttgttctgtttatcttgttttgactgagacttgttaggatagacggtttttcatcactgactagccaaagccttactctgcctgacatcgaccgtaaattgataatgaatttacatgcttccgcgacgatttacctcttgatcatcgatccgattgaagatcttcaattgttaattctcttgcctcgactcatagccatgatgagctcttgatcatgtttccttaaccctctattttttacggaagaATGATCAAGctgctgctcttgatcatcgtttc".upper()
	k = 8

	print("Most frequent 1-mers: {}".format(computeFrequenciesBySorting(genome, 1)))
	print("Most frequent 2-mers: {}".format(computeFrequenciesBySorting(genome, 2)))
	print("Most frequent 3-mers: {}".format(computeFrequenciesBySorting(genome, 3)))
	print("Most frequent 4-mers: {}".format(computeFrequenciesBySorting(genome, 4)))
	print("Most frequent 5-mers: {}".format(computeFrequenciesBySorting(genome, 5)))
	print("Most frequent 6-mers: {}".format(computeFrequenciesBySorting(genome, 6)))
	print("Most frequent 7-mers: {}".format(computeFrequenciesBySorting(genome, 7)))
	print("Most frequent 8-mers: {}".format(computeFrequenciesBySorting(genome, 8)))
	print("Most frequent 9-mers: {}".format(computeFrequenciesBySorting(genome, 9)))
	print("Most frequent 10-mers: {}".format(computeFrequenciesBySorting(genome, 10)))
	print("Most frequent 11-mers: {}".format(computeFrequenciesBySorting(genome, 11)))
	print("Most frequent 12-mers: {}".format(computeFrequenciesBySorting(genome, 12)))

	pattern = "ATGCAA"
	print("Pattern to Number of " + pattern + ":")
	print(patternToNumber(pattern))

	###### NumberToPattern ######
	number = 5437
	k = 8

	print("Number to Pattern for Pattern: {}, k: {}".format(number, k))
	print(numberToPattern(number, k))

	k = 2
	text = "AAGCAAAGGTGGG"

	print("Frequncy counts for Text {}\n and k = {}".format(text, k))
	print(computeFrequenciesBySorting(text, k))

	array = computeFrequenciesBySorting(text, k)
	string = " ".join(str(i) for i in array)









