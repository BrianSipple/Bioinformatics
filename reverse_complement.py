"""
Given a nucleotide p, we denote its complementary nucleotide as p. 
The reverse complement of a string Pattern = p1…pn is the string 
Pattern = pn … p1 formed by taking the complement of each nucleotide 
in Pattern, then reversing the resulting string. 
We will need the solution to the following problem throughout the book:

Reverse Complement Problem: Find the reverse complement of a DNA string.
     Input: A DNA string Pattern.
     Output: Pattern, the reverse complement of Pattern.

Sample Input:
     AAAACCCGGT

Sample Output:
     ACCGGGTTTT     
"""

complements = {
	'A': 'T',
	'T': 'A',
	'G': 'C',
	'C': 'G'
}

def reverseComplement(pattern):
	reverseComplement = ""
	for i in range( (len(pattern)-1), -1, -1):
		reverseComplement += complements.get(pattern[i], '').upper()

	return reverseComplement



if __name__ == "__main__":

	genome = "GGATTCAAGCTAGAAGTGTGCACATGGTAATATCTCATGGATGTAACGCCGCAAGACTGTTGAGTCTCCAGAAATGGGAGTAAGGCGGCTATCACAGGTCGGGTTTTGTACAGCATTCAATCGAGGAATGTACTCCCAAATTACCACACTAAAACTCGGACAGGGAGGTGACCTGTAATCGGGCAGTACCATACTCGTCCCACTGCATAAGAGTTTCAGTAAAGTGTCTTAAACATCCGTTCGCCGACGAGCACAAGGTGAATAGGCCGGGGTTTGTTTAACACCTGGGCGCTGCCGCCAGGACGCTCCGCCTAAACTTAGTCCCCGTGTTGATTACAGTGGTGCAATATACTAGCATGCCATTCCATTTACCGGTCAGGAAGGAACGTAATAATCCCCGATTTCGAGTAGATCCTGGAGATGGGGATCTACGGATAGCCTACCGGTACTACGATGATAGCTCTTCCTCACCTCGTCAACGAACCACTAGTGGCCAAGAATATTATGATTAACCCTTAACCATAGGTGACATGCACTCCGTTAGCTTGCGTCATCCGAAGAGATCCCGCGGTAGTATGAAAGGGCCTCAGCTAGTATTCCGTGGAGGTGTGGCCGTGCCAGCGATATGTAGAGCCTACGCGGTGCCGCGCCAGTGAGTGACCGCCTCGCCGATCTGGCCCCCAGACATCACAAGATATGGCTGTGTAGGTCTTTCCATGCTGTAGTTCCGTGAGCCTCAGTTACACCTGTAAGCTCGTGAAACTGTCTACTATTCCCCAAATGCACTCACATATGATTTAATTTGAACGTAGTGAAATGGAGGTTGTTTTGCTTCAATACACTCAACTTAATCACCTTAGTTTACACGACGGTTGTATGTACGGACACTAATACCATCCTGTTGTGCCAGCGATTCCGCTGTCTGTCCTCACGGGTAATCTATCTCACAGCCTCGATTACAGGTTTCGCCCTAGCCATGGTGTCCGCCCTTGGCGCGTTTGTGATAGGTTTCAGTCATAGTGGGGTGACATGGCGCACTTCAACAAAATTGGAGTACGATTGCCACGGAATATTGGTATATCGGCATACCTGAGCACGACGGTTTCGTCCCACTACATAGTCTGAGATCCCGGGCGAGAGCAATATTCGGGGGTGATTACTAGCTAATGCACTATGTTTTTACTTTGCTGCTGGGTATGAGGAGACCAAATTGGTCTCTACTCTTGTCGCTCTGTCTCCAGAGGACTAGCCTCGGCGACGAGACTTAGCAGAGATGCACATATGTTAGATAAGGCTTCCGTACAAGCAGCGTGGTCCGCGTCCGCTCCCGTTAGTCGACTTGTCCCACTGTCCTCTATGTGTAGCGCCTTATTGAGTGTAACAGGTTCACTATCTCGACTTTTTCGCCATGTTCTATCCATTTTGCTTTCATCAGCTTGGATGGAGATTCTGCGTTTCTCTATTGGACCCGATTCAAATATAAAGATGGACGGCGGTGATTCCGCGAAGACGTGCTTACACCAAAAGGATCGGCGTTCGCTCTGGGCGGGTAACGCCCCGGAGCTCCGTGGTGTGATAAACCCAGTGGCAAGGCACAGTCCCAGCTGGACAGCTTCTATTTCATAGTTAGCGTAGAGACCTGTACCCACATTACAAAGGTGACTCGGACCTTACGTCTTATTATGCCCTTTGACGCCAAATAATGGCGGATAGAGACCAATACTCCCTATATCGTAGTCTTATGTTAACATTCCCTACTTACTTCGCGTTCATCGAACATGCGCCCGTATCTGACGTAACGGTTCGAATACCTCCACAGATATGCGCCTCGCATCTTGTTCCTCCGTTACTGCTTGGAGGAATACGGAACCAGGTTATTTACTAGGTTTCAAAGATAGGTCGGCAGGCTCTCTTGGAAAATACTCCACCTAGACTGGTGCTGGCTGATACAACGATCATCGCGCCATCTCACCCGACCTTAGCCGCACAGTATTCGGTTAGCAGGGCCCCATTCTGATACGGAAACTTGCCAAGAGATGCTACCCTAACTATGGAGCAGGCCAGTTTTGGGGTCACACAGTTGCCTCTGGCGCGTTTACACCGCGTGCTACGGATAGGATGATAATCCTAAATACTTAGAAGGCATGTCTCGGGCCTTGTACGTCGTGCCGATACGCTCATAGAACCACCAGGCTAATATGCTGCTAGTGGCTAGGTCTGAGCGATTACGTAAATAAGTGTGCACCCCTCCGCGATACCTGAGTCTCATACCGCGCTCTGGAGTCGACGGTTTCAGGACGGTCTGAGGCGTGTTTCAACTCAATACTCTAAGACATCCTGAATCTTACTGCGCTTGTGGGCGTCCGAGTGCATTCTCGCTGCGTGAGAGAGCGGTCAATTGCTTGGCGGCAGAAGGACCGCGGTGGTGGTCTAGCTACACCCTACTTAGGCGACACCAAAGGTTCTGATAAGAATTGCACCGGGGATAGAAATTCTTGCGTACCACTTAGAAACATTGAGAAAGGGCATCATCTGACATGGGTGAGTTAACTGGCTATCTCAATTTGGTTGCCATCATTCCGAACTAGTTACTCGACTTTGGCTCGCAAGGTATCTTGTACGTAGTACGCATGTTTGGGATAGATTATTTAATTCCCTGGTAGACGATCGTCGCCATCTGATCCTCTCCGTGGCATCAAATTAGTTACAAGCCGGCGTATGCCCGTACATTTCTCGTTAACCAGCGGCGGGGATATGAAAGAAGGATATCCCCGTGTAGTCGCAAATCGTAAAGTTAGGGTATCCTCATCGACCAAATATTAACGGCCATTGTGTGTGGTAATAGCCAGAACAATGGGGCATGGTACCAATGTTAGATAACTGTAGAATAGTTCTCGTCGATGAGTGGGGGGCCAGATGTTCACCAAAGTGACCCATCAAAGAATTAACCCCCTCCGAAGGATAAACCGTGAAGGTGAGTGGCGAAAAGGTCGGATGCACTGCGGGACGACAATCATTCATATAGTTGTAATCTGCTAGCATTCATTCATGAGTCAAAAGAACGATCTTTTCGCTCGATCGGGATTACGTACATGTCCATAGCCCGGCGATCACTCTCACAGAACTTTCGATAGAGCCGCATCCGCTAGTATCCGTTTACCAGTATAATCAGGAAGTTACTCAGGTTCGAATTTACTCGGAGAGATTGGCCTGGTCGCCCCTTTTCCCGAATCCAGCCTGATCCACTCGTATGTTGCAGCCGCTCACCACGAAACCGTTCGTCACTGCCTAGCCCTCCGGCAGCTATTCCATAGTAACCTCTCATGGCACCCGAATATCGAACGTAGCTGCAATGAGGGTGGAGAGATTTTTCAACAAAATGCTACAGCTAATGTCTCATACTTTACCGATCCCTGATCATGATAACATTCGCATGACAATCAGGCGGGCTATGGAGGCATCCCCCGACTTTACGGGTACCTAACCGGAATGGGCGATTTCTTGGATTTTAGTCTCGTTGAGCTTCCGTTCGAACCTTGCAACATCGCATATATCAAGTACGACAGTTGGATCTTGGTATCACCCAAACTTCTTATTAAAGAGAGGGAGACAGTCTTACATCTTTTGAACCCTGGTCGCGACCCCCTACTCGTCATTATAACAGTCAGAATGCGTTTTTGTTCAACACTCTGAGGTGATGGAAGAACGCAACTTGACTCTGTACGTACCTACACTGTGCATAGTTCTTGAAGTTCGAAAGACAGGTTGACTCATACAGTGTTAGATGTGAGGAGCACGCCGCATCTTCTCCTGAGCGTTTACTTGATCTCACGGGTCTTGCGCGGCCGCATCGTCCTATAGCGATCCGTTGCAAAGAGTCCTTCTTTGATTCCAGCACCGAAGCCAGTCCAAGCATGGAGACCTGAGGCAGCCTCGTTAGAGTTGTCAACGAACGGTAAGTGCTGTAGAGAGGAACATTCCTAGCAGTTCACCGACGATTTATGGGACTTGCCAGTTCAAGCAGCGAGGTGTATAGACTTAGCGCAAGCTTACGCGTGACTCCTTTGACGCGCTCGGAGAGATGCATATATGGCAGGTCTCCAACGTACGGGTGTAAGATTCCTTTGTTTTCCGACTTTGATCCGATGGACCGCCGCATGTGCTCTAGGGTCGATAGCCTTAGAATATTCCGGAGACACCACACAGGAGGAATTAGGCCTTTCACGTGCCTATAACACCTATACCTGTCTTATACGGGTAGGACCCCACGCTTAGTCGCTTGCATGACTCGCTGCAAGCTTCCGTTTATTGGCTATAGCCTTTGCACTGCGGAGTGGTGGGAACGTCACCATCGGAACGAGCCTCACGTGAATCGTTGGGCCCCTCTCGTTTGCACAGTGGATAGGGGCTCCTGGAGGATAGACTTACCACGCGGAGCGAAGCCAGCACAAGTTTCCGGATCGGGGTGGATACGGCAAAAACCAAGGTATCGCACTGTTGGCGTCCACAGAGGAGTATAGGCTAATGAGGATGCTCTCGACTGCCCCTACTTAGCGTCATCGCCCGTGCCTGTGTACTCCTTATCTGTGGTTCGATTGCACTCCTGGGTTGGATGTCTTGATGGAAATGAGAGTGAAACTGCGATCCACTCACGGTGTAGAGGAAGATCTGCTCTGCTGAATTCTTTCTTCTGAGTAAGGACTTCGATCCGCCGGAGCTTCACCAATCGCGAGTGCCAAGTTGCCTATTGCCACAGGCCGAGTCGTACGGTACCCTGGCGCTAGGGAGCCAATCATAACTCTATAACGAAGAGAGTGTACACTAAGGACGGGCTCTCCAATAACCGAGAGGGCAGTATAAAGATTATCATATGCCGTTATTCTTTGAAAGTCTCATGCTCGAGAATATACATGATCGGTCCGGGCATAAGCCTGCTCAACAACCAGATTTTGTGTTATTGTCTGGCGCCCGACTAATTGCACACCTGCACATGCTTTCACCAGGCACGATTCACGCAAGGATATGCGGACAGCCTAGGTGCTCCGGAACCCGGCTGACCCGTTATCGGTAGCGGTTGGTACCATCACTGCTTACATAGGGCTGGCCATCTGGGTTCTCCCCCTTGTTGCTATGTTAACCACTACCGAAAATACCCGTGTTCCCTAAGTGCTATTGAGACACTATCTTGGATAGAGTATAACTGGTCATGAGCACCCGATGTTAATTTCGTATGTTACGAAGGCACGAGCTGTATTCGAAGTAGTCCTCTCGCGATGTGAAAAACGTCCAGGAGGAAGTAATGATAGTCAGGACGTAGCAGAGAGCCAAATCTGGTTTGCTAGACCGCACACTCAATGGACCACTGGTTGAGTCCAATGTTGGAACAATGTTAGTCCGCCAAGTTAATCAAAGCGTTGCAACACTGTCAAATCTTAAAATACGTAACACGGGTATTAGTAATTCATAATACGACGAAACGCGCAAATATCGGCGTGTTTTATTCGTCTGTGGGCGGCCGTCTTTTACCTTTACGACATATATTCCAGAACTAAGGGTAGTCCGGTCCAACTACTTGGAAGCGGTTGGATCGCCGTTAATCTACGTTGCACGAGGCGCGCAGCCCTCCCCCTAAATGGTCGGGCTACCTCGCTAGAATTTAACGTCCTAGGTCGTCTCCTCGTGGTCACTTAGTTATGGCGGCTTCTCTACGCACGTGCTCAAACTATTTAGAGGGTGCCGGTCGCCCGTATTTAACTCTTAGGGCCTTGAAGTACATATCGAGTGGCACTCAGGCTCAAAATTCCCAGTCCTTGGCATAAAGGGACTTTGCGTCAATTTTCAAGTCAATACCTGTAGTTCAATTAGGCGGGTAACACTCTGCTAAAGTATATATACTCGAAATTACAGGCGCTACCAACTAGTTGGGGAAAGTCCACATGTGGGATAATCGATCTTAGACGTGAGCCATCCCCTGGCGTGCTAAGAGGCCACTCTTTGCGGCCTCCACGGGCAAGATCGTGTATAACACGAAATGCCAGCGGACGGTCGGCGAATCCCTTGTCATATATGTATAGAATCCCTACACCTGCTGCTGGTTATCTGGGGTCTCGAATGGGGCGCGGGACGCCGTCGCAAAGGAGATAATGTCTAGATTCACGGAGTGTTGTAAGATGGTCGTGGTGGCACATCGAGAGAAATGTCAGCTTACCTGTCTCTTGGCTTGGATCACGGAGGGCACCGGCATGAGCTAGTGCATTACTGGCGATCGCACAAAAAGCCCGTCGAGTCCAACGCGGATGTGTTTTTAGGGGGACGTTGCCGGCAAAGTCCGGATGACATAAGGTGCGGGTCTTCCCGGCGCAATGATTCAAGCATGTGTTCGACCCTAAAATCTCGCAGGTGGTTCAAGTCATCAACACACTCAGTATTATCTTTTAGCCCGTGAATGGATTGACTGGGTCCCTCATAACCCGGCAACGTTCAGATGGACCGATCTACTGTTTTTCCCTAACGAGTTTTAAGACTGGGGTGTGCCCGGGCTTAATGCCTGGATTTACACCGATGGGCTAGAGATGCGAGTCTATAAAATGGCACCCAGCTTCTCGAGGATATCGGGCAGTGTCCTATTTCCCTTTTCGTGAAGGACGGACACGAATAGCCGTCAGATGCGATAACTGCACTTTAATGCATCACGCTCTAGACCCACGTTCTGGTGCCACACGGCTGAGCATTATGCACTCAAAGTCCAAAACCGTGGAGGATGCGGCGCAGATACCACGCCCTTGTGACGGGTTCCCGCGTTGCTATAGCTCCTGGTAGTGCGGAATATACGAACACCCACAGAATGGCTGTGCCTATGGTCACGTGTCAAACATGGTATGAAGGCGAATACAGACCTACTAGAGACAGCCCACAATCAGCTGTAACTATGAGTAGGCGAGTTGGCATACTCCAACTACTTGCTCTGGTCGTGACTTTTCGTGTTGGCCGTTCCTAGCCAGCCGTGCCTCTATGAAGGCCTATCCCCTCACGGCCACGCGGTAATCCAATACATCTCACTCACCGAATGGGTTCGGGAGCCGCGAGGATGAAATTTCTAACTTACTCTTTGGGGGGTTGACACCACAACAAACTCACGAAATGGGGTCACACGTTACACTCGGATTTCCTGGGTCATCATTCGTAAATGGCGCATTACGGAAATTTGTTTCGAAAGGATGTTAACTATGAGCGCCATCTGTTCGACCTTCAACTGCAGCCCTTGCGACAGTTCAATAAGTGAATTGGACCGGTGTCTGAGGAACAGTTGGTGCATTGTTTTAACAACTGCCGTTCTTAGTAAGGCCAATGAAGGAATCGGATGTATTAAATACATCTCACGCCACGAGTAATCCAGACGAACGCTGTAAAGTATTACTATCGACAATTAGATTTCTTCGTGTTACAAGTGGGCCGCGAGTTTATCGGCCCTCCGCTGCCAGTGATCGTATTCCACGTAAGTTGTCCCCACCGTGAGACTTGGAGCCGACAAACGTGTCGATGGGGGAGCCCGGCTCTCAATGTAATGTAAAGCTATGACTACAGCGTTGTCCCTCTAAATGATGGGCGTCGAACCTTCTGTCCTGACGCATCAGGGGTGTAAAATCAATGCTGTGCGAACGTCGCATCGACGAGGAAAGCCGCCGGGCCTCCAACTAATCACATACTTGGTCAATGTGTTGTGTTAGGTCGATCTGTCCAGCCATAGCGGGATCAGTACTAGTAAACGGACTAGGGCGTCGTACATGGCGGTCAGGTTGTGACCGGGATAAGGGGCTTTCAGTGCAACTAGGCCGGAGTGAGTGGAGACGATTGCGTTACCCTCACACTACGGTAGTTACGATTGACACTTCTTAGCGTTTCGCGCGGCCCCATGCACGCAATTTGTGAGTCTACTGAGATTTTGACATACGCACGAGTCGCCAAAGCCAGCCTATCACTCGCGACAATTACATATATTGAAAAAAGAAAGCCGGCGCACGATATGTAGCGAGCATCACTTGGTGGACCGAGCGGGCGTCATTCACTTGCGAACCCCCTTGAGTTGACTAAAGTAGAGACCGAGACTTACACAATAGTACTTATGTACTCGAAAAGCGATACCGGGTACCTCCAGAAGGCAAGAACGGACATGGACGGAGTCCAGATCTGATGGTGGGTATTCAAGGTGCATGGGCGCTCAGCCCTTCAACGAGCGCTCCTCCGCTGGCCGGCTTGGTAGATGGATTCTGGAACCATTGTCGAGACTGCTTAATCGTACGCACGACTTGGTGAAAGCTTTATGGACGTCGTAATCCAATCAACTCTGTGGCGAACGTGTTGGCCTTAATGCTTGCACACGATCGACTACCTTGTAAATCAGAAGTGCGCCTTTAACAATAAGAGCAATAATGTCGTGCAGTCGCGAAGAGGCCTGTAGCGCGGTCTATCAATGAATTGGCGCGCTTCTTGCCGTCGGGGTTCATCAGCGAGCCCACGAAGCTGCCTAGCGTGTATGCTAGAAGTCATAACTCCCTCGCGGAGCGCCGTTTCCTTTCCTGGCTTTCACCTAACAAAAAATGTTATATGGGCACTAGAGGCATAATTGCCTCTGACCACAGTGCCACGTCAACTGCCCTTGAGATATGGTATTCTAGGGTCATCTGAACCACAGGTCCCAGGCCACATGACGTAAACGGTAAATAAGGCTGCGGATTCCTGTCTGACGATCGCCCCTGTTCTAGACCAAAGCCATTGTGCGTTAACGTGGATTTTATGGTCCCCAGAATTGGACCTCAACCCATCGGATAAGGTAATTGCGCCTTTTTCCTCTACGGGTCGACCACCGGCTCCTTGAATAGATCACCAGTTTGGTAGATAACTATACGCGTCGCCCCTAGCTGTTCGGGACTACCTTTCGCTTCACCGTACCATGATGGTAGCGCATCCCCGGTTCCCAAAGATTGTGCTGCACTCTCCCCTGGGTGGTCTCGGCCCTTCCTAGGAGAGTTTTCGAAAAAACCTTTTGATCACTTCCCAA"

	print('Reverse complement for gemone of "{}":'.format(genome))
	print('\n\n')
	print(reverseComplement(genome))
