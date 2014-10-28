# coding=utf-8
_NUCLEOTIDES = ['A', 'C', 'G', 'T']




def patternToNumber(pattern):
    """
    Transforms a k-mer Pattern into an integer to mark its index location.

    Our set of locations will range from 0 to ((4^k)-1)
    """
    if not pattern:
        return 0

    symbol = pattern[-1]
    pattern = pattern[:-1]
    return (4 * patternToNumber(pattern)) + symbolToNumber(symbol)



def symbolToNumber(symbol):
    """
    Converts an A, C, G, or T into a number according to its lexicographic order
    """
    if symbol in _NUCLEOTIDES:
        return _NUCLEOTIDES.index(symbol)



def numberToSymbol(number):
    return _NUCLEOTIDES[int(number % 4)]


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

    prefix_pattern = numberToPattern(prefixIndex, k - 1)
    symbol = numberToSymbol(r)

    return prefix_pattern + symbol



def hammingDistance(p, q):
    """
    We say that a k-mer Pattern appears as a
    substring of Text with at most d mismatches if
    there is some k-mer substring Pattern' of Text
    having d or fewer mismatches with Pattern.
    """

    if len(p) == len(q):
        return sum([1 if c1 != c2 else 0 for c1, c2 in zip(p, q)])

    if len(p) > len(q):
        return sum([1 if c1 != c2 else 0 for c1, c2 in zip(p[:len(q)], q)])

    if len(p) < len(q):
        return sum([1 if c1 != c2 else 0 for c1, c2 in zip(p, q[:len(p)])])

    else:
        return 0



def quotient(n, m):
    return int(n / m)


def remainder(n, m):
    return n % m


def make_spaced_string_from_comma_separated_array(array):
    return " ".join([str(x) for x in array])
