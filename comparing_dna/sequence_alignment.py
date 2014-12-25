import numpy as np

def longet_common_subsequence(string1, string2):
    """
    Finds A longest common subsequence between two strings using
    insertion (string1) and deletion (string2) methods
    """
    pass


def min_number_of_coins(money, coin_choices):
    """
    Demonstrates use of dynamic programming to find the minimum
    number of coins needed to attain some some, "money"
    """

    min_num_coins = [0]
    max_array_length = max(coin_choices)  # TODO: Ensure that the array never exceeds this size

    for i in range(1, money + 1):
        array_index = i % max_array_length
        min_num_coins.append(money * money)  # initialize a maximum for this array position at the start of the iteration

        for j in range(len(coin_choices)):

            # Only in this case can a choice potentially meet our needs
            if i >= coin_choices[j]:

                if min_num_coins[i - coin_choices[j]] + 1 < min_num_coins[i]:

                    # Update the lookup list
                    min_num_coins[i] = min_num_coins[i - coin_choices[j]] + 1


    return min_num_coins[-1]


def south_or_east(i, j):
    """
    computes the length of the longest path to node (i, j) in
    a rectangular grid in the Manhattan Problem,
    abiding by the observation that the only way to reach node (i, j)
    is either by moving south (↓) from (i − 1, j) or east (→) from (i, j − 1).

    :return: max(i, j)
    """
    if i == 0 and j == 0:
        return 0

    x = -1 * (j ** j)  # Initialize an impossible minimum for x and y
    y = -1 * (i ** i)

    if i > 0:
        x = south_or_east(i - 1, j) + southward_edge_weight(i-1, j)

    if j > 0:
        y = south_or_east(i, j - 1) + eastward_edge_weight(i, j-1)

    return max(x, y)