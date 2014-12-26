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



def longest_path_in_graph(n, m, down, right):
    """
    dynamic programming algorithm for finding the
    length of a longest path in the Manhattan Tourist Problem.

    Conceptually, we can think of down-i-j and right-i-j as being
    the respective weights of the vertical and horizontal edges entering node (i, j).

    We denote the matrices holding down-i-j and right-i-j as `Down` and `Right`, respectively.
    """

    # Initialize the graph
    graph = [[0] * m] * n


    if n == 0 and m == 0:
        return 0

    y = -1 * (n ** n)  # Initialize an impossible minimum for y and x (rows and columns)
    x = -1 * (m ** m)

    for i in range(1, n):
        graph[i][0] = graph[i-1][0] + down

    for j in range(1, m):
        graph[0][j] = graph[0][j-1] + right

    # if n > 0:
    #     y = longest_path_in_graph(n - 1, m, down, right) + down
    #
    # if m > 0:
    #     x = longest_path_in_graph(n, m-1, down, right) + right

    for ii in range(n):
        for jj in range(m):
            graph[ii][jj] = max(graph[ii-1][j] + down, graph[ii][jj-1] + right)

    # The final value should be a result of computing the longest possible path!
    return graph[n][m]




