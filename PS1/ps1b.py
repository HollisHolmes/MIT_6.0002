###########################
# 6.0002 Problem Set 1b: Space Change
# Name:
# Collaborators:
# Time:
# Author: charz, cdenise

#================================
# Part B: Golden Eggs
# Part B: Golden Eggs
#================================

Problem 1
def dp_make_weight(egg_weights, target_weight, memo = {}):
    """
    Find number of eggs to bring back, using the smallest number of eggs. Assumes there is
    an infinite supply of eggs of each weight, and there is always a egg of value 1.

    Parameters:
    egg_weights - tuple of integers, available egg weights sorted from smallest to largest value (1 = d1 < d2 < ... < dk)
    target_weight - int, amount of weight we want to find eggs to fit
    memo - dictionary, OPTIONAL parameter for memoization (you may not need to use this parameter depending on your implementation)

    Returns: int, smallest number of eggs needed to make target weight
    !!!!
    State is the least number of moves given a remaining weight. No state can depend on another with a larger weight.
    Hence, the state map can be topologically sorted and we can run DP!!
    !!!!!
    """
    # memoize min eggs for each weight
    memo[0] = 0
    for weight in range(1, target_weight+1):
        # base case is it takes 0 moves to get to 0 from 0
        turns_at_this_weight = []
        # bottom up approach starting at weight 1
        for egg_weight in egg_weights:
            # weight remaing respective after subtracting each possible egg weight
            after_adding = weight - egg_weight
            # if we can get to this state
            if after_adding in memo:
                # we can get to 0 in 1 + how many it takes in the new state by optimal substructure
                turns_at_this_weight.append(1 + memo[after_adding])
        # we have # turns for each egg weight, only store the best option
        memo[weight] = min(turns_at_this_weight)
    # return min at the targen weight in O(n) time
    return memo[target_weight]

# (1, 5, 25)







# EXAMPLE TESTING CODE, feel free to  add more if you'd like
if __name__ == '__main__':
    egg_weights = (1, 5, 15, 20)
    n = 99
    print("Egg weights = (1, 5, 10, 25)")
    print("n = 99")
    print("Expected ouput: 9 (3 * 25 + 2 * 10 + 4 * 1 = 99)")
    print("Actual output:", dp_make_weight(egg_weights, n))
    print()

#### Greedy algorithm implementation that could return incorrect values   ####

# def dp_make_weight(egg_weights, target_weight, memo = {}):
#     """
#     Find number of eggs to bring back, using the smallest number of eggs. Assumes there is
#     an infinite supply of eggs of each weight, and there is always a egg of value 1.
#
#     Parameters:
#     egg_weights - tuple of integers, available egg weights sorted from smallest to largest value (1 = d1 < d2 < ... < dk)
#     target_weight - int, amount of weight we want to find eggs to fit
#     memo - dictionary, OPTIONAL parameter for memoization (you may not need to use this parameter depending on your implementation)
#
#     Returns: int, smallest number of eggs needed to make target weight
#     """
#
#     weightRemaining = target_weight
#     numEggs = 0
#
#     for i in range(-1, -len(egg_weights)-1, -1):
#         while weightRemaining >= egg_weights[i]:
#             weightRemaining -= egg_weights[i]
#             numEggs += 1
#     return numEggs
#
