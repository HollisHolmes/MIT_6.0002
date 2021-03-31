###########################
# 6.0002 Problem Set 1a: Space Cows
# Name:
# Collaborators:
# Time:

from ps1_partition import get_partitions
import time

#================================
# Part A: Transporting Space Cows
#================================

# Problem 1
def load_cows(filename):
    """
    Read the contents of the given file.  Assumes the file contents contain
    data in the form of comma-separated cow name, weight pairs, and return a
    dictionary containing cow names as keys and corresponding weights as values.

    Parameters:
    filename - the name of the data file as a string

    Returns:
    a dictionary of cow name (string), weight (int) pairs
    """

    #open file
    file = open(filename)
    #split into different lines

    lineList = file.read().splitlines()
    #dictionary for output
    store = {}

    #split each line about the comma
    #enter first value as key, second value as number
    for l in lineList:
        splitLines = l.split(',')
        store[splitLines[0]] = splitLines[1]

    return store



# Problem 2
def greedy_cow_transport(cows,limit=10):
    """
    Uses a greedy heuristic to determine an allocation of cows that attempts to
    minimize the number of spaceship trips needed to transport all the cows. The
    returned allocation of cows may or may not be optimal.
    The greedy heuristic should follow the following method:

    1. As long as the current trip can fit another cow, add the largest cow that will fit
        to the trip
    2. Once the trip is full, begin a new trip to transport the remaining cows

    Does not mutate the given dictionary of cows.

    Parameters:
    cows - a dictionary of name (string), weight (int) pairs
    limit - weight limit of the spaceship (an int)

    Returns:
    A list of lists, with each inner list containing the names of cows
    transported on a particular trip and the overall list containing all the
    trips
    """
    #store list of cows on each trip
    tripList = []

    #copy dictionary into a list sorted in decending order
    cowList = []   #cows sorted by wieght
    for x,y in cows.items():
        cowList.append([x,y])
    cowList.sort(key = lambda x : x[1], reverse=True)

    #calculate the cows taken for each trip, while cowList non-empty
    while cowList:
        nextCow = []
        trip = []
        tripCapacity = limit
        #loop through each cow still on earth
        for item in cowList:
            #if there's space for them, add to the trip
            if int(item[1]) <= tripCapacity:
                trip.append(item[0])   #add cow name to the current trip
                tripCapacity -= int(item[1])
            #else leave them on earth for the next trip
            else:
                nextCow.append(item)
        tripList.append(trip)
        cowList = nextCow[:]
    #return the list of tripsc/u
    return tripList


# Problem 3
def brute_force_cow_transport(cows,limit=10):
    """
    Finds the allocation of cows that minimizes the number of spaceship trips
    via brute force.  The brute force algorithm should follow the following method:

    1. Enumerate all possible ways that the cows can be divided into separate trips
        Use the given get_partitions function in ps1_partition.py to help you!
    2. Select the allocation that minimizes the number of trips without making any trip
        that does not obey the weight limitation

    Does not mutate the given dictionary of cows.

    Parameters:
    cows - a dictionary of name (string), weight (int) pairs
    limit - weight limit of the spaceship (an int)

    Returns:
    A list of lists, with each inner list containing the names of cows
    transported on a particular trip and the overall list containing all the
    trips
    """
    minPartition = None
    minPartitionLen = 10
    validTrips = []
    #all possible splits
    for partition in get_partitions(set(cows.keys())):
        #print("partition:", partition, "number of trips:", len(partition))
        #for each trip in each partition, check they are all under weight limits
        #if all under weight limit retun the number of trip
        underCapacity = 0
        for trips in partition:
            totalWeight = 0
            for cow in trips:
                totalWeight += int(cows[cow])
            #print("the trip is", trips, "the weight is", totalWeight)

            if totalWeight <= 10: underCapacity += 1
        if underCapacity == len(partition):
            validTrips.append(partition)
            if len(partition) < minPartitionLen:
                minPartition = partition
                minPartitionLen = len(partition)
    return validTrips
# Problem 4
def compare_cow_transport_algorithms():
    """
    Using the data from ps1_cow_data.txt and the specified weight limit, run your
    greedy_cow_transport and brute_force_cow_transport functions here. Use the
    default weight limits of 10 for both greedy_cow_transport and
    brute_force_cow_transport.

    Print out the number of trips returned by each method, and how long each
    method takes to run in seconds.

    Returns:
    Does not return anything.
    """
    #load cows from text file
    cows = load_cows("ps1_cow_data.txt")

    #time greedy algo
    start = float(time.time())
    greedyList = greedy_cow_transport(cows)
    end = float(time.time())
    greedyTime = end - start
    print("The greedy algorithm takes", greedyTime, "seconds to run")

    #time brute force algo
    start = time.time()
    bruteList = brute_force_cow_transport(cows)
    end = time.time()
    bruteTime = end - start
    print("The brute force algorithm takes", bruteTime, "seconds to run")

#ps1_cow_data.txt

compare_cow_transport_algorithms()

#Greedy algorithm runs faster as it is limited to the set of possible outcomes,
#and even then, does not evaluate each possiblility for this set.

#In some cases the greedy aglorithm will not return the best result.

#The brute force algorithm enumerates every possible compibation of trips for
#the cows (partitions), including unviable trips where the weight is greater
#than the maximum allowed. The number of partitions a set has is given by the
#catilan numbers, which for a set with cardinality n is: 2nCn / (n+1).

#The brute force algorithm will always return the best result as it evaluates
#every possible solution.
