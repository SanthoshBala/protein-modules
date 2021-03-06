#! /usr/bin/python

# statistics.py
# Author: Santhosh Balasubramanian
# Created: March 24, 2013
# Last Modified: March 24, 2013


# Python Imports
import math


# - - - - - - - - - - SAMPLE STATISTICS - - - - - - - - - - #


# mean: Computes mean of <numberList>, along with standard error of the mean.
# Returns <avg> and <sem> as a tuple.
def mean(numberList):    
    numValues = len(numberList)
    numerator = float(sum(numberList))
    denominator = float(numValues)

    average = numerator/denominator
    error = average/math.sqrt(numValues)

    return average, error

# median: Computes median of <numberList>.
def median(numberList):
    numValues = len(numberList)
    midpoint = numValues / 2

    numberList.sort()

    # Return Mean if <numValues> is Even
    if numValues % 2 == 0:
        return mean( [ numberList[midpoint], numberList[midpoint - 1] ] )[0]
    else:
        return numberList[midpoint]

# variance: Computes variance of <numberList>.
def variance(numberList):
    average = mean(numberList)[0]
    diffs = [(s - average) for s in numberList]
    squares = [s*s for s in diffs]
    
    numerator = float(sum(squares))
    denominator = float(len(numberList))
    
    return numerator/denominator

# stdDev: Computes standard deviation of <numberList>.
def stdDev(numberList):
    variation = variance(numberList)
    return math.sqrt(variation)

# poll: Return % of values in <numberList> that exceed <threshold>.
def poll(numberList, threshold):

    # Count # Upvotes
    upvotes = 0
    for number in numberList:
        if (number > threshold):
            upvotes = upvotes + 1

    # Compute % Upvotes
    numerator = float(upvotes)
    denominator = float(len(numberList))

    return numerator/denominator


# - - - - - - - - - - COMBINATORICS - - - - - - - - - - #

# getCombinatorialPairs: Returns n choose 2 pairs from <itemList>.
def getCombinatorialPairs(itemList):
    pairList = []
    numItems = len(itemList)

    # Iterate through All Pairs
    for i in range(numItems):
        for j in range(i, numItems):
            if i == j:
                continue
            pairList.append( ( itemList[i], itemList[j] ) )

    return pairList


# - - - - - - - - - - SET MANIPULATION - - - - - - - - - - #


# getJaccardSimilarity: Returns Jaccard index of <setA> and <setB>.
def getJaccardSimilarity(setA, setB):

    intersection = setA.intersection(setB)
    union = setA.union(setB)
    result = float(len(intersection))/float(len(union))

    return result


# - - - - - - - - - - MISCELLANEOUS - - - - - - - - - - #


# hoursMinutesSeconds: Converts <seconds> into # of hours, minutes, seconds.
def hoursMinutesSeconds(seconds):

    h = seconds / 3600
    m = (seconds % 3600) / 60
    s = (seconds % 3600) % 60

    return h, m, s

# hashNumberToHistogramBucket: Hashes a <value> to histogram with <bucketSize>
# buckets, and returns the bucket number.
def hashNumberToHistogramBucket(value, bucketSize):
    value = int(math.floor(value/float(bucketSize)))

    return value
