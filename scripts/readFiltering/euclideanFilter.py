"""
Runnable File containing a euclidean class:

python euclideanFilter.py --plot --test --param 1000,20,5,40

>> Still experimenting with params in this filter: numHashes,K,binWidth,bandLength
TODO: Banding parameters are basically random. We need a better way to choose them.

This filter is sensitive to read orientation in ways the other filters are not.
This makes sense theoretically, and explains the variance we were seeing where half of our reads didn't seem to group nicely.
"""

from test_filters import Filter, runFilter, readInputs
from collections import defaultdict 
import numpy as np

import re
import os
import subprocess

class euclideanFilter(Filter):
    '''
    numHashes:  Number of Hashes to use in signature matrix >= 1000?? 
        Very large, we rely on law of large numbers for correct approximation.
    K:          Length of kmers for characteristic matrix, >= 2
    threshold:  Desired threshold for grouping together sequences. between 0 and 1
        b:
        r:
    binWidth: The desired binWidth in a single random line for aligning.
    '''
    def __init__(self, param):
        params = param.split(",")
        self.numHashes, self.K, self.binWidth = int(params[0]), int(params[1]), float(params[2])
        self.bandLength = int(params[3])
        self.title = f"euclideanFilter_numHashes={self.numHashes}_K={self.K}_binWidth={self.binWidth}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix Using Kmer Counts
    2. Apply projection approximation for Euclidean distance to Characteristic Matrix
    3. Amplify the metric using band method
        TODO: Formal analysis of the effect of AND/OR gates used here.
    '''
    def preprocessReads(self):
        self.fracMinHashKmers()
        self.projectionSignature()
        self.band()

    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is the projection of the seqs characteristicMatrix
        represenation onto a random line.
    """
    def projectionSignature(self):
        self.signatureMatrix = np.full((self.numHashes, self.n), np.inf)
        for i in range(self.n): # iterate strings
            for j, (direction, offset) in enumerate(self.getRandomLines()):
                charVector = self.getCharacteristicVector(i)
                self.signatureMatrix[j,i] = self.projectAndBin(charVector, direction, offset)

    """
    This function yields a random set of unit lines to use with projection as hash functions
    Uniformly distributed lines on unit hypersphere: 
        https://math.stackexchange.com/questions/444700/uniform-distribution-on-the-surface-of-unit-sphere
    """
    def getRandomLines(self):
        seed = 1
        np.random.seed(seed)
        dimensions = self.m
        for _ in range(self.numHashes):
            s = np.random.normal(0, 1, dimensions)
            direction = s / np.sqrt(sum(s*s))
            offset = np.random.uniform(0,5,dimensions)
            yield direction, offset

    """
    Compute projection of point onto a random line and calculate which bin the point falls into.
    Returns pos or neg integer indicating which bin projected into. i.e.
        |-2|-1|0(contains origin at leftmost boundary)|1|2|
    """
    def projectAndBin(self, characteristicVector, line, offset):
        projection = np.dot(line,characteristicVector)*line + offset
        projectedWidth = line[0] * self.binWidth
        return int(projection[0] // projectedWidth)

    """
    Connect elements in any of the same bucket.
    """
    def connect(self, i, j):
        return self.adjacencyMatrix[i,j]

def main():
    saveFig, param, test = readInputs()
    filter:euclideanFilter = euclideanFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()