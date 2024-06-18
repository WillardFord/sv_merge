"""
Runnable File containing a sketch filter class:
    This approximates the cosine distance between two nodes.

python sketchFilter.py --plot --test --param 1000,3

>> Still experimenting with params in this filter: numHashes,K
TODO: Banding parameters are basically random. We need a better way to choose them.
"""

from test_filters import Filter, runFilter, readInputs
from collections import defaultdict 
import numpy as np

import re

class euclideanFilter(Filter):
    '''
    numHashes:  Number of Hashes to use in signature matrix >= 1000?? 
        Very large, we rely on law of large numbers for correct approximation.
    K:          Length of kmers for characteristic matrix, >= 2
    threshold:  Desired threshold for grouping together sequences. between 0 and 1
        b:
        r:
    binWidth: The desired binWidth in 
    '''
    def __init__(self, param):
        params = param.split(",")
        self.numHashes, self.K= int(params[0]), int(params[1])
        self.title = f"sketchFilter_numHashes={self.numHashes}_K={self.K}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix Using Kmer Counts
    2. Apply projection approximation for Euclidean distance to Characteristic Matrix
    3. Amplify the metric using band method
        TODO: Formal analysis of the effect of AND/OR gates used here.
    '''
    def preprocessReads(self):
        self.m = np.power(4, self.K)
        self.sketchSignature()
        self.band()

    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is --
    """
    def sketchSignature(self):
        self.signatureMatrix = np.zeros((self.numHashes, self.n))
        hashes = self.getRandomPlanes()
        for j in range(self.numHashes): # iterate planes
            for i in range(self.n): # iterate strings
                self.signatureMatrix[j,i] = self.sketch(self.getCharacteristicVector(i), hashes[j])

    """
    Return the characteristic vector for a single read, i
    """
    def getCharacteristicVector(self, i):
        characteristicVector = np.zeros(self.m)
        for j, kmer in enumerate(self.generateKmers(self.K)):
            characteristicVector[j] = len(re.findall(f'(?={kmer})', self.seqs[i]))
        return characteristicVector

    """
    Generate all possible kmers without storing them in memory.
    """
    def generateKmers(self, k:int):
        alphabet = "ACTG"
        if k == 1:
            for char in alphabet:
                yield char
        else:
            for mer in self.generateKmers(k-1):
                for char in alphabet:
                    yield mer + char

    """
    This function returns a random set of binary vectors in [-1, 1]. 
    They function as normal vectors to separating planes used as hashes
    """
    def getRandomPlanes(self):
        dimensions = self.m
        lines = np.zeros((self.numHashes, dimensions))
        for i in range(self.numHashes):
            lines[i,:] = np.random.choice([1, -1], size=dimensions)
        return lines

    """
    Returns 1 if dot product is positive, else -1. (0's are randomly assigned)
    This approximates the cosine distance though is not exact.
    """
    def sketch(self, characterisicVector, plane):
        value = np.dot(characterisicVector, plane)
        if value > 0:
            return 1
        if value < 0:
            return -1
        return np.random.choice([1, -1])

    """
    TODO: choose bands to be more optimal for filter and context
    Use banding technique to bin reads by similarity
        -- Choose threshold, t between 0 and 1 to determine our false positive rate. 
            Higher threshold indicates less false postives
        -- b*r = n
            w/ b = number of bands, r = band length, n = number of hashes
        -- t ~= (1/b)^(1/r)
    """
    def band(self, threshold = 0.85):
        def connectBucket(bucket):
            for i in range(len(bucket)-1):
                for j in range(i+1, len(bucket)):
                    self.adjacencyMatrix[bucket[i],bucket[j]] = True
                    self.adjacencyMatrix[bucket[j],bucket[i]] = True

        numBands, bandLength = self.getNumBands(threshold)
        for i in range(numBands):
            buckets = defaultdict(list[int])
            for j in range(self.n):
                buckets[self.signatureMatrix[(i*bandLength):((i+1)*bandLength), j].tobytes()].append(j)
            for bucket in buckets.values():
                connectBucket(bucket)

    """
    Determine the correct number of bands and their length given a desired threshold maybe??
    TODO: Need much better method for determine the number of bands.

    numBands and bandLength define each other given numHashes. So just need to decide which to optimize.
    """
    def getNumBands(self, threshold):
        numBands = min(40, self.numHashes) 
        bandLength = self.numHashes // numBands
        return numBands, bandLength

    '''
    Connect elements in any of the same bucket.
    '''
    def connect(self, i, j):
        self.adjacencyMatrix[i,j]

def main():
    saveFig, param, test = readInputs()
    filter:euclideanFilter = euclideanFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()