"""
Runnable File containing a minHash class:

TODO 
    1) The numBands and bandLength parameters currently only work for a small range of peaks ~ 300k - 800k
            numBands = 40
            bandLength = 25
        We could fix this by using a tighter threshold as we get to larger k's, 
            but to reduce variance we're better off adjusting numBands and bandLength dynamically based on number of peaks.
    
    2) The numHashes parameter maybe should also be a constant. We should be optimizing for a desired threshold.
    
    3) We can also compute the kmer counts on the fly, skipping the characteristic matrix all together,
        but it's likely not neccessary for just the python implementation.



python minHashFilter.py --plot --param 1000,13
"""

from test_filters import Filter, runFilter, readInputs
from collections import defaultdict 
import numpy as np

import re

class minHashFilter(Filter):
    '''
    numHashes:  Number of Hashes to use in signature matrix >= 1000?? Very large, we rely on law of large numbers for correct approximation.
    K:          Length of kmers for characteristic matrix, >= 2
    threshold:  Desired threshold for grouping together sequences. between 0 and 1
        b:
        r:
    '''
    def __init__(self, param):
        params = param.split(",")
        self.numHashes, self.K = int(params[0]), int(params[1])
        self.title = f"minHashFilter_numHashes={self.numHashes}_K={self.K}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix
        Using Counts is very strict
        Using presence is very lenient for our sequences that are already very similar
        Maybe we should bin counts? i.e. take count // (len(seq)*p) for some p in (0,1)
    2. Apply MinHashes to Characteristic Matrix
    3. Bin using band method
    '''
    def preprocessReads(self):
        self.minHashSignature()
        self.band()

    """
    Generate signatureMatrix of hashs x seqs
        In practice use hashmap not permutation, but simpler in Python
    """
    def minHashSignature(self):
        self.signatureMatrix = np.full((self.numHashes, self.n), np.inf)
        rng = np.random.default_rng()
        kmers = self.fracMinHashKmers()
        hashes = [rng.permutation(self.m) for _ in range(self.numHashes)]
        for j, kmer in enumerate(kmers): # iterate kmers
            for i in range(self.n): # iterate strings
                if kmer in self.seqs[i]: # Characteristic Matrix Value
                    for k, perm in enumerate(hashes):
                        self.signatureMatrix[k,i] = min(self.signatureMatrix[k,i], perm[j])


    """
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

    """
    Using a subset of the total kmer set allows us to compute much faster. We get this set by using
        Using s = 1000 as suggest default. This allows large k's to be used. (21, 31, 51)
    
        FRACs(W ) = { h(w) ≤ H/s ∣ ∀w ∈ W } 
        from https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.full.pdf
    """
    def fracMinHashKmers(self, s = 0.5e-2):
        maxVal = 0xFFFFFFFFFFFFFFFF # maximum hashable value on base 64 system.
        limit = maxVal*s
        kmers = []
        for kmer in self.generateKmers(self.K):
            if abs(hash(kmer)) < limit:
                kmers.append(kmer)
        self.m = len(kmers)
        return kmers

    """
    Yield all possible kmers without storing them in memory.
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

    '''
    Connect elements in any of the same bucket.
    '''
    def connect(self, i, j):
        return self.adjacencyMatrix[i,j]

def main():
    saveFig, param, test = readInputs()
    filter:minHashFilter = minHashFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()