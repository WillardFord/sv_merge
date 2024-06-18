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
        self.buildCharacteristicMatrix()
        self.minHashSignature()
        self.band()

    """
    Generate Matrix of kmers x seqs 
        -- In practice this matrix is too large to store in memory. But it's unclear how to operate instead.
        -- Binary values indicating presence or not

        TODO: Test whether or not specific subsets of kmers efficiently provide sufficient information for identity
    """
    def buildCharacteristicMatrix(self):
        kmers  = self.getKmers()
        self.characteristicMatrix = np.zeros((self.m, self.n), np.bool_)
        for i, seq in enumerate(self.seqs):
            for j, kmer in enumerate(kmers):
                self.characteristicMatrix[j,i] = kmer in seq

    """
    Generate signatureMatrix of hashs x seqs
        In practice use hashmap not permutation, but simpler in Python
    """
    def minHashSignature(self):
        self.signatureMatrix = np.full((self.numHashes, self.n), np.inf)
        rng = np.random.default_rng()
        hashes = [rng.permutation(self.m) for _ in range(self.numHashes)]
        for j in range(self.m): # iterate kmers
            for i in range(self.n): # iterate strings
                if self.characteristicMatrix[j,i] > 0:
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
        def swap_keys_values(dictionary):
            new_dictionary = {}
            for key, values in dictionary.items():
                for value in values:
                    new_dictionary[value] = key
            return new_dictionary
        def get_numBands_bandLength(threshold):
            numBands = min(50, self.numHashes) 
            bandLength = self.numHashes // numBands
            print(numBands, bandLength)
            print(f"Threshold {np.power(1/numBands, 1/bandLength)}")
            return numBands, bandLength

        # TODO: Automate how b, r are chosen. 
        # These values give about the desired threshold 
        numBands, bandLength = get_numBands_bandLength(threshold)
        numBands = 40
        bandLength = 25
        self.bucketsSet = [defaultdict(list[int]) for _ in range(numBands)]
        for i in range(numBands):
            for j in range(self.n):
                self.bucketsSet[i][self.signatureMatrix[(i*bandLength):((i+1)*bandLength), j].tobytes()].append(j)
            self.bucketsSet[i] = swap_keys_values(self.bucketsSet[i])

    """
    TODO: I can make this a yield function as well. 
        Then I just have to define m to be large enough without throwing errors
        And then set M to it's minimum value as soon as we find out what it is.
    """
    def getKmers(self):
        kmers = self.fracMinHashKmers()
        return kmers

    """
    Using a subset of the total kmer set allows us to compute much faster. We get this set by using
        Using s = 1000 as suggest default. This allows large k's to be used. (21, 31, 51)
    
        FRACs(W ) = { h(w) ≤ H/s ∣ ∀w ∈ W } 
        from https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2.full.pdf
    """
    def fracMinHashKmers(self, s = 0.5e-2):
        kmers = []
        maxVal = 0xFFFFFFFFFFFFFFFF # maximum hashable value on base 64 system.
        limit = maxVal*s
        print(limit)
        for kmer in self.generateKmers(self.K):
            #print(abs(hash(kmer)))
            if abs(hash(kmer)) < limit:
                kmers.append(kmer)
        self.m = len(kmers)
        print("Downsampled to:", self.m)
        return kmers


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

    '''
    Connect elements in any of the same bucket.
    '''
    def connect(self, i, j):
        for buckets in self.bucketsSet:
            if buckets[i] == buckets[j]:
                return True
        return False

def main():
    saveFig, param, test = readInputs()
    filter:minHashFilter = minHashFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()