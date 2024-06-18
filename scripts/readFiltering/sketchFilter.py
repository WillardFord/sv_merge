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
        self.buildCharacteristicMatrixCounts()
        self.sketchSignature()
        self.band()

    """
    Generate Matrix of kmers x seqs
        -- In practice this matrix is too large to store in memory. But it's unclear how to operate instead.
        -- count of each kmer by each sequence.
        -- this a very strict definition
    """
    def buildCharacteristicMatrixCounts(self):
        kmers  = self.getKmers()
        self.characteristicMatrix = np.zeros((self.m, self.n), np.uint32)
        for i, seq in enumerate(self.seqs):
            for j, kmer in enumerate(kmers):
                self.characteristicMatrix[j,i] = len(re.findall(f'(?={kmer})', seq))


    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is --
    """
    def sketchSignature(self):
        self.signatureMatrix = np.zeros((self.numHashes, self.n))
        hashes = self.getRandomPlanes()
        for j in range(self.numHashes): # iterate planes
            for i in range(self.n): # iterate strings
                self.signatureMatrix[j,i] = self.sketch(self.characteristicMatrix[:,i], hashes[j])

    """
    This function returns a random set of binary sketches. 
    They function as normal vectors to separating planes used as hashes
    """
    def getRandomPlanes(self):
        dimensions = self.m
        lines = np.zeros((self.numHashes, dimensions))
        for i in range(self.numHashes):
            lines[i,:] = np.random.choice([1, -1], size=dimensions)
        return lines

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
        def swap_keys_values(dictionary):
            new_dictionary = {}
            for key, values in dictionary.items():
                for value in values:
                    new_dictionary[value] = key
            return new_dictionary
        def get_numBands_bandLength(threshold):
            numBands = min(50, self.numHashes) 
            bandLength = self.numHashes // numBands
            #print(numBands, bandLength)
            #print(f"Threshold {np.power(1/numBands, 1/bandLength)}")
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
    Returns set of all kmers to use for generating the characteristicMatrix
        For large K we should use a random subset of all kmers but for small k we can use them all.
    """
    def getKmers(self):
        kmers = self.allKmers()
        self.m = len(kmers)
        return kmers

    """
    Return set of all kmers
    """
    def allKmers(self):
        kmers = []
        for kmer in self.generateKmers(self.K):
            kmers.append(kmer)
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
            if abs(hash(kmer)) < limit:
                kmers.append(kmer)
        print("Downsampled to:", len(kmers))
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
    filter:euclideanFilter = euclideanFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()