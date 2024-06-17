"""
Runnable File containing a minHash class:

TODO wip
python minHashFilter.py ____
"""

from test_filters import Filter, testFilter, readInputs
from collections import defaultdict 
import numpy as np
import sys

import re

class minHashFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, param):
        self.numHashes, self.numKmers = param
        self.numKmers
        self.title = f"minHashFilter_numHashes_{self.numHashes}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix
        Using Counts is very strict
        Using presence is very lenient for our sequences that are already very similar
        Maybe we should bin counts? i.e. take count // (len(seq)*p) for some p in (0,1)
    2. Apply MinHashes to Characteristic Matrix
    3. Bin 
    '''
    def preprocessReads(self):
        self.buildCharacteristicMatrix()
        self.minHash()

    """
    Matrix of seqs x kmers containing the count of each kmer by each sequence.
        Effciently this would just mark the presence as 0 or 1,
        We also can just use a subset of kmers which would be more efficient.
    """
    def buildCharacteristicMatrix(self, k = 6):
        # count all frequencies of kmers.
        m = np.power(4,6)
        self.characteristicMatrix = np.zeros((m,self.n), np.uint32)
        kmers = self.allKmers(k, [""])
        for i, seq in enumerate(self.seqs):
            for j, kmer in enumerate(kmers):
                self.characteristicMatrix[i,j] = len(re.findall(f'(?={kmer})', seq))
                """
                If we wanted to use original formulation of presence or no presence:
                This is a much more lenient definition of identity.
                self.characteristicMatrix[i,j] = len(seq.find(kmer))
                """

    """
    
    """
    def minHash(self):
        
        pass

    """
    Appends all possible kmers to each sequence of list and returns.
    i.e. [""] will return the list of all possible kmers.
    """
    def allKmers(self, k:int, jmers):
        if k == 0:
            return jmers
        kmers = []
        alphabet = "ACTG"
        for mer in jmers:
            for char in alphabet:
                kmers.append(mer + char)
        return self.allKmers(k-1, kmers)

    '''
    Connect elements in the same bucket.
    '''
    def connect(self, i, j):
        pass

def main():
    saveFig, param = readInputs()
    filter:minHashFilter = minHashFilter(param)
    testFilter(filter, saveFig = saveFig)

if __name__ == "__main__":
    main()