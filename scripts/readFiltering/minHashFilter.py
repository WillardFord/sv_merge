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

TODO:
    1. Fix hashing to not rely on combinations and actually use hashing
    2. Implement binary representation of sequences so that inputs are uniformly distributed across the input space

python minHashFilter.py --plot --param 1000,13,40
                                    numHashes,K,bandLength
"""

from test_filters import Filter, runFilter, readInputs

import numpy as np

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
        self.bandLength = int(params[2])
        self.title = f"minHashFilter_numHashes={self.numHashes}_K={self.K}_r={self.bandLength}"
        self.hashes = 0

    '''
    1. Generate Signaure Matrix (Hashed data for all reads across all hashes)
    2. Bin using band method
    '''
    def processReads(self, test = False):
        self.minHashSignature(test)
        self.band()
        return self.signatureMatrix

    """
    Generate signatureMatrix of hashs x seqs
    """
    def minHashSignature(self, test):
        self.signatureMatrix = np.full((self.numHashes, self.n), np.inf)
        for k, univ_hash in enumerate(self.hashes):
            if k == self.numHashes: break
            for i in range(self.n):
                characteristicVector = self.getCharacteristicVector(i)
                for kmer in self.kmer_dict.keys():
                    if characteristicVector[self.kmer_dict[kmer]] == 0: continue
                    self.signatureMatrix[k,i] = min( self.signatureMatrix[k,i], univ_hash(self.kmer_dict[kmer]) )

def main():
    saveFig, param, test = readInputs()
    filter:minHashFilter = minHashFilter(param)
    runFilter(filter, saveFig, test = test, loadMinHash = True)

if __name__ == "__main__":
    main()