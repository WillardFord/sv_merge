"""
Runnable File containing a sketch filter class:
    This approximates the cosine distance between two nodes.

python sketchFilter.py --plot --test --param 1000,3,40

>> Still experimenting with params in this filter: numHashes,K,bandLength
TODO: Banding parameters are basically random. We need a better way to choose them.
"""

from test_filters import Filter, runFilter, readInputs
import numpy as np
import math

class sketchFilter(Filter):
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
        self.numHashes, self.K, self.bandLength = int(params[0]), int(params[1]), int(params[2])
        self.title = f"sketchFilter_numHashes={self.numHashes}_K={self.K}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix Using Kmer Counts
    2. Apply projection approximation for Euclidean distance to Characteristic Matrix
    3. Amplify the metric using band method
        TODO: Formal analysis of the effect of AND/OR gates used here.
    '''
    def processReads(self):
        self.sketchSignature()
        self.band()
        return self.signatureMatrix

    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is a plane.
        Here we load in a random vector but we'd likely generate one in practice.
    """
    def sketchSignature(self):
        self.signatureMatrix = np.zeros((self.numHashes, self.n))
        for i in range(self.n): # iterate reads
            characteristicVector = self.getCharacteristicVector(i)
            for j, line in enumerate(self.getRandomPlanes(characteristicVector.shape[0])):
                if j == self.numHashes:
                    break
                self.signatureMatrix[j,i] = math.copysign(1, np.dot(characteristicVector, line))

    """
    This function returns a random set of binary vectors in [-1, 1] that describe a plane.
    """
    def getRandomPlanes(self, length):
        planeFile = "../../output/randomStorage/randPlanes"
        for row in range(20000):
            yield np.loadtxt(planeFile, skiprows = row, max_rows = 1, usecols = np.arange(0, length))


def main():
    saveFig, param, test = readInputs()
    filter:sketchFilter = sketchFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()