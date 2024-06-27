"""
Runnable File containing a sketch filter class:
    This approximates the cosine distance between two nodes.

python sketchFilter.py --plot --test --param 1000,3,40

>> Still experimenting with params in this filter: numHashes,K,bandLength
TODO: Banding parameters are basically random. We need a better way to choose them.
"""
import time
from test_filters import Filter, runFilter, readInputs
import numpy as np

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
        self.sketchSignature()
        print("Woot")
        self.band()

    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is --
    """
    def sketchSignature(self):
        self.signatureMatrix = np.zeros((self.numHashes, self.n))
        for i in range(self.n): # iterate strings
            for j, line in enumerate(self.getRandomPlanes()):
                print("Built random line\t\t\t", time.time())
                characteristicVector = self.getCharacteristicVector(i)
                print("Built characteristic vector\t\t", time.time())
                self.signatureMatrix[j,i] = self.sketch(characteristicVector, line)
                print("Completed one hash\t\t\t", time.time())
            print("Completed one read")
        print("Built signatureMatrix!")

    """
    This function returns a random set of binary vectors in [-1, 1] that describe a plane.

    TODO: takes ~ 6 seconds for k=20 and s=1e-3 --> dimension=(4^20)*s*1.05, 
        this needs to be basically instant.
    """
    def getRandomPlanes(self):
        seed = 1
        np.random.seed(seed)
        for _ in range(self.numHashes):
            orthonormal_line = np.random.choice([1, -1], size=self.m)
            yield orthonormal_line

    """
    Returns 1 if dot product is positive, else -1. (0's are randomly assigned)
    This approximates the cosine distance though is not exact.

    TODO: takes ~1.6 for k=20 and s=1e-3 --> dimension=(4^20)*s*1.05, 
        this needs to be basically instant
    """
    def sketch(self, characterisicVector, plane):
        value = np.dot(characterisicVector, plane)
        if value > 0:
            return 1
        if value < 0:
            return -1
        return np.random.choice([1, -1])

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