"""
Runnable File containing a euclidean class:

python euclideanFilter.py --plot --test --param 1000,3,5

>> Still experimenting with params in this filter: numHashes,K,binWidth
TODO: Banding parameters are basically random. We need a better way to choose them.

Notes:
    python euclideanFilter.py --test --param 1000,2,10

    This filter has a huge degree of variance. Certain reads seem to never cluster together very nicely,
        except on certain runs with the same parameters where it all of sudden works.
        TODO Adding lines that don't just go through the origin would generate a farther 
            distance from most of the sequence vectors and perahaps reduce the variance that we are seeing.
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
    binWidth: The desired binWidth in a single random line for aligning.
    '''
    def __init__(self, param):
        params = param.split(",")
        self.numHashes, self.K, self.binWidth = int(params[0]), int(params[1]), float(params[2])
        self.title = f"euclideanFilter_numHashes={self.numHashes}_K={self.K}_binWidth={self.binWidth}"

    '''
    Preprocess all reads, required for some filters

    1. Build Characteristic Matrix Using Kmer Counts
    2. Apply projection approximation for Euclidean distance to Characteristic Matrix
    3. Amplify the metric using band method
        TODO: Formal analysis of the effect of AND/OR gates used here.
    '''
    def preprocessReads(self):
        self.m = np.power(4,self.K)
        self.projectionSignature()
        self.band()

    """
    Generate signatureMatrix of hashs x seqs
        A hash, in this case, is the projection of the seqs characteristicMatrix
        represenation onto a random line.
    """
    def projectionSignature(self):
        self.signatureMatrix = np.full((self.numHashes, self.n), np.inf)
        hashes = self.getRandomLines()
        for j in range(self.numHashes): # iterate lines
            for i in range(self.n): # iterate strings
                self.signatureMatrix[j,i] = self.projectAndBin(self.getCharacteristicVector(i), hashes[j])

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
    This function returns a random set of unit lines to use with projection as hash functions
        Only generate lines through the origin for now, but some amount of shift may be worth testing later
        - The smaller distances from the projected lines without it may cause most strings to appear too similar.
    Uniformly distributed lines on unit hypersphere: 
        https://math.stackexchange.com/questions/444700/uniform-distribution-on-the-surface-of-unit-sphere
    """
    def getRandomLines(self):
        dimensions = self.m
        lines = np.zeros((self.numHashes, dimensions))
        for i in range(self.numHashes):
            s = np.random.normal(0, 1, dimensions)
            lines[i,:] = s / np.sqrt(sum(s*s))
        return lines

    """
    Compute projection of point onto a random line and calculate which bin the point falls into.
    Returns pos or neg integer indicating which bin projected into. i.e.
    |-2|-1|0(contains origin at leftmost boundary)|1|2|
        Bins are computed as number of binwidths from the origin.
        - Because we are centered at the origin we can simply compare first dimensions.
        - A 0 value in our vector in the first dimension is rare enough we don't care
        - Use binWidth projected to first dimension as our comparison width
        Don't need to normalize projections because we use unit vectors
    """
    def projectAndBin(self, characteristicVector, line):
        projection = np.dot(line,characteristicVector)*line
        projectedWidth = line[0] * self.binWidth
        return int(projection[0] // projectedWidth)

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
        return self.adjacencyMatrix[i,j]

def main():
    saveFig, param, test = readInputs()
    filter:euclideanFilter = euclideanFilter(param)
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()