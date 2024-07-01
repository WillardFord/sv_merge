"""
Runnable File containing a hashFilter class:

python hashFilter.py --param 0.05
"""

from test_filters import Filter, runFilter, readInputs
from collections import defaultdict 
import numpy as np

class hashFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, percentHashed):
        self.percentHashed = float(percentHashed)
        self.title = f"hashFilter_percentHashed_{percentHashed}"

    '''
    Preprocess all reads, required for some filters
    '''
    def processReads(self):
        def connectBucket(bucket):
            for i in range(len(bucket)-1):
                for j in range(i+1, len(bucket)):
                    self.adjacencyMatrix[bucket[i],bucket[j]] = True
                    self.adjacencyMatrix[bucket[j],bucket[i]] = True

        buckets = defaultdict(list[int])

        l = self.getMinLength()
        freq = int(l * self.percentHashed)
        indicies = np.random.randint(0, high=(l-1), size=freq)

        for i, seq in enumerate(self.seqs):
            hashSeq = np.frombuffer(seq.encode(), dtype=np.uint8)[indicies].tobytes()
            buckets[hashSeq].append(i)

        for bucket in buckets.values():
                connectBucket(bucket)

    def getMinLength(self):
        n = np.inf
        for seq in self.seqs:
            if len(seq) < n:
                n = len(seq)
        return n

def main():
    saveFig, percentHashed, test = readInputs()
    filter:hashFilter = hashFilter(float(percentHashed))
    runFilter(filter, saveFig, test)

if __name__ == "__main__":
    main()