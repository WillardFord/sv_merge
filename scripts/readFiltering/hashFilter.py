"""
Runnable File containing a hashFilter class:

python hashFilter.py [thresh]
"""

from test_filters import Filter, testFilter, readInputs
from collections import defaultdict 
import numpy as np
import sys

class hashFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, percentHashed):
        self.percentHashed = percentHashed
        self.title = f"hashFilter_percentHashed_{percentHashed}"
        self.buckets = defaultdict(list[int])

    '''
    Preprocess all reads, required for some filters
    '''
    def preprocessReads(self):
        # TODO, this should be min length of all seqs. 
        # Can be estimated by length of first minus threshold from length filter
        n = len(self.seqs[0])
        freq = int(n * self.percentHashed // 100) # Lower Bound
        if freq < 1:
            print("Percent too low, using 1 base.")
            freq = 1
            self.title = f"hashFilter_percentHashed_1base"
        indicies = np.random.randint(0, high=(n-1), size=freq)
        for i, seq in enumerate(self.seqs):
            hashSeq = np.frombuffer(seq.encode(), dtype=np.uint8)[indicies].tobytes()
            self.buckets[hashSeq].append(i)
        
        def swap_keys_values(dictionary):
            new_dictionary = {}
            for key, values in dictionary.items():
                for value in values:
                    new_dictionary[value] = key
            return new_dictionary
        
        self.newLabels = swap_keys_values(self.buckets)

    '''
    Connect elements in the same bucket.
    '''
    def connect(self, i, j):
        return self.newLabels[i] == self.newLabels[j]

def main():
    saveFig, percentHashed = readInputs()
    filter:hashFilter = hashFilter(percentHashed)
    testFilter(filter, saveFig = saveFig)

if __name__ == "__main__":
    main()