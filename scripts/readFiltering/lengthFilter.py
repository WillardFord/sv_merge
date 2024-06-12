"""
Runnable File containing a lengthFilter class:

python lengthFilter.py [thresh]
"""
from test_filters import Filter, testFilter, readInputs
import sys


class lengthFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, seqs, labels, threshold):
        self.threshold = threshold
        super().__init__(seqs, labels)

    def __init__(self, threshold):
        self.threshold = threshold
        self.title = f"lengthFilter_threshold_{threshold}"

    '''
    Preprocess all reads, required for some filters
    '''
    def preprocessReads(self):
        pass

    '''
    Connect elements with dif(lengths) <= self.threshold.
    '''
    def connect(self, i, j):
        if abs(len(self.seqs[i]) - len(self.seqs[j])) > self.threshold:
            return 0
        return 1

def main():
    saveFig, lengthThreshold = readInputs()
    filter:lengthFilter = lengthFilter(lengthThreshold)
    testFilter(filter, saveFig = saveFig)

if __name__ == "__main__":
    main()