from test_filter import Filter, testFilter
from abc import ABC

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

    '''
    Fill out classes based on size difference.
    '''
    def connect(self, i, j):
        if abs(len(self.seqs[i]) - len(self.seqs[j])) > self.threshold:
            return 0
        return 1

def main():
    lengthThreshold = 2
    filter:lengthFilter = lengthFilter(lengthThreshold)
    testFilter(filter)

if __name__ == "__main__":
    main()