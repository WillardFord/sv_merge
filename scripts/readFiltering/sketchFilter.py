"""
Runnable File containing a sketchFilter class:

TODO wip
python sketchFilter.py ____
"""

from test_filters import Filter, testFilter, readInputs
from collections import defaultdict 
import numpy as np
import sys

class sketchFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, param):
        self.param = param
        self.title = f"sketchFilter_param_{param}"
        self.buckets = defaultdict(list[int])

    '''
    Preprocess all reads, required for some filters
    '''
    def preprocessReads(self):
        pass

    '''
    Connect elements in the same bucket.
    '''
    def connect(self, i, j):
        pass

def main():
    saveFig, param = readInputs()
    filter:sketchFilter = sketchFilter(param)
    testFilter(filter, saveFig = saveFig)

if __name__ == "__main__":
    main()