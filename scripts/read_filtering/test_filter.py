"""
Testing File for read filtering methods.
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
from abc import ABC, abstractmethod

"""
Generic Abstract Class for filtering methods.
"""
class Filter(ABC):
    '''
    seqs: length n list of sequences
    labels: length n list of corresponding sequence labels
    adjacencyMatrix: n by n matrix indicating if two sequences are sufficiently similar
    self.n
    self.firstOcc: First occurrence of each label in labels. Used during plotting and calculating metrics.
    '''
    def __init__(self, seqs, labels):
        self.fill(seqs, labels)

    '''
    Load all class fields
    '''
    def fill(self, seqs, labels):
        try:
            print("Seqs already exists:", self.seqs)
        except:
            self.seqs = seqs
            self.labels = labels
            self.n :int = len(seqs)
            self.adjacencyMatrix: list[list[int]] = [[0 for _ in range(self.n)] for _ in range(self.n)]
            labelSet = set(self.labels)
            self.firstOcc = [0 for _ in range(len(labelSet))]
            for i, item in enumerate(sorted(labelSet)):
                self.firstOcc[i] = self.labels.index(item)

    '''
    Fill out classes based on class filtering method
    '''
    @abstractmethod
    def connect(self, i, j):
        pass

    '''
    Build adjacencyMatrix
    '''
    def buildAdjacencyMatrix(self):
        for i in range(self.n):
            for j in range(self.n):
                if i == j:
                    self.adjacencyMatrix[i][j] = 1
                    continue
                self.adjacencyMatrix[i][j] = self.connect(i,j)

    '''
    Plot adjacencyMatrix
    '''
    def plotAdjacencyMatrix(self, filterName, seqsName, outDir = "../../output/plots"):
        plt.imshow(self.adjacencyMatrix, cmap = 'binary_r', interpolation='nearest')
        plt.title(f"Adjacency Matrix, {filterName} on {seqsName}")
        plt.xticks(np.arange(self.n), np.arange(1, self.n + 1))
        plt.yticks(np.arange(self.n), np.arange(1, self.n + 1))

        plt.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)
        plt.xlabel(f"White in pos i,j indicates i connects j under {filterName}")
        for x in self.firstOcc[1:]:
            plt.axhline(x - 0.5, color='red', linestyle='-', linewidth=2)
            plt.axvline(x - 0.5, color='red', linestyle='-', linewidth=2)

        saveLocation = join(outDir, f"{filterName}_{seqsName}.jpg")
        plt.savefig(saveLocation) 

    '''
    Return tuple (true positive rate, false positive rate)
    tpr - sum of all connections in the clusters where connections should exist over that total space
    fpr - equivalently for connections that shouldn't exist
    '''
    def tpr_fpr(self):
        trueMask = np.zeros((self.n, self.n), dtype = bool)
        for i in range(0, len(self.firstOcc)-1):
            trueMask[self.firstOcc[i]:self.firstOcc[i+1], self.firstOcc[i]:self.firstOcc[i+1]] = 1
        # Add Last Region
        trueMask[self.firstOcc[-1]:, self.firstOcc[-1]:] = 1

        tpr = (np.sum(self.adjacencyMatrix, where = trueMask) - self.n) / (np.sum(trueMask) - self.n) # Don't include y=x line.
        fpr = np.sum(self.adjacencyMatrix, where = ~trueMask) / np.sum(~trueMask)

        return tpr, fpr

'''
Filter Testing method that prints out a heatmap Adjecency Matrix and 
'''
def testFilter(filter : Filter):
    # Load ground truth set of sequences and haplotypes
    # Seqs and labels should be ordered by labels for nice plotting and accurate metric counting.
    seqs, labels = readInSeqs()
    seqs, labels = zip(*sorted(zip(seqs, labels), key = lambda x : x[1]))
    seqsName = "chr10_756193-756593"

    filter.fill(seqs, labels)

    # Filter sequence set
    filter.buildAdjacencyMatrix()

    # Test filtered sequence set and haplotypes
    filterName = filter.__class__.__name__
    filter.plotAdjacencyMatrix(filterName, seqsName)

    # Return a set of metrics
    scores = filter.tpr_fpr()
    print(scores)

def readInSeqs(seqPath:str = None):
    if seqPath == None:
        seqPath = "/Users/wford/Documents/sv_merge/output/test_extract/chr10_756193-756593.fastq"

    seqs :list[str] = []
    labels :list[str] = []
    with open(seqPath,"r") as f:
        i = 0
        curHap = -1
        curSeq = ""
        for line in f.readlines():
            if i % 4 == 0:
                for segment in line.strip().split(" "):
                    if segment.startswith("HP"):
                        curHap = segment.split(":")[-1]
            elif i % 4 == 1:
                curSeq = line.strip()
                seqs.append(curSeq)
                labels.append(curHap)
            i += 1
    return seqs, labels

def main():
    print("\nWhat are you doing here?\n\nThis is a resources file called by other filters.")

if __name__ == "__main__":
    main()
