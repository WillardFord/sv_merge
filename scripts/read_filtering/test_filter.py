"""
Testing File for read filtering methods.

"""
from collections import defaultdict
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
    classes: n by n matrix indicating if two sequences are sufficiently similar
    '''
    def __init__(self, seqs, labels):
        self.seqs:list[str] = seqs
        self.labels :list[int] = labels
        self.n :int = len(seqs)
        self.adjacencyMatrix: list[list[int]] = [[0] * self.n]*self.n

    def fill(self, seqs, labels):
        try:
            print(self.seqs)
        except:
            self.seqs = seqs
            self.labels = labels
            self.n :int = len(seqs)
            self.adjacencyMatrix: list[list[int]] = [[0] * self.n]*self.n
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
    def plotAdjacencyMatrix(self, filterName, seqsName, outDir = "plots"):
        print(self.adjacencyMatrix)

        colors = plt.cm.tab10.colors[:len(np.unique(self.labels))]
        cmap = plt.cm.colors.ListedColormap(colors)

        plt.imshow(self.adjacencyMatrix, cmap=cmap)
        plt.title(f"Adjacency Matrix, {filterName} on {seqsName}")
        #plt.axis('off')

        saveLocation = join(outDir, f"{filterName}_{seqsName}.jpg")
        plt.savefig(saveLocation) 
        plt.show()
        

    '''
    Return tuple (true positive rate, false positive rate)
    Note: If 1 group is split into 2 clusters, whichever cluster appears first 
        will count as tpr and the rest will be counted in fpr
    This should change. TODO
    '''
    def tpr_fpr(self):
        tpc = 0
        fpc = 0
        classes = defaultdict(list)
        for i, row in enumerate(self.adjacencyMatrix):
            classes[tuple(row)].append(i) # This only works if the connection is mutual, i.e. x~y -> y~x
        clusterLabels = [-1]*self.n
        for i, cluster in enumerate(classes.values()):
            for j, item in enumerate(cluster):
                if j == 0 and clusterLabels[self.labels[item]] == -1:
                    clusterLabels[self.labels[item]] = i
                if clusterLabels[self.labels[item]] == i:
                    tpc += 1
                else:
                    fpc += 1
        return tpc/self.n , fpc/self.n

'''
Filter Testing method that prints out a heatmap Adjecency Matrix and 
'''
def testFilter(filter : Filter):
    # Load ground truth set of sequences and haplotypes
    # Seqs and labels should be ordered by labels
    seqs, labels = readInSeqs()

    seqs, labels = zip(*sorted(zip(seqs, labels), key = lambda x : x[1]))

    seqsName = "chr10_756193-756593"

    filter.fill(seqs, labels)

    # Filter sequence set
    filter.buildAdjacencyMatrix()

    filterName = filter.__class__.__name__

    # Test filtered sequence set and haplotypes
    #scores = filter.tpr_fpr()
    # TODO add outDir for saving files somewhere, needed during testing.
    filter.plotAdjacencyMatrix(filterName, seqsName)
    print("yay")
    # Return a set of metrics
    #print(scores)


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
    print("\nWhat are you doing here?\n\nThis is a resources file called by other filtering files.")

if __name__ == "__main__":
    main()
