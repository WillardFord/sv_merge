"""
Testing File for read filtering methods.
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
from abc import ABC, abstractmethod
import sys

"""
Generic Abstract Class for filtering methods.
"""
class Filter(ABC):
    '''
    Load all class fields.

    seqs: length n list of sequences
    labels: length n list of corresponding sequence labels
    adjacencyMatrix: n by n matrix indicating if two sequences are sufficiently similar
    self.n
    self.firstOcc: First occurrence of each label in labels. Used during plotting and calculating metrics.
    self.title: Title for use in plots
    '''
    def fill(self, seqs, labels):
        seqs, labels = zip(*sorted(zip(seqs, labels), key = lambda x : x[1]))
        self.seqs = seqs
        self.labels = labels
        self.n :int = len(seqs)
        self.adjacencyMatrix: list[list[int]] = [[0 for _ in range(self.n)] for _ in range(self.n)]
        labelSet = set(self.labels)
        self.firstOcc = [0 for _ in range(len(labelSet))]
        for i, item in enumerate(sorted(labelSet)):
            self.firstOcc[i] = self.labels.index(item)
        try:
            print(self.title, file = sys.stderr)
        except:
            self.title = "Unspecified Filter"

    '''
    Preprocess all reads, required for some filters
    '''
    @abstractmethod
    def preprocessReads(self):
        pass

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
    TODO: indicate which reads come from which samples. 
    TODO: Add tpr_fpr labels in key??
    '''
    def plotAdjacencyMatrix(self, outDir = "../../output/plots"):
        plt.imshow(self.adjacencyMatrix, cmap = 'binary_r', interpolation='nearest')
        plt.title(f"Adjacency Matrix, {self.title}")
        plt.xticks(np.arange(self.n), np.arange(1, self.n + 1))
        plt.yticks(np.arange(self.n), np.arange(1, self.n + 1))

        plt.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)
        plt.xlabel(f"White in pos i,j indicates i connects j")
        for x in self.firstOcc[1:]:
            plt.axhline(x - 0.5, color='red', linestyle='-', linewidth=2)
            plt.axvline(x - 0.5, color='red', linestyle='-', linewidth=2)

        saveLocation = join(outDir, f"{self.title}.jpg")
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

        # Don't include y=x line.
        tpr = f"{int((np.sum(self.adjacencyMatrix, where = trueMask) - self.n)/2)}:{int((np.sum(trueMask) - self.n)/2)}"
        fpr = f"{int(np.sum(self.adjacencyMatrix, where = ~trueMask)/2)}:{int(np.sum(~trueMask)/2)}"

        return tpr, fpr

'''
Filter Testing method that saves a heatmap Adjacency Matrix and 
'''
def testFilter(filter : Filter, saveFig = False):
    # Load ground truth set of sequences and haplotypes
    # Randomly chosen region: output/chr1/HG002/chr1_676787-678886.fastq

    region = ["chr1", 676787, 678886]
    seqs, labels = loadSamplesSeqs(region)
    filter.fill(seqs, labels)

    # Filter sequence set
    filter.preprocessReads()
    filter.buildAdjacencyMatrix()

    # Test filtered sequence set and haplotypes
    if saveFig:
        filter.plotAdjacencyMatrix()

    # Return a set of metrics
    tpr, fpr = filter.tpr_fpr()
    print(f"{tpr}\t{fpr}")

def loadSamplesSeqs(region:str):
    chrom, start, end = region
    seqs = []
    labels = []
    samples = ["HG002", "HG00733"]
    dataDir = join("../../output", chrom)
    for sample in samples:
        fastqFile = join(dataDir, sample, f"{chrom}_{start}-{end}.fastq")
        sampleSeqs, sampleLabels = readFastq(fastqFile)
        sampleLabels = (np.array(sampleLabels, dtype = np.uint) + len(set(labels))).tolist() # Different samples have different haplotypes.
        seqs += sampleSeqs
        labels += sampleLabels
    return seqs, labels

def readFastq(fastqFile:str = None):
    if fastqFile == None:
        fastqFile = "/Users/wford/Documents/sv_merge/output/test_extract/chr10_756193-756593.fastq"

    seqs :list[str] = []
    labels :list[str] = []
    with open(fastqFile,"r") as f:
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

def readInputs():
    saveFig = False
    param = 10
    for i, item in enumerate(sys.argv):
        if item == "--param":
            param = int(sys.argv[i+1])
        elif item == "--plot":
            saveFig = bool(sys.argv[i+1])
    return saveFig, param

def main():
    print("\nWhat are you doing here?\n\nThis is a resources file called by other filters.")

if __name__ == "__main__":
    main()
