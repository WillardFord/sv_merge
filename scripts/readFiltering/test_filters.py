"""
Testing File for read filtering methods.
"""
import matplotlib.pyplot as plt
import numpy as np
from os.path import join
import os
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
        self.adjacencyMatrix: list[list[bool]] = [[0 for _ in range(self.n)] for _ in range(self.n)]
        labelSet = set(self.labels)
        self.firstOcc = [0 for _ in range(len(labelSet))]
        for i, item in enumerate(sorted(labelSet)):
            self.firstOcc[i] = self.labels.index(item)

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
                    self.adjacencyMatrix[i][j] = True
                    continue
                self.adjacencyMatrix[i][j] = self.connect(i,j)

    '''
    Plot adjacencyMatrix
    TODO: indicate which reads come from which samples. 
    TODO: Add tpr_fpr labels in key??
    '''
    def plotAdjacencyMatrix(self, outDir = "../../output/plots"):
        plt.imshow(self.adjacencyMatrix, 
                   cmap = 'binary_r', 
                   interpolation='nearest', 
                   vmin=0, vmax=1
        )
        plt.title(f"Adjacency Matrix, {self.title}")
        plt.xticks(np.arange(self.n), np.arange(1, self.n + 1))
        plt.yticks(np.arange(self.n), np.arange(1, self.n + 1))

        plt.tick_params(axis='x', top=True, bottom=False, labeltop=True, labelbottom=False)
        plt.xlabel(f"White in pos i,j indicates i connects j")
        for x in self.firstOcc[1:]:
            plt.axhline(x - 0.5, color='red', linestyle='-', linewidth=2)
            plt.axvline(x - 0.5, color='red', linestyle='-', linewidth=2)

        saveLocation = join(outDir, f"{self.title}.jpg")
        print(f"Saving at {saveLocation}", file = sys.stderr)
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
def testFilter(filter : Filter, saveFig):
    print(filter.title)

    fastqFile = "/Users/wford/Documents/sv_merge/output/test_extract/chr10_756193-756593.fastq"
    seqs, labels = readFastq(fastqFile)

    filter.fill(seqs, labels)

    tpr, fpr = runLoadedFilter(filter, saveFig)
    print(f"{tpr}\t{fpr}")

"""
Get Results from loaded filter
"""
def runLoadedFilter(filter, saveFig):
    # Filter sequence set
    filter.preprocessReads()
    filter.buildAdjacencyMatrix()

    # Test filtered sequence set and haplotypes
    if saveFig:
        filter.plotAdjacencyMatrix()

    # Return a set of metrics
    return filter.tpr_fpr()

"""
General Callable Function to handle average results and test runs.
"""
def runFilter(filter : Filter, saveFig, test):
    if test:
        return testFilter(filter, saveFig)
    
    runAllSamples(filter, saveFig)

"""
Load filter on all samples
"""
def runAllSamples(filter, saveFig):
    samplePath = "../../output/"
    TruePos, TotPos, FalsePos, TotNeg = 0, 0, 0, 0
    for directory in [x for x in os.listdir(samplePath) if x.startswith("HG")]:
        sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg = runAllRegions(filter, saveFig, join(samplePath, directory))
        TruePos += sampleTruePos
        TotPos += sampleTotPos
        FalsePos += sampleFalsePos
        TotNeg += sampleTotNeg
        print(f"Completed\t{directory}")

    print(f"{TruePos}:{TotPos}\t{FalsePos}:{TotNeg}")
    print(f"TPR:\t{TruePos/TotPos}\nFPR:\t{FalsePos/TotNeg}" ,file = sys.stderr)

"""
Test filter on all regions for a single sample and return the total number of false positives, 
    maximum possible, and equivalent for num false negatives

    minSeqsThresh is the minimum required number of phased sequences for us to use the region.
    Arbitrarily set to 10.
"""
def runAllRegions(filter, saveFig, directory):
    minSeqsThresh = 10
    sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg  = 0, 0, 0, 0
    for chrom in os.listdir(directory):
        for fastqFile in os.listdir(join(directory, chrom)):
            seqs, labels = readFastq(join(directory, chrom, fastqFile))
            if len(seqs) < minSeqsThresh:
                continue
            filter.fill(seqs, labels)
            tpr, fpr = runLoadedFilter(filter, saveFig)

            truePos, regionTotPos = tpr.split(":")
            falsePos, regionTotNeg = fpr.split(":")

            sampleTruePos += int(truePos)
            sampleTotPos += int(regionTotPos)
            sampleFalsePos += int(falsePos)
            sampleTotNeg += int(regionTotNeg)
        print(f"Completed\t{join(directory, chrom)}")

    return sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg

"""
Read in reads and labels form a phased fastq file
"""
def readFastq(fastqFile:str):

    seqs :list[str] = []
    labels :list[str] = []
    with open(fastqFile,"r") as f:
        i = 0
        curHap = -1
        for line in f.readlines():
            if i % 4 == 0:
                for segment in line.strip().split(" "):
                    if segment.startswith("HP"):
                        curHap = segment.split(":")[-1]
            elif i % 4 == 1:
                if curHap == -1:
                    continue
                curSeq = line.strip()
                seqs.append(curSeq)
                labels.append(curHap)
                curHap = -1
                curSeq = ""
            i += 1
    return seqs, labels

"""
Read in parameters for loading a filter
"""
def readInputs():
    saveFig = False
    param = 10
    test = False
    for i, item in enumerate(sys.argv):
        if item == "--param":
            param = sys.argv[i+1]
        elif item == "--plot":
            saveFig = True
        elif item == "--test":
            test = True
    return saveFig, param, test

def main():
    print("\nWhat are you doing here?\n\nThis is a resources file called by other filters.")

if __name__ == "__main__":
    main()
