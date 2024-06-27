"""
Testing File for read filtering methods.
"""
from os.path import join
import os
from abc import ABC, abstractmethod
import sys
import subprocess
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

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
    def fill(self, seqs, labels, fastqFile):
        seqs, labels = zip(*sorted(zip(seqs, labels), key = lambda x : x[1]))
        self.seqs = seqs
        self.labels = labels
        self.n :int = len(seqs)
        self.adjacencyMatrix = np.zeros((self.n, self.n), bool)
        labelSet = set(self.labels)
        self.firstOcc = [0 for _ in range(len(labelSet))]
        for i, item in enumerate(sorted(labelSet)):
            self.firstOcc[i] = self.labels.index(item)
        self.fastqFile = fastqFile

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

    """
    Return the characteristic vector for a single read, i
    """
    def getCharacteristicVector(self, i):
        """
        Generate all possible kmers without storing them in memory.
        """
        def generateKmers(k:int):
            alphabet = "ACGT"
            if k == 1:
                for char in alphabet:
                    yield char
            else:
                for mer in generateKmers(k-1):
                    for char in alphabet:
                        yield mer + char

        characteristicVector = np.zeros(self.m, dtype=np.int16)

        kmerDir = f"/Users/wford/Documents/sv_merge/output/kmerTables/{int(self.K)}mers"
        chrom = os.path.basename(self.fastqFile).split("_")[0]
        regionDir = os.path.join(kmerDir, chrom, os.path.basename(self.fastqFile)[:-6])
        tableFile = os.path.join(regionDir, f"read{i+1}.ktab") # Files are 1-indexed
        #print(tableFile)

        # Some regions are too small, don't have haplotypes, or too few reads and led to errors. So table doesn't exist
        if not os.path.isfile(tableFile):
            raise FileNotFoundError(tableFile)

        command = ["tabex", tableFile, "LIST"]
        table_output = subprocess.check_output(command)
        table_output_str = str(table_output)
        del table_output

        # Base 4 to base 10 conversion
        alphabet = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
        seq_to_index = lambda s: sum(alphabet[c] * (4 ** p) for p, c in enumerate(s[::-1]))

        for row in  table_output_str.split("\\n")[1:]:
            kmer, count = row.strip().split("=")
            count = int(count.strip())
            kmer = kmer.split(":")[1].strip()
            characteristicVector[seq_to_index(kmer)] = count

        return characteristicVector

    """
    Use banding technique to bin reads by similarity
    """
    def band(self):
        def connectBucket(bucket):
            for i in range(len(bucket)-1):
                for j in range(i+1, len(bucket)):
                    self.adjacencyMatrix[bucket[i],bucket[j]] = True
                    self.adjacencyMatrix[bucket[j],bucket[i]] = True

        for i in range(self.signatureMatrix.shape[0]//self.bandLength):
            buckets = defaultdict(list[int])
            startBucket = i*self.bandLength
            endBucket = (i+1)*self.bandLength
            for j in range(self.n):
                buckets[self.signatureMatrix[startBucket:endBucket, j].tobytes()].append(j)
            for bucket in buckets.values():
                connectBucket(bucket)

    '''
    Plot adjacencyMatrix
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

"""
General Callable Function.
    Calls testFiler or runAllSamples
"""
def runFilter(filter : Filter, saveFig = False, test = False, verbose = True):
    if test:
        return testFilter(filter, saveFig)
    
    return runAllSamples(filter, saveFig, verbose)

'''
Filter Testing method that saves a heatmap Adjacency Matrix 
'''
def testFilter(filter : Filter, saveFig):
    print(filter.title)

    testDir = "/Users/wford/Documents/sv_merge/output/test_extract/chr10"
    regions = [
        "chr10_756193-756593.fastq",
        #"chr10_764127-765357.fastq",
        #"chr10_774190-774696.fastq",
        #"chr10_775244-776398.fastq",
        #"chr10_776350-776890.fastq",
        #"chr10_777367-779252.fastq",
        #"chr10_781668-782528.fastq",
        #"chr10_787900-788394.fastq",
        #"chr10_790342-790742.fastq"
    ]

    totTruePos, totPos, totFalsePos, totNeg  = 0, 0, 0, 0

    for region in regions:
        fastqFile = join(testDir, region)

        tpr,fpr = runFilterOnFastq(filter, fastqFile)
        if tpr == fpr and fpr == 0:
            continue

        truePos, regionTotPos = tpr.split(":")
        falsePos, regionTotNeg = fpr.split(":")

        totTruePos += int(truePos)
        totPos += int(regionTotPos)
        totFalsePos += int(falsePos)
        totNeg += int(regionTotNeg)
    
    print(f"{totTruePos}:{totPos}\t{totFalsePos}:{totNeg}")
    print(f"TPR:\t{totTruePos/totPos}\nFPR:\t{totFalsePos/totNeg}" ,file = sys.stderr)

    return totTruePos, totPos, totFalsePos, totNeg

"""
Load filter on all samples
"""
def runAllSamples(filter, saveFig = False, verbose = True):
    """
    Test filter on all regions for a single sample and return the total number of false positives, 
        maximum possible, and equivalent for num false negatives
    """
    def runAllRegions(filter, saveFig, directory, verbose):
        badRegions = open("./tmp/badRegions.txt", "w+")

        sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg  = 0, 0, 0, 0
        for chrom in os.listdir(directory):
            for fastqFile in os.listdir(join(directory, chrom)):
                fastqPath = join(directory, chrom, fastqFile)
                try:
                    tpr,fpr = runFilterOnFastq(filter, fastqPath)
                except Exception as e:
                    badRegions.write(fastqPath + "\n")
                    continue

                truePos, regionTotPos = tpr.split(":")
                falsePos, regionTotNeg = fpr.split(":")

                sampleTruePos += int(truePos)
                sampleTotPos += int(regionTotPos)
                sampleFalsePos += int(falsePos)
                sampleTotNeg += int(regionTotNeg)
                if verbose:
                    print(f"{sampleTruePos}:{sampleTotPos}\t{sampleFalsePos}:{sampleTotNeg}", file = sys.stderr)
            if verbose:
                print(f"Completed\t{join(directory, chrom)}", file=sys.stderr)
        badRegions.close()
        return sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg

    samplePath = "../../output/"
    TruePos, TotPos, FalsePos, TotNeg = 0, 0, 0, 0
    for directory in [x for x in os.listdir(samplePath) if x.startswith("HG")]:
        if "733" in directory: continue # Data has duplicated reads so skip for now. TODO

        sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg = runAllRegions(filter, saveFig, 
                                                                                  join(samplePath, directory),
                                                                                  verbose)

        TruePos += sampleTruePos
        TotPos += sampleTotPos
        FalsePos += sampleFalsePos
        TotNeg += sampleTotNeg
        print(f"Completed\t{directory}", file=sys.stderr)

    if verbose:
        print(f"{TruePos}:{TotPos}\t{FalsePos}:{TotNeg}")
        print(f"TPR:\t{TruePos/TotPos}\nFPR:\t{FalsePos/TotNeg}" ,file = sys.stderr)

    return TruePos, TotPos, FalsePos, TotNeg

"""
Run initialized Filter on fastqFile input
"""
def runFilterOnFastq(filter, fastqFile):
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

    seqs, labels = readFastq(fastqFile)

    if len(seqs) == 0:
        raise Exception("Not enough seqs")

    filter.fill(seqs, labels, fastqFile)
    return runLoadedFilter(filter, saveFig = False)

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
    print("Testing all other filters using default parameters.")
    filters = [ ("sketchFilter.py",     "1000,3"),
                ("minHashFilter.py",    "1000,12"),
                ("euclideanFilter.py",  "1000,3,5"),
                ("lengthFilter.py",     "10"),
                ("hashFilter.py",       "5")
    ]

    # Sketch Filter
    command = ["python", "", "--test" ,"--param", ""]
    for filterFile, params in filters:
        command[1] = filterFile
        command[4] = params
        subprocess.run(command)

if __name__ == "__main__":
    main()
