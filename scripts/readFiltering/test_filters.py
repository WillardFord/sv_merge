"""
Testing File for read filtering methods.
"""
from os.path import join
import os
from abc import ABC, abstractmethod
import sys
import subprocess
from collections import defaultdict
from multiprocessing import Pool
from itertools import repeat

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
        self.adjacencyMatrix = np.eye(self.n, dtype=bool)
        labelSet = set(self.labels)
        self.firstOcc = [0 for _ in range(len(labelSet))]
        for i, item in enumerate(sorted(labelSet)):
            self.firstOcc[i] = self.labels.index(item)
        self.fastqFile = fastqFile
        self.minKmerRead = 750

        # Don't need 
        #s = 1e-3
        #maxVal = 0xFFFFFFFFFFFFFFFF # maximum hashable value on base 64 system.
        #self.limit = maxVal*s

        # Largest range of bed region in HG002 dataset
        self.m = 14156

        # Assign each seen kmer a unique index
        self.kmer_dict = dict()
        self.next_kmer_index = 0

        self.characteristicVectors = []

    '''
    Preprocess all reads, required for some filters
    '''
    @abstractmethod
    def processReads(self):
        pass

    """
    Return the characteristic vector for a single read, i
    """
    def getCharacteristicVector(self, i):
        # We can store all counts vectors in memory for a given region, may not be useful in practice
        if len(self.characteristicVectors) > i:
            if len(self.characteristicVectors[i] < self.minKmerRead):
                # Extend vector if needed
                self.characteristicVectors[i] = np.concatenate(
                    (
                        self.characteristicVectors[i],
                        np.zeros(
                            self.minKmerRead - len(self.characteristicVectors[i])
                        )
                    ),
                    axis=0
                )
            return self.characteristicVectors[i]

        characteristicVector = np.zeros(self.minKmerRead, dtype=np.int16)
        kmerDir = f"../../output/kmerTables/{int(self.K)}mers" # convert to int to eliminate any trailing 0s
        chrom = os.path.basename(self.fastqFile).split("_")[0]
        regionDir = os.path.join(kmerDir, chrom, os.path.basename(self.fastqFile)[:-6])
        tableFile = os.path.join(regionDir, f"read{i+1}.ktab") # Files are 1-indexed

        # Some regions are too small, don't have haplotypes, or too few reads and led to errors. So table doesn't exist
        if not os.path.isfile(tableFile):
            raise FileNotFoundError(tableFile)

        command = ["Tabex", tableFile, "LIST"]
        table_output = str(subprocess.check_output(command))

        for row in  table_output.split("\\n")[1:-1]: # First and last line don't contain kmers
            kmer, count = row.strip().split("=")
            kmer = kmer.split(":")[1].strip().upper()
            #if abs(hash(kmer)) < self.limit: Don't need frac min hash if lazily generating vectors
            count = int(count.strip())

            # Get kmer index
            if kmer not in self.kmer_dict.keys():
                self.kmer_dict[kmer] = self.next_kmer_index
                self.next_kmer_index += 1
            
            localIndex = self.kmer_dict[kmer]
            if localIndex >= self.minKmerRead:
                characteristicVector = np.concatenate((characteristicVector, np.zeros(self.minKmerRead)), axis=0)
                self.minKmerRead *= 2
            characteristicVector[localIndex] = count

        self.characteristicVectors.append(characteristicVector)
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

        for i in range(self.numHashes//self.bandLength):
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
def runFilter(filter : Filter, saveFig = False, test = False, verbose = True, inplace = True):
    if test:
        return testFilter(filter, saveFig)
    
    return runAllSamples(filter, saveFig, verbose, inplace)

'''
Filter Testing method that saves a heatmap Adjacency Matrix 
'''
def testFilter(filter : Filter, saveFig):
    print(filter.title)

    testDir = "../../output/HG002/chr8/"
    regions = [
        "chr8_593480-594653.fastq",
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
def runAllSamples(filter, saveFig = False, verbose = True, inplace = True):
    samplePath = "../../output/"

    for sample in [join(samplePath, x) for x in os.listdir(samplePath) if x.startswith("HG") and "733" not in x]:
        sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg = runAllRegions(filter, sample, verbose, inplace)
        TruePos += sampleTruePos
        TotPos += sampleTotPos
        FalsePos += sampleFalsePos
        TotNeg += sampleTotNeg

    if verbose:
        print(f"{TruePos}:{TotPos}\t{FalsePos}:{TotNeg}")
        print(f"TPR:\t{TruePos/TotPos}\nFPR:\t{FalsePos/TotNeg}" ,file = sys.stderr)

    return TruePos, TotPos, FalsePos, TotNeg

"""
Function to run across all directories in parallel
"""
def runAllRegions(filter, directory, verbose, inplace = True):
    sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg  = 0, 0, 0, 0
    for chrom in os.listdir(directory):
        chromPath = join(directory, chrom)
        fastqFiles = [join(chromPath,x) for x in os.listdir(chromPath)]
        inputs = zip(repeat(filter), fastqFiles, repeat(inplace))
        numCores = 4
        with Pool(numCores) as p:
            outputs = p.starmap(runFilterOnFastq, inputs)

        for tpr, fpr in outputs:
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
    return sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg

"""
Run initialized Filter on fastqFile input
"""
def runFilterOnFastq(filter, fastqFile, inplace = True):
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
    def runLoadedFilter(filter, saveFig, inplace = True, fastqFile = ""):
        # Filter sequence set
        print(fastqFile)
        characteristicMatrix = filter.processReads()

        # Save characteristic matrix if inplace = "Output Directory for characteristic matricies"
        if type(inplace) == str:
            print("Saving characteristic matrix")
            saveFile = join(inplace, fastqFile.replace(".fastq",""))
            np.save(saveFile, characteristicMatrix)
        print("Did this save?")

        # Test filtered sequence set and haplotypes
        if saveFig:
            filter.plotAdjacencyMatrix()

        # Return a set of metrics
        return filter.tpr_fpr()

    seqs, labels = readFastq(fastqFile)

    if len(seqs) == 0:
        raise Exception("Not enough seqs")

    filter.fill(seqs, labels, fastqFile)
    return runLoadedFilter(filter, saveFig = False, inplace = inplace, fastqFile = fastqFile)

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
    print("Testing all filters using default parameters.")
    filters = [ ("sketchFilter.py",     "10,20,4"),
                ("minHashFilter.py",    "10,20,10"),
                ("euclideanFilter.py",  "10,20,5,10"),
                ("lengthFilter.py",     "10,1"),
                ("hashFilter.py",       "5")
    ]

    # Sketch Filter
    command = ["python", "", "--test" ,"--param", ""]
    for filterFile, params in filters:
        print("Testing:", filterFile)
        command[1] = filterFile
        command[4] = params
        subprocess.run(command)

if __name__ == "__main__":
    main()
