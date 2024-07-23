"""
Testing File for read filtering methods.
"""
from os.path import join
import os
from abc import ABC, abstractmethod
import sys
import subprocess
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from itertools import repeat
import time

from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import numpy as np

M = 24000

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
    def fill(self, seqs, labels, fastqFile, randLines = None, randOffsets = None, randPlanes = None, hashes = None):

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

        # fracMinHash stuff
        #s = 1e-3
        #maxVal = 0xFFFFFFFFFFFFFFFF # maximum hashable value on base 64 system.
        #self.limit = maxVal*s

        # Largest range of bed region in HG002 dataset
        # This is probably still too small as the multiple reads in that location will have different sequences, 
        #   leading to a higher number of total kmers
        #   Defined below too, don't assign one without the other.
        self.m = M

        # Assign each seen kmer a unique index
        self.kmer_dict = dict()
        self.next_kmer_index = 0

        # Load Euclidean random vectors
        self.randLines = randLines
        self.randOffsets = randOffsets

        # Load sketch random sketch planes
        self.randPlanes = randPlanes

        # Load random minHash functions
        self.hashes = hashes

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
        characteristicVector = np.zeros(self.minKmerRead, dtype=np.int16)
        kmerDir = f"../../output/kmerTables_20bp/{int(self.K)}mers" # convert to int to eliminate any trailing 0s
        chrom = os.path.basename(self.fastqFile).split("_")[0]
        regionDir = os.path.join(kmerDir, chrom, os.path.basename(self.fastqFile)[:-6])
        tableFile = os.path.join(regionDir, f"read{i+1}.ktab") # Files are 1-indexed

        # Some regions are too small, don't have haplotypes, or too few reads and led to errors. So table doesn't exist
        if not os.path.isfile(tableFile):
            with open("./tmp/bad_reads", "a") as f:
                f.write(tableFile + "\n")
            return characteristicVector
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
                if self.minKmerRead * 2 > self.m:
                    addValue = self.m - localIndex
                else: addValue = self.minKmerRead
                characteristicVector = np.concatenate((characteristicVector, np.zeros(addValue)), axis=0)
                self.minKmerRead *= 2
            characteristicVector[localIndex] = count

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
        os.makedirs(outDir, exist_ok = True)
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
def runFilter(filter : Filter, saveFig = False, test = False, verbose = True, inplace = True, 
              loadEuclidean = False, loadSketch = False, loadMinHash = False):
    if test:
        return testFilter(filter, saveFig, loadEuclidean = loadEuclidean, loadSketch = loadSketch, loadMinHash = loadMinHash)
    
    return runAllSamples(filter, saveFig, verbose, inplace, 
                         loadEuclidean = loadEuclidean, loadSketch = loadSketch, loadMinHash = loadMinHash)

'''
Filter Testing method that saves a heatmap Adjacency Matrix 
'''
def testFilter(filter : Filter, saveFig, loadEuclidean = False, loadSketch = False, loadMinHash = False):
    print(filter.title)

    m = M
    if loadEuclidean:
        lineFile = "../../output/randomStorage/randLines"
        randLines = np.loadtxt(
            lineFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = np.arange(0, m)
        )
        offsetFile = "../../output/randomStorage/randOffsets"
        randOffsets = np.loadtxt(
            offsetFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = 0
        )
    else:
        randLines = None
        randOffsets = None
    
    if loadSketch:
        planeFile = "../../output/randomStorage/randPlanes"
        randPlanes = np.loadtxt(
            planeFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = np.arange(0, m)
        )
    else:
        randPlanes = None

    if loadMinHash:
        numHashes = 2000
        large_prime = 7919
        hashes = ([0 for _ in range(numHashes)],[0 for _ in range(numHashes)], large_prime)
        for i in range(numHashes):
            hashes[0][i] = np.random.randint(1, large_prime - 1)
            hashes[1][i] = np.random.randint(0, large_prime - 1)
    else:
        hashes = None

    print("Loaded hashes")

    testDir = "../../output/HG002/chr8/"
    regions = [
        "chr8_593480-594653.fastq",
    ]

    totTruePos, totPos, totFalsePos, totNeg  = 0, 0, 0, 0

    for region in regions:
        fastqFile = join(testDir, region)

        tpr, fpr = runFilterOnFastq(filter, fastqFile, verbose = True, test = True, 
                                  inplace = False, 
                                  randLines = randLines, randOffsets = randOffsets, 
                                  randPlanes = randPlanes, 
                                  hashes = hashes
        )
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
def runAllSamples(filter, saveFig = False, verbose = True, inplace = True, 
                  loadEuclidean = False, loadSketch = False, loadMinHash = False):

    m = M
    if loadEuclidean:
        lineFile = "../../output/randomStorage/randLines"
        randLines = np.loadtxt(
            lineFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = np.arange(0, m)
        )
        offsetFile = "../../output/randomStorage/randOffsets"
        randOffsets = np.loadtxt(
            offsetFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = 0
        )
    else:
        randLines = None
        randOffsets = None

    if loadSketch:
        planeFile = "../../output/randomStorage/randPlanes"
        randPlanes = np.loadtxt(
            planeFile, 
            skiprows = 0, max_rows = filter.numHashes, 
            usecols = np.arange(0, m)
        )
    else:
        randPlanes = None

    # TODO this universal hash construction leads to some really weird results.
    if loadMinHash:
        numHashes = 2000
        large_prime = 7919
        num_buckets = 10000
        hashes = [0 for _ in range(numHashes)]
        for i in range(numHashes):
            a = np.random.randint(1, large_prime - 1)
            b = np.random.randint(0, large_prime - 1)
            hashes[i] = lambda x: ((a * x + b) % large_prime) % num_buckets
    else:
        hashes = None

    TruePos = TotPos = FalsePos = TotNeg = 0

    samplePath = "../../output/"
    for sample in [join(samplePath, "HG002_20bp")]:
        sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg = runAllRegions(filter, sample, verbose, 
                                                                                  inplace = inplace,
                                                                                  randLines = randLines, randOffsets = randOffsets,
                                                                                  randPlanes = randPlanes, 
                                                                                  hashes = hashes)
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
def runAllRegions(filter, directory, verbose, inplace = True, 
                  randLines = None, randOffsets = None, randPlanes = None, hashes = None):
    sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg  = 0, 0, 0, 0
    totalTime = 0
    for chrom in os.listdir(directory):

        # Build paths
        chromPath = join(directory, chrom)
        fastqFiles = [join(chromPath,x) for x in os.listdir(chromPath)]

        inputs = zip(repeat(filter), fastqFiles, repeat(False), repeat(inplace), repeat(randLines), 
                     repeat(randOffsets), repeat(randPlanes), repeat(hashes), repeat(False))

        start = time.time()
        print(f"Starting {chrom}\t", start)
        print(f"Num regions:\t{len(fastqFiles)}")
        num_cores = cpu_count()
        outputs = Parallel(n_jobs=num_cores, verbose=1)(delayed(runFilterOnFastq)(i, j, k, l, m, n, o, p, q) for \
                                                        i, j, k, l, m, n, o, p, q  in inputs)

        end = time.time()
        print(f"Completed {chrom}\t", end)
        totalTime += end-start
        print(f"{chrom} time (s):\t", end-start)

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

    print("Total compute time:\t", totalTime)
    return sampleTruePos, sampleTotPos, sampleFalsePos, sampleTotNeg

"""
Run initialized Filter on fastqFile input
inplace indicates to save the signature matrix

rand___ and hashes are all preloaded values reused across every thread and passed between filters. Only generated once.
"""
def runFilterOnFastq(filter, fastqFile, test = False, inplace = True, randLines = None, randOffsets = None, randPlanes = None, hashes = None, verbose = False):
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
    def runLoadedFilter(filter, saveFig, inplace = True, fastqFile = "", test = False):
        # Filter sequence set
        characteristicMatrix = filter.processReads(test = test)

        # Save characteristic matrix if inplace = "Output Directory for characteristic matricies"
        if type(inplace) == str:
            saveFile = join(inplace, os.path.basename(fastqFile.replace(".fastq","")))
            np.save(saveFile, np.insert(characteristicMatrix, 0, filter.labels, axis = 0))
            #print(saveFile)

        # Test filtered sequence set and haplotypes
        if saveFig:
            filter.plotAdjacencyMatrix()

        # Return a set of metrics
        return filter.tpr_fpr()

    seqs, labels = readFastq(fastqFile)
    if len(seqs) == 0:
        with open("./tmp/bad_regions", "a") as f:
            f.write(fastqFile + "\n")
        return "0:0", "0:0"
        raise Exception("Not enough seqs")

    filter.fill(seqs, labels, fastqFile, randLines = randLines, randOffsets = randOffsets, randPlanes = randPlanes, hashes = hashes)

    if verbose: print("Filter loaded, running filtering process")

    return runLoadedFilter(filter, saveFig = False, inplace = inplace, fastqFile = fastqFile, test = test)

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
    filters = [ 
    #    ("sketchFilter.py",     "10,20,4"),
        ("minHashFilter.py",    "1000,20,1"),
    #    ("euclideanFilter.py",  "10,20,5,10"),
    #    ("lengthFilter.py",     "10,1"),
    #    ("hashFilter.py",       "5")
    ]

    # Sketch Filter
    command = ["python", "___", "--test" ,"--param", "____", "--plot"]
    for filterFile, params in filters:
        print("Testing:", filterFile)
        command[1] = filterFile
        command[4] = params
        subprocess.run(command)

if __name__ == "__main__":
    main()
