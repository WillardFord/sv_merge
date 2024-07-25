'''
Generate a ROC curve for each filter across all data, assuming that we are using default arguments.
Generating the curve over number of filters required to test similarity.

TODO: We'd rather save the characteristic matrix for each filter given some parameters so that we can test several different combination methods at once.
    - minHash
    - euclidean
    - sketch
'''

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import csv

from test_filters import runFilter
from lengthFilter import lengthFilter
from euclideanFilter import euclideanFilter
from sketchFilter import sketchFilter
from minHashFilter import minHashFilter

"""
Calculate area under curve using middle estimate Riemann sum.
"""
def auc(fprs, tprs):
    n = len(tprs)
    area = 0
    for i in range(1, n):
        area += (fprs[i]-fprs[i-1])*((tprs[i]+tprs[i-1])/2)
    return area

"""
Plot the ROC curve with area under curve labeled.
"""
def plotROC(tprs, fprs, location):
    roc_auc = auc(fprs, tprs)
    # Plot the ROC curve
    plt.figure()  
    plt.plot(fprs, tprs, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], 'k--', label='No Skill')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve for {filter.title}')
    plt.legend()
    plt.savefig(location)

"""
Write fpr,tpr, and parameters to csv
"""
def toCSV(location, rows):
    with open(location+"tsv", "w") as f:
        writer = csv.writer(f, delimiter='\t')
        header = ("fpr","tpr","param")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)

"""
Calculate the tprs and fprs of lengthFilter
"""
def rocLength(location):
    tprs = []
    fprs = []
    params = []
    for i in range(0,100,2):
        inputParams = f"{i},0" # 0 indicates absolute threshold, not percentage based 
        localFilter = lengthFilter(inputParams)

        truePos, allPos, falsePos, allNeg, outFilter = runFilter(localFilter, verbose = False)

        tpr = truePos/allPos
        fpr = falsePos/allNeg

        tprs.append(tpr)
        fprs.append(fpr)
        params.append(i)

        print(f"Completed param:{i}")
        print(f"{truePos}:{allPos}\t{falsePos}:{allNeg}")

        if i > 0 and (tpr == 1 or tprs[-1] - tprs[-2] < .015):
            break

    rows = sorted(zip(tprs, fprs, params), key=lambda x: x[1])
    toCSV(location, rows)

    tprs, fprs, params = zip(*rows)
    plotROC(tprs, fprs, location)

"""
Generate the list of regions that are effectively seperated by lengthFilter with given param
For both absolute length and percentage length.
"""
def rocLengthTest():
    # Absolute Difference
    for i in range(88, 500):
        inputParams =  f"{i},0"
        print("Absolute Dif", inputParams)
        localFilter = lengthFilter(inputParams)
        truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = False)
        print(f"{truePos}:{allPos}\t{falsePos}:{allNeg}")
        if truePos == allPos:
            break

    # Percentage of size based difference
    for i in np.arange(0, 100, 1/5):
        inputParams =  f"{i},i"
        print("Percentage Dif", inputParams)
        localFilter = lengthFilter(inputParams)
        truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = False)
        print(f"{truePos}:{allPos}\t{falsePos}:{allNeg}")
        if truePos == allPos:
            break

"""
Geneate the characteristic matricies of euclideanFilter
"""
def rocEuclidean():
    characteristicDirectory = "../../output/signatureMtxs_20bp/euclidean"
    # Test a range of bin sizes for assigning projected values
    for i in range(1,100):
        inputParams =  f"1000,21,{i},1"
        outputDirectory = os.path.join(characteristicDirectory,inputParams)

        if os.path.isdir(outputDirectory):
            print(f"Path already exists\t{outputDirectory}")
            continue

        os.makedirs(outputDirectory, exist_ok=True)
        localFilter = euclideanFilter(inputParams)
        truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = False, inplace = outputDirectory, loadEuclidean = True)

"""
Geneate the characteristic matricies of sketchFilter
"""
def rocSketch():
    characteristicDirectory = "../../output/signatureMtxs_20bp/sketch"

    inputParams =  f"1000,21,1"
    outputDirectory = os.path.join(characteristicDirectory,inputParams)

    if os.path.isdir(outputDirectory):
        print(f"Path already exists\t{outputDirectory}")
        return

    os.makedirs(outputDirectory, exist_ok=True)

    localFilter = sketchFilter(inputParams)
    truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = False, inplace = outputDirectory, loadSketch = True)

"""
Geneate the characteristic matricies of minHashFilter
"""
def rocMinHash():
    # TODO First one was to charMtxs/sketch on accident, fix when completed running
    characteristicDirectory = "../../output/signatureMtxs_20bp/minHash" 

    inputParams =  f"1000,21,1"
    outputDirectory = os.path.join(characteristicDirectory,inputParams)

    if os.path.isdir(outputDirectory):
        print(f"Path already exists\t{outputDirectory}")
        return

    os.makedirs(outputDirectory, exist_ok=True)

    localFilter = minHashFilter(inputParams)
    truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = False, inplace = outputDirectory, loadMinHash = True)

def main():
    filterFile = sys.argv[1]
    if filterFile == "lengthTest":
        rocLengthTest()
    elif "length" in filterFile:
        location = sys.argv[2]
        rocLength(location)
    elif "euclidean" in filterFile :
        rocEuclidean()
    elif "sketch" in filterFile :
        rocSketch()
    elif "minhash" in filterFile.lower():
        rocMinHash()
    else:
        print("Check filterFile name")

if __name__ == "__main__":
    main()