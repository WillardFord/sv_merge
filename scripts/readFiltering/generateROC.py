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
Calculate the tprs and fprs of lengthFilter
"""
def rocEuclidean():
    characteristicDirectory = "../../output/charMtxs/euclidean"

    inputParams =  f"1000,10,5,1"
    outputDirectory = os.path.join(characteristicDirectory,inputParams)
    os.makedirs(outputDirectory, exist_ok=True)

    localFilter = euclideanFilter(inputParams)
    truePos, allPos, falsePos, allNeg = runFilter(localFilter, verbose = True, inplace = outputDirectory)

def main():
    filterFile = sys.argv[1]

    if "length" in filterFile:
        location = sys.argv[2]
        rocLength(location)
    elif "euclidean" in filterFile :
        rocEuclidean()
    else:
        print("Check filterFile name")

if __name__ == "__main__":
    main()