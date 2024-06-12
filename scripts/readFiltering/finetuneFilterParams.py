'''
Given filter file name, finetune the parameter over all regions and all samples

python finetuneFilterParams.py hashFilter.py 1 20

TODO Future steps should consider that different regions have different properties --> Only optimize over certain labeled regions
    i.e. highly repetitive regions, [types of structural variants] ask Fabio and Ryan for a competant list.
'''
import subprocess
import sys
from os.path import join
import matplotlib.pyplot as plt
import numpy as np

def main():
    file, start_param, end_param = loadArguments()

    filterType = file.split("Filter")[0]

    tprs, fprs, params = optimizeParam(file, start_param, end_param)

    # Sort by fpr
    tprs, fprs, params = zip(*sorted(zip(tprs, fprs, params), key = lambda x: x[1]))

    plotTPRbyFPR(tprs, fprs, params, filterType)


def loadArguments():
    file = sys.argv[1]
    start_param = sys.argv[2]
    end_param = sys.argv[3]
    return file, int(start_param), int(end_param)

def optimizeParam(file, start_param, end_param, n = 20):
    tprs = []
    fprs = []
    params = []
    for p in range(start_param, end_param):
        tpr, fpr = getTPR_FPR(file, p)
        tprs.append(tpr)
        fprs.append(fpr)
        params.append(p)
    return tprs, fprs, params

def getTPR_FPR(file, param):
    workingDir = "/Users/wford/Documents/sv_merge/scripts/readFiltering"
    script = join(workingDir, file)
    command = ["python", script, "--param", str(param), "--plot", "False"]
    tpr_fpr = subprocess.run(command, 
                             capture_output = True, text = True)
    tpr, fpr = tpr_fpr.stdout.strip().split("\t")
    return float(tpr), float(fpr)

def plotTPRbyFPR(tprs, fprs, params, filterType, outDir = "../../output/plots"):
    plt.figure(figsize=(8, 6))

    # Create a color gradient based on the params
    colors = np.arange(len(params))
    plt.scatter(fprs, tprs, c=colors, cmap='viridis', marker='o', label='Data Points', edgecolors='black', linewidths=0.5)

    # Add labels for each point
    #for i, txt in enumerate(params):
    #    plt.annotate(txt, (fprs[i], tprs[i]), textcoords="offset points", xytext=(0,10), ha='center')

    x_ticks = np.linspace(min(fprs), max(fprs), num=5)
    plt.xticks(x_ticks, [round(x, 2) for x in x_ticks])

    # Customize y-axis ticks
    y_ticks = np.linspace(min(tprs), max(tprs), num=5)
    plt.yticks(y_ticks, [round(y, 2) for y in y_ticks])

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'TPR vs FPR of {filterType}')
    plt.colorbar(label='Parameter')
    plt.legend()

    saveLocation = join(outDir, f"paramOptimization_{filterType}.jpg")
    plt.savefig(saveLocation) 


if __name__ == "__main__":
    main()