import os
from multiprocessing import cpu_count

import pandas as pd
from joblib import Parallel, delayed
import numpy as np

def save_lengths():
    sample_dir = "../../output/HG002_20bp"
    save_dir = "../../output/analysis"
    for chrom in os.listdir(sample_dir):
        chrom_dir = os.path.join(sample_dir, chrom)
        input = os.listdir(chrom_dir)


        num_cores = cpu_count()
        outputs = Parallel(n_jobs=num_cores, verbose=1)(delayed(get_lengths)(i) for i in input)

        df = pd.DataFrame(data=outputs, columns=)



        print(f"Saved {chrom}")



def get_lengths(fastq_file):
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
    
    seqs, labels = readFastq(fastq_file)
    seqs, labels = zip(*sorted(zip(seqs, labels), key = lambda x : x[1]))
    
    n = len(seqs)
    lengths = np.zeros((n), dtype = int)
    for i in range(n):
        lengths[i] = len(seqs)
    np.save



def main():
    pass


if __name__ == "__main__":
    main()