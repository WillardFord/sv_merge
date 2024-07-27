"""
Runnable File containing a lengthFilter class:

python lengthFilter.py --param 10,1

Second param [0|1] indicates [absolute threshold|percent length threshold]
"""
from test_filters import Filter, runFilter, readInputs

class lengthFilter(Filter):
    '''
    seqs: length n list of sequences
    classes: n by n matrix indicating if two sequences are sufficiently similar
    threshold: maximum size difference before connection becomes impossible
    '''
    def __init__(self, seqs, labels, threshold):
        self.threshold = threshold
        super().__init__(seqs, labels)

    def __init__(self, param):
        params = param.split(",")
        self.threshold = float(params[0])
        self.percent = bool(params[1])
        self.title = f"lengthFilter_threshold_{self.threshold}"

        if not self.percent:
            self.threshold = int(self.threshold)

    '''
    Preprocess all reads, required for some filters
    '''
    def processReads(self, test = False):
        # Set max length
        if self.percent:
            self.maxLength = 0
            for i in range(self.n):
                if len(self.seqs[i]) > self.maxLength:
                    self.maxLength = len(self.seqs[i])

        for i in range(self.n-1):
            for j in range(i, self.n):
                if i == j:
                    self.adjacencyMatrix[i,j] = 1
                else:
                    self.adjacencyMatrix[i,j] = self.run_connect(i,j)
        return self.adjacencyMatrix

    """
    Conect Edges
    """
    def run_connect(self, i, j):
        distance = abs(len(self.seqs[i]) - len(self.seqs[j]))
        if self.percent: 
            return  distance <= self.maxLength * self.threshold
        return distance <= self.threshold

    '''
    Legacy function from first implementation. Can probably remove this from every single filter but will do later.
    '''
    def connect(self, i, j):
        return self.adjacencyMatrix[i,j]
        distance = abs(len(self.seqs[i]) - len(self.seqs[j]))
        if self.percent: 
            return  distance <= self.maxLength * self.threshold
        return distance <= self.threshold

def main():
    saveFig, param, test = readInputs()
    filter:lengthFilter = lengthFilter(param)
    runFilter(filter, saveFig, test, verbose = False)

if __name__ == "__main__":
    main()