import numpy as np

def genRandomPlanes():
    m = 24000
    with open("../../output/randomStorage/randPlanes", "w") as f:
        for _ in range(20000):
            orthonormal_line = np.random.choice(2, size = m)
            orthonormal_line[orthonormal_line == 0] = -1
            np.savetxt(f, orthonormal_line, fmt="%i", newline="\t")
            f.write('\n')

def genRandomLines():
    m = 24000
    with open("../../output/randomStorage/randLines", "w") as f:
        with open("../../output/randomStorage/randOffsets", "w") as g:
            for _ in range(20000):
                s = np.random.normal(0, 1, m)
                direction = s / np.sqrt(sum(s*s))
                np.savetxt(f, direction, newline="\t")
                f.write('\n')

                offset = np.random.uniform(1,5,m)
                np.savetxt(g, offset, newline="\t")
                g.write('\n')

def main():
    genRandomPlanes()
    genRandomLines()

if __name__ == "__main__":
    main() #Don't want to accidentally overwrite data
    pass