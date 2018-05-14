import numpy as np
# makes dict of strings, of size k

def kmerGenerator(inputStrings,k):
    numberOfStrings = len(inputStrings)
    kmers = dict.fromkeys(np.arange(numberOfStrings))

    for j, string in enumerate(inputStrings):
        extendedString = string + string[:k-1]
        kmers[j] = [extendedString[i:i + k] for i in range(0, len(extendedString) - k + 1)]
    return kmers
