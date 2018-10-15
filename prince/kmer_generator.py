import numpy as np
from Bio.Seq import Seq
# makes dict of strings, of size k

def kmer_generator(inputStrings,k):
    numberOfStrings = len(inputStrings)
    kmers = dict.fromkeys(np.arange(numberOfStrings))

    for j, string in enumerate(inputStrings):
        extendedString = string + string[:k-1]
        reverseString = Seq(extendedString).reverse_complement()
        kmers[j] = [extendedString[i:i + k] for i in range(0, len(extendedString) - k + 1)] +\
            [reverseString[i:i + k] for i in range(0, len(extendedString) - k + 1)]
    return kmers
