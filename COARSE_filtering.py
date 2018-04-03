from Bio.Seq import Seq

def coarse_filtering(reads, k, template_kmers):
    numberOfReads = len(reads) * 2
    lengthOfReads = len(reads[0])
    leftoverPiece = lengthOfReads % k
    halfOfReads = numberOfReads/2
    recruitedReads = []
    template_set = set([x for v in template_kmers.values() for x in v])


    for j, read in enumerate(reads):
        # Forward
        readSeq = str(read.seq)
        remainingRead = readSeq[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in template_set for r in listOfSplits]

        if True in recruited:
            recruitedReads.append(readSeq)

        # Reverse
        reverseRead = str((read.reverse_complement()).seq)
        remainingRead = reverseRead[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in template_set for r in listOfSplits]

        if True in recruited:
            recruitedReads.append(reverseRead)

    return(recruitedReads)

