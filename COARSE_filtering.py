from Bio.Seq import Seq

def coarse_filtering(reads1, reads2, k, template_kmers):
    numberOfReads = len(reads1) * 2
    lengthOfReads = len(reads1[0])
    leftoverPiece = lengthOfReads % k
    halfOfReads = numberOfReads/2

    recruitedForwardReads = [""] * numberOfReads
    recruitedReverseReads = [""] * numberOfReads



    for j, read in enumerate(reads1):
        # Forward
        readSeq = str(read.seq)
        remainingRead = readSeq[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in [x for v in template_kmers.values() for x in v]for r in listOfSplits]

        if True in recruited:
            recruitedForwardReads[j] = readSeq



        # Reverse
        reverseRead = str((read.reverse_complement()).seq)
        remainingRead = reverseRead[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in [x for v in template_kmers.values() for x in v] for r in listOfSplits]

        if True in recruited:
            recruitedReverseReads[j] = reverseRead


    for p, read in enumerate(reads2):
        # Forward
        readSeq = str(read.seq)
        remainingRead = readSeq[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in [x for v in template_kmers.values() for x in v] for r in listOfSplits]

        if True in recruited:
            recruitedForwardReads[halfOfReads + p] = readSeq

        # Reverse
        reverseRead = str((read.reverse_complement()).seq)
        remainingRead = reverseRead[:(lengthOfReads - leftoverPiece)]
        listOfSplits = [remainingRead[i:i + k] for i in range(0, len(remainingRead), k)]
        recruited = [r in [x for v in template_kmers.values() for x in v] for r in listOfSplits]

        if True in recruited:
            recruitedReverseReads[halfOfReads + p] = reverseRead


    recruitedForwardReads = filter(None, recruitedForwardReads)
    recruitedReverseReads = filter(None, recruitedReverseReads)


    return recruitedForwardReads + recruitedReverseReads


