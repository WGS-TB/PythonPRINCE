from Bio.Seq import Seq

def coarse_filtering(reads, k, template_kmers):
    recruitedReads = []
    template_set = set([x for v in template_kmers.values() for x in v])

    for j,record in enumerate(reads):
	sequence  = record.seq
	if j ==0 :
	    read_length = len(sequence)

	#Generate k-splits
	reads = [sequence[n:n+25] for n in range(0,read_length-25+1,25)]
	
	if j == 0:
	    lengthOfReads = len(reads[0])
	    leftoverPiece = lengthOfReads % k
	    split_length = lengthOfReads - leftoverPiece
	    k_splits_per_read = len(reads)

	for read in reads:	
            # Forward
            readSeq = str(read)
            listOfSplits = [readSeq[i:i + k] for i in range(0, split_length, k)]

            recruited = [r in template_set for r in listOfSplits]

            if True in recruited:
                 recruitedReads.append(readSeq)
            
	    # Reverse
            reverseRead = str((read.reverse_complement()))
            listOfSplits = [reverseRead[i:i + k] for i in range(0, split_length, k)]

	    recruited = [r in template_set for r in listOfSplits]

	    if True in recruited:
                 recruitedReads.append(reverseRead)

    return((j+1)*k_splits_per_read, recruitedReads)

