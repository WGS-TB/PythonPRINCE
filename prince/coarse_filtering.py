from Bio.Seq import Seq

def coarse_filtering(reads, k, template_kmers, quality_filter=20):
    recruitedReads = []
    template_set = set([x for v in template_kmers.values() for x in v])
    skipped_reads = 0
    
    for j,record in enumerate(reads):
        sequence  = record.seq
        if j == 0:
            read_length = len(sequence)
        #Generate k-splits
        reads = [sequence[n:n+25] for n in range(0,read_length-25+1,25)]
        
        quality = record.letter_annotations["phred_quality"]
        
        if j == 0:
            lengthOfReads = len(reads[0])
            leftoverPiece = lengthOfReads % k
            split_length = lengthOfReads - leftoverPiece
            k_splits_per_read = len(reads)
        
        for i,read in enumerate(reads):
            
            #Throw away low quality sections of the read and record the number of disposals.
            avg_qual = sum(quality[i*25:(i*25)+25])/25.0
            if avg_qual < quality_filter:
                skipped_reads +=1
                continue
            
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
    #Nucleotides seen takes into account number of reads, read length and skipped read sections
    nucleotides_seen = (j+1)*k_splits_per_read-skipped_reads
    return(nucleotides_seen, recruitedReads)
