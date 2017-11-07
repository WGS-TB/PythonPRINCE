import itertools
import numpy as np


#Find the incomplete repeat overhang
def find_length_of_incomplete_repeat(templates, k, reads):
    incompleteRepeats = [0] * len(templates)
    for temp, template in enumerate(templates):
        t = template + template[0:k-1]
        templateKmersOriginal = [t[i:i + k] for i in range(0, len(t) - k + 1)]

        seen = set()
        seen_add = seen.add
        templateKmers = [x for x in templateKmersOriginal if not (x in seen or seen_add(x))]

        numberOfTemplateKmers = len(templateKmers)
        fortyPercentMatch = numberOfTemplateKmers * 0.4
        eightyPercentMatch = numberOfTemplateKmers * 0.8

        endsOfReads = []

        for read in reads:
            print "read", read
            match = [0] * numberOfTemplateKmers
            readKmers = [read[i:i + k] for i in range(0, len(read) - k + 1)]
            # print readKmers

            for readkmer in readKmers:
               # print templateKmers
               tempMatch = [readkmer == km for km in templateKmers]
               match = [m + t for m, t in zip(match, tempMatch)]
               print match

            onesMatch = [1 if x > 0 else 0 for x in match]
            longestOnes= max(sum(1 for _ in l) for n, l in itertools.groupby(onesMatch))

            if longestOnes > fortyPercentMatch and longestOnes < eightyPercentMatch:
                endsOfReads.append(longestOnes)

        if len(endsOfReads) > 0:
            incompleteRepeats[temp] = np.bincount(endsOfReads).argmax()
        else:
            incompleteRepeats[temp] = 0


    print incompleteRepeats
    return incompleteRepeats

#
# temps = ["1234"]
# reads = ["4123", "2341", "6712", "1234", "2345", "3456"]
#
# find_length_of_incomplete_repeat(temps, 3, reads)
#
#
#




