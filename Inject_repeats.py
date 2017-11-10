from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




def injectRepeats(genomes, maxCNV):

    #Extract templates, they're in a list now!
    templates = list(SeqIO.parse("templates.fasta", "fasta"))
    templateNames = [t.id for t in templates]
    templates = [str(t.seq) for t in templates]
    print templateNames

    for genome in genomes:
        for copyNumber in range(1,maxCNV+1):

            currentGenome = str((SeqIO.read(genome + ".fasta", "fasta")).seq)

            for i, template in enumerate(templates):

                index = currentGenome.find(template)

                if index != -1:
                    currentGenome = currentGenome[:index - 500] + template * copyNumber + currentGenome[index + 500:]

                if index == -1:
                    print "Template " + templateNames[i] + " was not found."

            currentGenome = SeqRecord(Seq(currentGenome), id=str(genome) + "_" + "copyNumber:" + str(copyNumber))

            SeqIO.write(currentGenome, (str(genome) + "_" + "copyNumber:" + str(copyNumber) + ".fasta"), "fasta")





injectRepeats(["Beijing-391"], 10)

