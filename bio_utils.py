from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_rna

class Orf:    
    def __init__(self, start, stop, frame, value, length, name):
        self.start = start
        self.stop = stop
        self.value = value
        self.frame = frame
        self.length = length
        self.name = name
    
    def __str__(self):
        return '{} <Start: {}, Stop: {}, Length: {}, Frame: {}>'.format(self.name, self.start, self.stop, self.length, self.frame)
    
class Comparison:
    def __init__(self, firstOrf, secondOrf, alignment):
        self.firstOrf = firstOrf
        self.secondOrf = secondOrf
        self.alignment = alignment
    def __str__(self):
        return 'Comparison <Score: {}>'.format(self.alignment.score)
    
def toRNAs(orfs):
    rnas = []
    for orf in orfs:
        rnas.append(toRNA(orf.value))
    return rnas

def toRNA(dna):
    conv = {"T" : "A", "A" : "U", "C" : "G", "G" : "C"}
    rna = ""
    for i in dna:
        rna += str(conv[i])
    return rna

def toProteins(rnas):
    proteins = []
    for rna in rnas:
        r = Seq(rna, generic_rna)
        proteins.append(r.translate())
    return proteins

def writeToFile(filepath, data):
    with open(filepath, 'w') as the_file:
        for i in data:
            the_file.write("{}\n\n".format(i))