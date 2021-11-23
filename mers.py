import difflib
import sqlite3
from difflib import SequenceMatcher, Differ

class kmerParser:
    #constructor()
    def __init__(self, input, inputType): #arg is sting or file for sequence (use for match)
        if inputType == 'string':
            self.data = [input]
            self.size = len(self.data[0])
        elif inputType == 'fastq':
            self.data = self.fastqCleaner(input)
            self.size = len(self.data[0])
        

    def getUniqueKmers(self, length): #fix to include counts of each kmer, object? hashmap or dict?
        kmersDict = {} #hashmap or dict?
        #fix for string
        for seq in self.data:
            for i in range(self.size - length + 1):
                target = seq[i:i+length]
                if(target not in kmersDict):
                    kmersDict[target] = 1
                else:
                    kmersDict[target] += 1
        return kmersDict

    def fastqCleaner(self, fastqFile):
        with open(fastqFile) as f:
            lines = f.read().splitlines()
        data = [line for line in lines if line[0] in 'ACGT' and line.isalpha()]
        return data

def toDB(dd):
    try:
        con = sqlite3.connect('data.db')
        cur = con.cursor()

        cur.execute('''CREATE TABLE kmer (
                        sequence text,
                        count integer
                        )''')
    except sqlite3.OperationalError:
        pass
    for entry in dd:
        cur.execute("INSERT INTO kmer VALUES ('{}','{}')".format(entry, dd[entry]))
    con.commit()
    return

test = kmerParser('SP1.fastq', "fastq")
dd = test.getUniqueKmers(21)
toDB(dd)

#not sure if i even need this :( since both are same length
#sequence alignment algorithm using dynamic programing 
def seqeuncMatchingCost(a,b,gapCost,matchCost):
    M = [] #2d memoization array
    #initialize array
    sub = []
    for i in range(len(a) + 1):
        for j in range(len(b) + 1):
            sub.append(0)
        M.append(sub)
        sub = []

    #top right to bottom left
    for i in range(len(a) + 1):
        for j in range(len(b) + 1):
            #base cases
            if(i == 0):
                M[i][j] = gapCost * j
            elif(j == 0):
                M[i][j] = gapCost * i
            else:
                #check if 2 chars are same to apply matchCost or not
                isMatch = a[i-1] == b[j-1]
                if(not isMatch):
                    M[i][j] = min(matchCost + M[i-1][j-1], gapCost + M[i][j-1], gapCost + M[i-1][j])
                else:
                    M[i][j] = min(M[i-1][j-1], gapCost + M[i][j-1], gapCost + M[i-1][j])
    return M[i][j]

#match(kmer: string, seq: string): list[String]
def match(kseq, seq): 
    #targetKmers = kmerParser(seq, 'string').getUniqueKmers(len(kseq))
    targetKmers = kmerParser(seq, 'fastq').getUniqueKmers(len(kseq))
    retSet = set()
    #brute force
    for seqence in targetKmers:
        if(seqeuncMatchingCost(seqence, kseq, len(seqence) ,1) <= 2):
            #just compare the strings directly, allignment is not important
            retSet.add(seqence)
    return retSet
    #use kmerParser functionality to get unique of lenth - change so that it can take in a string as arg (alt constructor)
    #matches with 2 or more...
    #1- get all unique kmers of size len(kseq) using above infastructure


#Each sequence is 31 nt long sp has 10 21 mers??each
#maybe create dict or hash map to store each and update cound during parsing

print(match('ACGTTTCACC','SP1.fastq'))




