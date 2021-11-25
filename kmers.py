import sqlite3

class kmerParser:
    '''
    class used to parse DNA sequences
    ---
    Attributes:
        data (list): DNA sequence data stored as a list
    '''
    def __init__(self, inputData, inputType):
        '''
        Constructor
        ---
        Parameters:
            inputData (fastq file or string): DNA data input, currently the class supports fastq file parsing and DNA string parsing
            inputType (string): type of input data, as of now either string or fastq
        '''
        #if the input is a DNA sequence in form of a string 
        if inputType == 'string':
            self.data = [inputData]
            self.size = len(self.data[0]) #size of the string
        #if input is fastq file, use the class method fastqCleaner to process
        elif inputType == 'fastq':
            self.data = self.fastqCleaner(inputData)
            self.size = len(self.data[0]) #size of the first sequence, in the case of SP1.fastq, all sequences are the same length
        
    #getUniqueKmers(length: number): dictionary
    def getUniqueKmers(self, length):
        '''
        function which gets and counts all unique kmers of a spesified length
        ---
        Parameters:
            length(number): length of the target kmers
        ---
        Returns:
            (dictionary): dictionary with all unique kmers of spesified length as keys and counts as values
        '''
        kmersDict = {} #dict to store
        for sequence in self.data:
            for i in range(self.size - length + 1):
                target = sequence[i:i+length]
                if(target not in kmersDict): #not seen before
                    kmersDict[target] = 1
                else: #seen before +1 count
                    kmersDict[target] += 1
        return kmersDict

    #fastqCleaner(fastqFile: filepath): list
    def fastqCleaner(self, fastqFile):
        '''
        parses fastq file for DNA data
        ---
        Parameters:
            fastqFile(filepath): file to be parsed
        ---
        Returns:
            (list): parsed data
        '''
        lines = open(fastqFile, 'r').read().splitlines()
        data = [line for line in lines if line[0] in 'ACGT' and line.isalpha()]
        return data

    #toDatabase(kmerDictionary: dictionary): void
    def toDatabase(self, kmerDictionary, output):
        '''
        takes kmer dicitonary data and adds it to a sqlite3 database
        ---
        Parameters:
            kmerDictionary(dictionary): kmer data in the form of a dicitonary, can be aquired with getUniqueKmers method
            output(filepath): .db to output the data
        '''
        try: #create db and table if it does not exist
            con = sqlite3.connect(output)
            
            cur = con.cursor()

            cur.execute('''CREATE TABLE kmer (
                            sequence text,
                            count integer
                            )''')
        except sqlite3.OperationalError:
            print(output + 'already exitsts\n')
            pass
        for entry in kmerDictionary:
            cur.execute("INSERT INTO kmer VALUES ('{}','{}')".format(entry, kmerDictionary[entry]))
        con.commit()
        con.close()
        return


#sequence alignment algorithm using dynamic programing NOT USED
#turns out not to be needed for the match function since all sequences being compared are the same length
#regardless, it is useful for analyzing other mismatched sequences 
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


def seqMatchSameLength(s1,s2):
    '''
    simple sequence matching algorithm
    ---
    Parameters:
        s1(string): sequence 1
        s2(string): sequence 2
    ---
    Returns:
        (number): matching cost
    '''
    misMatchCount = 0
    for i in range(len(s1)):
        if(s1[i] != s2[i]):
            misMatchCount += 1
            if(misMatchCount > 2):
                return misMatchCount
    return misMatchCount

#match(kmer: string, seq: string): list[String]
def match(kseq, seq): 
    '''
    takes in target sequence and returns matching kmers in seq that differ by no more than 2 letters
    ---
    Parameters:
        kseq (string): a k-mer target
        seq (string): DNA sequence longer than k
    ---
    Returns:
        (list): list of k-mers in seq that match kseq with at most 2 letters different
    '''
    targetKmers = kmerParser(seq, 'string').getUniqueKmers(len(kseq))
    retSet = set()
    for seqence in targetKmers:
        if(seqMatchSameLength(seqence, kseq) <= 2):
            #just compare the strings directly, allignment is not important
            retSet.add(seqence)
    return retSet

def main():
    print('Welcome\n')
    print('Here is a sample of the match funciton in action: match(\'ACGT\', \'ACACACGT\') outputs:\n')
    print(match('ACGT', 'ACACACGT'))
    parser = kmerParser('SP1.fastq','fastq')
    kmerDict = parser.getUniqueKmers(21)
    parser.toDatabase(kmerDict, 'data.db')
    print('\nA data.db file should now show up in this program\'s directory containing all unique 21-mers and their counts as generated from SP1.fastq\n')

if __name__ == "__main__":
    main()

