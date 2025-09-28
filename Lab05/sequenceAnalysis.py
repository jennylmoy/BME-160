#!/usr/bin/env python3
# Name: Jennifer Moy (jelmoy)
# Group Members: Avirral Agarwal (avsagarw), Sharvari Bulbule (smbulbul), Kalpita Balu (kbalu)

class ProteinParam :
    '''
    Calculate the physical-chemical properties of the inputed protein sequence through methods.
    '''
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
        Instantiate parameters. Create an amino acid dictionary counter in constructor method and track amino acid characters.
        '''
        self.protein = protein.upper() # capitalize every character in input string
        self.aaDict = {'A': 0, 'C': 0, 'D': 0, 'E': 0,
                       'F': 0, 'G': 0, 'H': 0, 'I': 0,
                       'L': 0, 'K': 0, 'M': 0, 'N': 0,
                       'P': 0, 'Q': 0, 'R': 0, 'S': 0,
                       'T': 0, 'V': 0, 'Y': 0, 'W': 0
                      } # create a dictionary of each of the 20 amino acids 
        for key in self.protein: # keep track of how many amino acids each in input sequence
            if key in self.aaDict:
                self.aaDict[key] += 1
        

    def aaCount (self):
        '''
        In input amino acid sequence, count the number of valid amino acid characters.
        '''
        count = 0
        for i in self.protein: # iterate through input protein sequence and count the amount of valid amino acids
            if i in self.aaDict:
                count += 1
        return count 

    def pI (self):
        '''
        Find and return the pH that results in the lowest net charge from the _charge_ method. *Note: did a binary search for the extra credit.
        '''
        # Extra Credit: binary search 
        leftMost = 0.0
        rightMost = 14.0
        while rightMost - leftMost > 0.01:
            midMost = (leftMost + rightMost) / 2
            # if the net charge is positive, go to the upper bound of pH range
            if ProteinParam._charge_(self, midMost) > 0: 
                leftMost = midMost
            # if the net charge is negative, go to the lower bound of pH range
            elif ProteinParam._charge_(self, midMost) < 0:
                rightMost = midMost
        return (leftMost + rightMost)/ 2
            

    def aaComposition (self):
        '''
        Using instantiated dictionary from __init__ constructor method, which counted and tracked how many of each amino acid 
        characters are in the input sequence, return the updated dictionary.
        '''
        return self.aaDict

    
    def _charge_ (self, pH):
        '''
        Using the formula provided above and assigned names above, iterate through the charge 
        dictionaries to find the net charge.
        '''
        Nterminus = 0
        Cterminus = 0
        for key, value in self.aaDict.items(): # iterate through dictionary to compute net charge using provided formula
            if key in self.aa2chargePos:
                Nterminus += value * ((10**self.aa2chargePos[key])/(10**self.aa2chargePos[key] + 10**pH))
            if key in self.aa2chargeNeg:
                Cterminus += value * ((10**pH)/(10**self.aa2chargeNeg[key] + 10**pH))
        Nterminus += (10**self.aaNterm)/(10**self.aaNterm + 10**pH)
        Cterminus += (10**pH)/(10**self.aaCterm + 10**pH) 
        
        return Nterminus - Cterminus
             

    def molarExtinction (self, Cystine = True):
        '''
        Using the provided formula above, use the dictionary aa2abs280 for the molar extinction
        coefficients. Multiply the extinction coefficients by how many times they appear in the input
        sequence. Return the molar extinction value. *Note: did extra credit.
        '''
         # Extra Credit: Cystine default value for evaluation of molar extinction coefficient
        if Cystine == True: # if Cystine is True, compute molar extinction coefficient as normal
            molarExtinctionValue = 0
            for key, value in self.aa2abs280.items():
                molarExtinctionValue += self.aa2abs280[key] * self.aaDict[key]
            return molarExtinctionValue
        else: # if Cystine is False, exclude Cystine from the computation in molar extinction coefficient
            molarExtinctionValue = 0
            for key, value in self.aa2abs280.items():
                if key == 'C':
                    continue
                molarExtinctionValue += self.aa2abs280[key] * self.aaDict[key]
            return molarExtinctionValue

    
    def massExtinction (self, Cystine = True):
        '''
        Calculate the mass extinction coefficient by using the molar extinction method to call a value
        and the molecular weight method to call a value, then divide the values. *Note: did extra credit.
        '''
        # Extra Credit: Cystine default value for evaluation of mass extinction coefficient
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
        Find the total molecular weight of the input sequence using the provided formula above, molecular weight dictionary, and molecular weight of 
        water. Return the resulting total molecular weight.
        '''
        totalMolecularWeight = self.mwH2O
        for key,value in self.aaDict.items(): # iterate through dictionary to find molecular weight
            for aa, mwNum in self.aa2mw.items():
                if key == aa:
                    totalMolecularWeight += value * (mwNum - self.mwH2O)
        return totalMolecularWeight
   
import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

class NucParams:

    '''
    Extract amino acid, nucleotide, and codon composition from input sequence and return these dictionaries. 
    '''
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''
        Instantiate the parameters given in the method. Create empty dictionaries for tracking the compositions.
        '''
        self.inString = inString.upper() # capitlize the input string sequence
    
        self.validNuc = {'A', 'C', 'G', 'T', 'U', 'N'} # create a dictionary of valid nucleotides to filter out through the input sequence
        
        # create empty dictionaries for relevant compositions
        self.nucComp = {} 
        self.codonComp = {}
        self.aaComp = {}
        
        self.addSequence(inString) # call the addSequence() method
        
    def addSequence (self, inString):
        '''
        Count the number of nucleotides and number of codons in input RNA sequence. 
        '''
        # count the number of nucleotides
        for nucleotide in inString: # iterate through input sequence
            if nucleotide in self.validNuc: # filters only valid nucleotides from input sequence
                if nucleotide in self.nucComp: # adds nucleotides to new dictionary and tracks how many in sequence
                    self.nucComp[nucleotide] += 1 
                else:
                    self.nucComp[nucleotide] = 1


        # count the number of codon in input RNA sequence
        rnaNucString = inString.replace('T','U')
        for codon in range(0, len(rnaNucString), 3): # counts every three letters in input sequence
            codonSeq = rnaNucString[codon:codon + 3]
            if len(codonSeq) == 3 and codonSeq in self.rnaCodonTable:
                if codonSeq in self.codonComp: # adds codons to the new dictionary and tracks how many in the sequence
                    self.codonComp[codonSeq] += 1
                else:
                    self.codonComp[codonSeq] = 1

            aminoAcid = self.rnaCodonTable[codonSeq] # obtains the letter for the amino acid from the codon
            if aminoAcid in self.validNuc: # filters out only valid amino acids 
                if aminoAcid in self.aaComp: # add amino acids to the new dictionary and tracks how many in the sequence
                    self.aaComp[aminoAcid] += 1
                else:
                    self.aaComp[aminoAcid] = 1


    def aaComposition(self):
        '''
        Return the amino acid composition dictionary.
        '''
        return self.aaComp
    def nucComposition(self):
        '''
        Return the nucleotide composition dictionary.
        '''
        return self.nucComp
    def codonComposition(self):
        '''
        Return the codon composition dictionary.
        '''
        return self.codonComp
    def nucCount(self):
        '''
        Return the sum of the values in the nucleotide composition dictionary.
        '''
        return sum(self.nucComp.values())
    
class OrfFinder:
    '''
    Find the open reading frames (ORFs), including dangling start and stop fragments, in a given DNA sequence. Keep track of its start index position, stop index position, and frame. 
    Calculate the ORF length. 
    '''
    def __init__(self, minLength = 300, startCodon = {'ATG'}, stopCodons = {'TAG', 'TAA','TGA'}, longestOnly = True): 
        '''
        Initialize the parameter of the DNA sequence. Create dictionaries of the start and stop codons to be accessible within the class, and of the complimentary bases. 
        Initialize a name that reverses the forward DNA sequence. 
        '''
        self.minLength = minLength
        self.startCodon = startCodon
        self.stopCodons = stopCodons
        self.longestOnly = longestOnly
        self.complimentaryBases = {'T': 'A', 'A': 'T', 'G': 'C', 'C': 'G'} # complimentary bases to convert forward to reverse strand

        # initialized finalized lists without duplicates 
        self.finalizedFinalizedORFsForwardOutput = [] 
        self.finalizedFinalizedORFsReverseOutput = []

    def findForwardGenes(self, sequence):
        '''
        Find every ORF in the DNA sequence in the forward frames (+1, +2, +3). Return the information of ORF length, start index position, stop index position, and frame.
        '''
        ORFsForwardOutput = [] # create a list to append tuples of ORF information 
    
        for frame in range(3): # frame shift (+1 = 0), (+2 = 1), (+3 = 2)
        
            startIndex = None # initialize the startIndex as None
            startCodonFound = False # initialize the startCodonFound as False (as a flag in program)
            stopCodonFound = False # initialize stopCodonFound as False (flag)
            startCodonList = [] # keeps track of start codon index positions

            for codonIndex in range(frame, len(sequence)- 2, 3): # counting the codons in the forward strand, for each frame shift
                codon = sequence[codonIndex:codonIndex + 3] # find the nucleotides of the codon
                
                # Case 1: Beginning a new ORF and found a start codon
                if codon in self.startCodon:
                        startIndex = codonIndex + 1 # assign the start index as the codon index plus 1 (due to python's indexing differing from DNA indexing)
                        startCodonFound = True # switch the flagging name to True
                        startCodonList.append(startIndex) # append start codon position to list 

                # Case 2: Already begun ORF counting after finding a start codon, and found a stop codon
                if startIndex is not None and codon in self.stopCodons: 
                        stopIndex = codonIndex + 3 # assign the stop codon index as the codon index plus 3
                        ORFLength = stopIndex - startIndex + 1 # calculate the ORF length
                        ORFsForwardOutput.append((f'+{frame+1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                        startIndex = None # reset the start index position
                        stopCodonFound = True # found a stop codon

                        ORFsForwardOutput.append((f'+{frame+1}', 1, len(sequence), len(sequence)))

                        if self.longestOnly is False: # when longest gene only is False
                            for start in startCodonList: # iterate through startCodonList for each ORF
                                ORFsForwardOutput.append((f'+{frame+1}', start, stopIndex, stopIndex - start + 1))
                        else: # when longest gene only is True
                            ORFsForwardOutput.append((f'+{frame+1}', startCodonList[0], stopIndex, stopIndex - startCodonList[0] + 1))

                # Case 3: When no start codon was found, but a stop codon was found (dangling stop fragment)
                if not startCodonFound and codon in self.stopCodons: 
                    startIndex = 1 # start codon position is at the beginning of DNA sequence
                    stopIndex = codonIndex + 3 
                    ORFLength = stopIndex - startIndex + 1 # +1 is inclusive of the last position of the stop codon
                    if self.longestOnly is False: # when longest gene only is False
                        ORFsForwardOutput.append((f'+{frame+1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                    startIndex = None # reset the start index position
                    stopCodonFound = True # found a stop codon

             # Case 4: After beginning a ORF, no stop codons have been found and reached end of DNA sequence (dangling start fragment)
            if startIndex is not None:
                stopIndex = len(sequence) 
                ORFLength = stopIndex - startIndex + 1 # +1 is inclusive of the last position of the stop codon
                ORFsForwardOutput.append((f'+{frame + 1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                if self.longestOnly is False: # when longest gene only is False
                    for start in startCodonList: # iterate through startCodonList for each ORF
                        ORFsForwardOutput.append((f'+{frame+1}', start, stopIndex, ORFLength))
                else: # when longest gene only is True
                    ORFsForwardOutput.append((f'+{frame+1}', startCodonList[0], stopIndex, stopIndex - startCodonList[0] + 1))

            # Case 5: If no start codon was found and no stop codon was found
            if startIndex is None and stopCodonFound is False:
                ORFsForwardOutput.append((f'+{frame+1}', 1, len(sequence), len(sequence))) # this appends the ORF that spans the entire sequence

        finalizedORFsForwardOutput = list(set(ORFsForwardOutput)) # remove any duplicates 
        for ORFList in finalizedORFsForwardOutput:
            if ORFList[3] > self.minLength: # iterating through list to remove any lengths smaller than the minimum length
                self.finalizedFinalizedORFsForwardOutput.append(ORFList) # add lengths larger than minimum length to a new list 
            
        return self.finalizedFinalizedORFsForwardOutput

    def findReverseGenes(self, sequence):
        '''
        Find every ORF in the DNA sequence in the reverse frames (-1, -2, -3). Return the information of ORF length, start index position, stop index position, and frame.
        '''
        reverseStrand = ''.join([self.complimentaryBases[nucleotide] for nucleotide in sequence[::-1]]) # reverses the forward DNA sequence and switches the complimentary bases
        ORFsReverseOutput = [] # create a list to append tuples of ORF information (for reverse)
        for frame in range(3): # frame shift (0 = -1), (1 = -2), (2 = -3)

            startIndex = None # initialize the startIndex as None
            startCodonFound = False # initialize the startCodonFound as False (as a flag in program)
            stopCodonFound = False # initialize stopCodonFound as False (flag)
            startCodonList = [] # keeps track of start codon index positions

            for codonIndex in range(frame, len(reverseStrand)- 2, 3): # counting the codons in the reverse strand, for each frame shift
                codon = reverseStrand[codonIndex:codonIndex + 3] # find the nucleotides of the codon
                
                # Case 1: Beginning a new ORF and found a start codon
                if codon in self.startCodon:
                        startIndex = codonIndex + 1 # convert to forward positionings, +1 due to python index starting at 0, not 1
                        startCodonFound = True # switch the flagging name to True
                        startCodonList.append(startIndex) # append start codon position to list
                        
                # Case 2: Already begun ORF counting after finding a start codon, and found a stop codon 
                if startIndex is not None and codon in self.stopCodons: 
                        stopIndex = codonIndex + 3 # substract the length of strand minus 'index' to get the forward positionings
                        ORFLength = stopIndex - startIndex + 1 # +1 is inclusive of the last position of the stop codon
                        ORFsReverseOutput.append((f'-{frame+1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                        startIndex = None # reset the start index position
                        stopCodonFound = True # found a stop codon

                        ORFsReverseOutput.append((f'-{frame+1}', 1, len(sequence), len(sequence)))
                        
                        if self.longestOnly is False: # when longest gene only is False
                            for start in startCodonList: # iterate through startCodonList for each ORF
                                ORFsReverseOutput.append((f'+{frame+1}', start, stopIndex, stopIndex - start + 1))
                        else: # when longest gene only is True
                            ORFsReverseOutput.append((f'+{frame+1}', startCodonList[0], stopIndex, stopIndex - startCodonList[0] + 1))

                # Case 3: When no start codon was found, but a stop codon was found (dangling stop fragment)
                if not startCodonFound and codon in self.stopCodons: 
                    startIndex = 1 # start codon position is at the beginning of DNA sequence
                    stopIndex = codonIndex + 3  # convert to forward positionings
                    ORFLength = stopIndex - startIndex + 1 # +1 is inclusive of the last position of the stop codon
                    if self.longestOnly is False: # when longest gene only is False
                        ORFsReverseOutput.append((f'-{frame+1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                    startIndex = None # reset the start index position
                    stopCodonFound = True # found a stop codon

            # Case 4: after beginning a ORF, no stop codons have been found and reached end of DNA sequence (dangling start fragment)
            if startIndex is not None: 
                stopIndex = len(reverseStrand) 
                ORFLength = stopIndex - startIndex + 1 # +1 is inclusive of the last position of the stop codon
                ORFsReverseOutput.append((f'-{frame + 1}', startIndex, stopIndex, ORFLength)) # append the ORF information into a list

                if self.longestOnly is False: # when longest gene only is False
                    for start in startCodonList: # iterate through startCodonList for each ORF
                        ORFsReverseOutput.append((f'+{frame+1}', start, stopIndex, ORFLength))
                else: # when longest gene only is True
                    ORFsReverseOutput.append((f'+{frame+1}', startCodonList[0], stopIndex, stopIndex - startCodonList[0] + 1))

            # Case 5: If no start codon was found and no stop codon was found
            if startIndex is None and stopCodonFound is False:
                ORFsReverseOutput.append((f'-{frame+1}', 1, len(sequence), len(sequence))) # this appends the ORF that spans the entire sequence

        finalizedORFsReverseOutput = list(set(ORFsReverseOutput)) # remove any duplicates 
        for ORFList in finalizedORFsReverseOutput:
            if ORFList[3] > self.minLength: # iterating through list to remove any lengths smaller than the minimum length
                self.finalizedFinalizedORFsReverseOutput.append(ORFList) # add lengths larger than minimum length to a new list 
                
        return self.finalizedFinalizedORFsReverseOutput

    def sortedORFs(self, sequence):
        '''
        Sorts the lists of forward and reverse ORFs into order by longest to shortest length. If the lengths are the same for different ORFs,
        it will next sort by startIndex.
        '''
        forwardORFs = self.findForwardGenes(sequence) # call the findForwardGenes method and initialize it
        reverseORFs = self.findReverseGenes(sequence) # call the reverseForwardGenes method and initialize it
        
        totalORFs = forwardORFs + reverseORFs # append the lists together
        sortedORFs = sorted(totalORFs, key=lambda x: (-x[3], x[1], x[0])) # sort the list 

        return sortedORFs 


