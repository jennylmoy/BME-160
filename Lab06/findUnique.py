#!/usr/bin/env python3
# Name: Jennifer Moy (jelmoy)
# Group Members: Avirral Agarwal, Kalpita Balu, Sharvari Bulbule

import sys
class FastAreader :
    
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
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

class tRNAFinder:
    '''
    Find the essential tRNA substring sets from 22 tRNA sequences by creating power sets of each tRNA sequence and finding the unique substrings through all the sequences. 
    Ouput the header, sequence, and respective essential substrings for each tRNA.
    '''

    def __init__(self):
        '''
        Initialize dictionaries to be used within the methods.
        '''        
        self.tRNAHeaderSequenceDict = {} # dictionary {header: sequence}
        self.tRNAPowerSetDict = {} # dictionary {header: all substrings (power set) of the specific sequence}
        self.tRNAUniqueDict = {} # dictionary {header: all unique substrings of the specific sequence}
        self.tRNAEssentialDict = {} # dictionary {header: all essential substrings of the specific sequence}
        
    def powerSet(self, dictionary): # dictionary = self.tRNAHeaderSequenceDict {header: sequence}
        '''
        Find the power set for each tRNA sequence by itearting through a dictionary of headers and sequences. Append each header and its power set to a new dictionary. 
        '''
        for header, sequence in dictionary.items(): # iterating through the header and sequence of self.tRNAHeaderSequenceDict dictionary
            tRNAsubstringPowerSet = set() # temporary set for all the substrings of a specific header sequence
            for i in range(len(sequence)): # outer loop
                for j in range(i + 1, len(sequence) + 1): # inner loop
                    substring = sequence[i:j] # creates substring of tRNA set
                    tRNAsubstringPowerSet.add(substring) # add substring to list
              
                self.tRNAPowerSetDict[header] = tRNAsubstringPowerSet # append the header and all substrings of the sequence into  a new dictionary
        return self.tRNAPowerSetDict

    def uniqueSet(self, dictionary): # dictionary = self.tRNAPowerSetDict {header: all substrings (power set) of the specific sequence}
        '''
        Find the unique substrings of a power set by comparing to other tRNA power sets. Append each header and its unique substring set to a new dictionary.
        '''
        tRNAunique = set() # set of the unique substrings of a specific sequence

        for header, powerset in dictionary.items(): # iterating through the header and powerset of self.tRNAPowerSetDict
            othertRNASubstrings = set() # set of substrings of not the specific sequence (other tRNAs)
            for comparingHeader, comparingPowerset in dictionary.items(): # iterating through the comparing header and powerset of same dictionary (self.tRNAPowerSetDict)
                if comparingHeader != header: # if the headers don't match
                    othertRNASubstrings.update(comparingPowerset) # add the substrings into the othertRNASubstrings set 
            tRNAunique = powerset - othertRNASubstrings # substract
            self.tRNAUniqueDict[header] = tRNAunique # append the unique substrings set into a dictionary
        return self.tRNAUniqueDict

    def essentialSet(self):
        '''
        Sort dictionaries by length or alphabetically. Find the essential tRNA substrings of a certain tRNA sequence, 
        contains the substrings with mimimal character length and no repetition of starting nucleotides. Append each header and its essential substring set to a new dictionary. 
        '''

        tempSortedtRNAunique = {} # placeholder dictionary {header: sorted unique sets by length}

        for header, uniqueSets in self.tRNAUniqueDict.items(): # iterating through the header and unique substring sets of self.tRNAUniqueDict
            sortedElements = sorted(uniqueSets, key=len) # sort the unique substring sets by their length
            tempSortedtRNAunique[header] = sortedElements # append the header and sorted (by length) unique substring sets into a new dictionary

        sortedtRNAUniqueDict = dict(sorted(tempSortedtRNAunique.items())) # sorts headers of placeholder dictionary (tempSortedtRNAunique) alphabetically, {sorted headers: sorted unique sets}

        for header, unique in sortedtRNAUniqueDict.items(): # iterating through the header and unique substring sets in dictionary (sortedtRNAUniqueDict)
            essentialSubstrings = set() # create a set of essential substrings 
            nonEssentialSubstrings = set() # create a set of non-essential substrings ("trash")

            for index, substring in enumerate(unique): # iterating through the index and substring of the unique substring set 
                if substring in nonEssentialSubstrings: # if substring in non-essential
                    nonEssentialSubstrings.add(substring) # "trash" substring
                else:
                    essentialSubstrings.add(substring) # if substring not in non-essential, keep it
                for comparingSubstring in unique[index+1:]: # comparing the substrings with the rest in the unique set
                    if comparingSubstring.startswith(substring): # if the comparing substrings start with another substring
                        nonEssentialSubstrings.add(comparingSubstring) # "trash" it

            self.tRNAEssentialDict[header] = essentialSubstrings # append all the essential substrings of a specific sequence into a new dictionary
            essentialSubstrings = set() # empty the set to start on a new sequence

        return self.tRNAEssentialDict
    
    def tRNAOutput(self):
        '''
        Sort respective dictionaries alphabetically and numerically. Output each header and sequence with its respective substrings and position represented by periods. 
        '''
        # self.tRNAEssentialDict = {header: all essential substrings of the specific sequence}; essential set = already sorted alphabetically but not by position of substrings in sequence
        
        tRNADict = {} # dictionary {substring: position}
        tRNADictWithHeader = {} # dictionary {header: its specific substrings tRNADict}
        
        sortedHeaderSequenceDict = dict(sorted(self.tRNAHeaderSequenceDict.items()))  # dictionary {alphabetically sorted header: sequence}

        for header, substringSet in self.tRNAEssentialDict.items(): # iterating through the header and substring in self.tRNAEssentialDict 
            for comparingHeader, sequence in sortedHeaderSequenceDict.items(): # iterating through a comparing header and sequence in sortedHeaderSequenceDict
                if header == comparingHeader: # if the headers mathc
                    for substring in substringSet: # iterating through each substring in a specific set 
                        position = sequence.find(substring) # find the position of the substring in the sequence
                        tRNADict[substring] = position # append the substring and its position to a new dictionary
            tRNADictWithHeader[header] = tRNADict
        
        

        tRNAsubstringDict = {} # dictionary {substring: position}
        sortedtRNAsubstringDict = {} # dictionary {substring: position sorted numerically}

        for header, dictionary in tRNADictWithHeader.items(): # iterating through header and dictionary of {substring: position} in tRNADictWithHeader
            for tRNAsubstring, position in dictionary.items(): # iterating through the substring and its position in that dictionary
                tRNAsubstringDict[tRNAsubstring] = position # append substring and position into a new dictionary
                sortedtRNAsubstringDict[header] = dict(sorted(tRNAsubstringDict.items(), key=lambda x: x[1])) # sort that dictionary by position index

    
        for header, sequence in sortedHeaderSequenceDict.items(): # iterating through the header and sequence in self.tRNAHeaderSequenceDict
            print(header) # print the corresponding header
            print(sequence) # print the corresponding sequence
            for comparingheader, dictionary in sortedtRNAsubstringDict.items(): # iterating through the header and dictionary of sortedtRNAsubstringDict
                if header == comparingheader: 
                    for tRNAsubstring, position in dictionary.items(): # iterating through substring and sorted position of the dictionary
                        if tRNAsubstring in sequence:
                            period = '.' * position # multiply the period by the position value
                            print(f'{period}{tRNAsubstring}') # output the substrings with period values 

def main(inCL=None):
    ''' 
    Creates a dictionary of headers and sequences utilizing the FastAreader class. Removes unreadable characters (',' ; '_' ; '.', '-') from the header and sequence.
    Call the tRNAFinder class to output all the headers, sequences, and respective substrings. 
    '''
    reader = FastAreader() 
    finder = tRNAFinder()
  
    # create the self.tRNAHeaderSequenceDict {header: sequence}
    for header, sequence in reader.readFasta(): 
        cleanedHeader = header.replace(' ', '') # clean header
        cleanedSequence = sequence.replace('.','').replace('_','').replace(',','').replace('-','') # clean sequence
        finder.tRNAHeaderSequenceDict[cleanedHeader] = cleanedSequence # add cleaned header and cleaned sequence to a dictionary

    # call the methods
    finder.powerSet(finder.tRNAHeaderSequenceDict)
    finder.uniqueSet(finder.tRNAPowerSetDict)
    finder.essentialSet()
    finder.tRNAOutput()
        
if __name__ == "__main__":
    main()  

