from sequenceAnalysis import FastAreader

class CommandLine():

    """
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was

    set as an option using add_argument, then myCommandLine.args.requiredbool will

    name that option.
    """
    def __init__(self, inOpts=None):
        """
        Implement a parser to interpret the command line argv string using argparse.
            --> Argparse adjusted so that default is minGene = 100
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), action='store',
                                 help='minimum Gene length', default = 100)
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class OrfFinder:
    '''This class will find ORFs in the given sequence as well as find the reverse compliment strand of DNA and also find orfs inthat strand.'''
    stopCodons = ['TGA', 'TAG', 'TAA'] # create a list of stop codons 
    startCodons = ['ATG'] # create a list of start codons 
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', } # replace the nucleotide bases to form the compliment 

    def __init__(self, seq):
        self.seq = seq # initilize the variable seq 
        self.orfs = [] # create an empty list for the orfs to be put into 

    def findOrfs(self): 
        '''In this method, orfs will be found in the DNA sequence and stored in a list.
        This will be done by iterating through codons (sets of three nucleotides) to look for start and stop codons and finding the open reading frames in the sequence. '''
        startPositions = [] # create a list to store start positions of ORFS 
        theStart = 0 # this will indicate if a start codon has been found 
        theCodon = 0# this will indicate if any codon has been found 

        for frame in range(0, 3): # loop through the reading frames at 0, 1, and 2
            theStart = 0  # reset the notifyer for start codons 
            theCodon = 0 # reset the notifyer for any codon 
            startPositions = []  # reset the list of start positions for each reading frame 
            for i in range(frame, len(self.seq), 3): # loop through the sequence with the current reading frame 
                codon = self.seq[i:i + 3] # extract a codon (3 nucleotides) from the given sequence 
                if codon == 'ATG':  #check if the codon is a start codon (ATG)
                    startPositions.append(i) # if the codon is ATG then the index will be stored
                    theStart = 1 # if ATG found, the notifyer will be 1 which indicates that a start codon has been found 
                    theCodon = 1 # if ATG is found, the notifyer will be 1 which indicates that any codon has been found 

                if codon in OrfFinder.stopCodons and theStart: # check if the codin is a stop or start 
                    start = startPositions[0] + 1 - frame # calculate the length of the ORF from start and stop positions 
                    stop = i + 3 
                    length = stop - start + 1
                    self.saveOrf((frame % 3) + 1, start, stop, length) # save ORF info (frame, start position, stop position, length)
                    startPositions = [] # reset the start positions list 
                    theStart = 0 # reset the start codon notifyer 
                    theCodon = 1 # set any codon notifyer to 1 which means a stop codon was found after a start codon 

                if not theCodon and codon in OrfFinder.stopCodons: # if no codon has been found and the current codon is a stop codon this loop will run 
                    start = 1 # define start and stop positions 
                    stop = i + 3
                    length = stop - start + 1
                    self.saveOrf((frame % 3) + 1, start, stop, length)
                    # save the ORF information 
                    startPositions = [] # reset start positions list
                    theCodon = 1 # set codon notifyer to indicate a stop codon was found 

            if theStart:  # if a start codon has been found this loop will run 
                start = startPositions[0] + 1 # define start and stop positions 
                stop = len(self.seq)
                length = stop - start + 1
                self.saveOrf((frame % 3) + 1, start, stop, length)
                # save ORF information 

        return self.orfs # end the loop and method, returning the list of identified ORFs 

    def findRevOrfs(self): # method to get the reverse complement of the given DNA sequence 
        '''In this method, the reverse compliment will be found by reversing and swapping bases.
        Using the reverse compliment, I will then find the start and stop codons and ORFs and store them in a list. '''
        reverseComp = self.reverseComp() # initilze variables 
        startPositions = [] # similar to finding forward orfs, create list and set notifyers to 0 
        foundStart = 0
        foundCodon = 0

        for frame in range(0, 3): # loop through frames 0, 1, 2
            foundStart = 0 # reset notifyers for start and stop as well as reset the list of start positions 
            foundCodon = 0
            startPositions = []  
            for i in range(frame, len(reverseComp), 3): # loops through the reverse complement with the current reading frame 
                codon = reverseComp[i:i + 3]  #extract a codon (3 nucleotides) from the reverse compliment sequence 

                if codon == 'ATG': # check if the codon is a start codon (ATG)
                    startPositions.append(i) # if ATG then store the index of the start codon 
                    foundStart = 1 # indicates a start codon has been found 
                    foundCodon = 1 # indicates any codon has been found 

                if codon in OrfFinder.stopCodons and foundStart: # check if the codon is a stop or start codon 
                    stop = len(reverseComp) - startPositions[0] # calculate start and stop positions, and the length of the ORF 
                    start = len(reverseComp) - (i + 2) # calcualte start position based on the stop codon 
                    if frame == 1: stop += 1
                    elif frame == 2: stop += 2
                    length = stop - start + 1 # calculates the length of the ORF 
                    self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length) # save ORF information 
                    startPositions = [] # reset start position list 
                    foundStart = 0 # reset start codon notifyer 
                    foundCodon = 1 # indicates a stop codon was found after start codon 

                if not foundCodon and codon in OrfFinder.stopCodons:  # check if no start codon has been found and the current codon is a stop codon 
                    start = len(reverseComp) - i - 2 # define start and stop positions assuming the end of the sequence 
                    stop = len(reverseComp)
                    length = stop - start + 1
                    self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length) # save ORF information 
                    startPositions = [] # reset start positions list 
                    foundCodon = 1 # indicates a stop codon was found 

            if foundStart: # if a start codon is found 
                start = startPositions[0] + 1 # define start and stop positions assuming the start codon to the beggining of the sequence 
                stop = 1
                length = stop - start + 1
                self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length) # save ORF information 

        return self.orfs # return list of identified ORFs in the reverse compliment strand  

    def saveOrf(self, frame, start, stop, length): # define a method to save ORF information 
        self.orfs.append([frame, start, stop, length]) # append the ORF info to the list of ORFs 

    def reverseComp(self): # define a method to generate the reverse compliment of the DNA sequence 
            '''This method is where the reverse compliment will be created and used in the findRevOrfs method to find orfs in the compliment string. '''
            complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'} # dictionary that stores the compliment bases 
            reverseCompSeq = '' # create an empty string to store the reverse compliment string 
            for base in reversed(self.seq): # reverse the sequence and find bases in the reverse seq
                if base in complement: # if the bases are in the compliment 
                    reverseCompSeq += complement[base] # increment the new bases to the reverseCompSeq to create the reverse compliment strand 
                else:
                    reverseCompSeq += base # if the base is not found in the dictionary then skip it
            return reverseCompSeq


def main(inCommandLine=None): # define the main function 
    '''In the main function I will use the command line class as well as fasta to parse the file and properly read the fasta file. '''
    if inCommandLine is None: # check if command line arguments are provided 
        theCommandLine = CommandLine() # if there are no command line arguments, create a CommandLine instance to parse arguments 
        if theCommandLine.args.longestGene: # check if the 'longestgene' flag is set 
            fastaFile = FastAreader() # if flag is set, use fasta reader to read given fasta sequence 
            for header, sequence in fastaFile.readFasta(): # iterate over each header and sequence in the fasta file 
                print(header) # print the header 
                orfData = OrfFinder(sequence) # create orffinder object for the current sequence 
                orfData.findOrfs() # find orfs in the sequence and its reverse complement 
                orfData.findRevOrfs()
                filteredList = filter(lambda orf: orf[3] > theCommandLine.args.minGene, orfData.orfs)
                # filter the orfs based on the minimum gene length  
                for frame, start, stop, length in sorted(filteredList, key=lambda orf: orf[3], reverse=True):
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length))
                    # print the filtered ORFs sorted by length in descending order 

    else: # If command-line arguments are provided, create a CommandLine instance with those arguments
        theCommandLine = CommandLine(inCommandLine)
    print(theCommandLine.args) # Print the parsed command-line arguments

if __name__ == "__main__": # call the main function 
    main() 


# Python3 SierraDonnORFfinder.py  < tass2.fa > testoutputFile -mG 300 -lG -s ATG -s GTG
# type this into comand line to run 


# Python3 SierraDonnORFfinder.py  < lab5test.fa > extratestFile -mG 0 -lG -s ATG -s GTG
