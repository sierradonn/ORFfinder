class ProteinParam :

    ''' The following program will calculate the pysical-chemical properties of a protein sequence'''

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

    def __init__ (self, protein): #added the protein attribute to be used for all the other functions
        self.protein = protein
        
        aaGroup = ''.join(protein).split() #concatenates characters in protein and splits them this removes whitespace
        self.proteinInput = ''.join(aaGroup).upper() #concatenates aa and converts to uppercase 
        self.aaComp = {} #this creates a dictionary that will store the counts of each amino acid entered from the input sequence 
        
        for aa in self.aa2mw.keys(): #loops through all the aa in the dictionary aa2mw
            self.aaComp[aa] = float(self.proteinInput.count(aa)) #will count each aa and store it in aaComp dictionary 
            
    '''the count function will count the number of amino acids in the input string 
        if you input VLSPADKTNVKAAW
        then the print result will be 14 '''

    def aaCount (self): #added the len method to count number of amino acids
        self.aaNum = 0 #initilizese the variable aaNum 
        for aa in self.protein:  #for the amino acids in the protein, loop through each one
            if aa in self.aa2mw: #if an amino acid is in the dictionary aa2wm do something to it 
                self.aaNum += 1 #what will be done is one will be added for each charcter in the dictionary 
        return self.aaNum #stores the count value in the variable aaNum to be used later 
    
    
    '''the composition function will call the function written in _init_ 
    this function will output percentages of aa found in the sequence ''' 

    def aaComposition (self) :
        return self.aaComp #returns dictionary set up in __init__ function 
   
    '''the pI and charge functions will be used to find the  theoretical isolelectric point 
    this will cahnge depending on ipnut and aa that are in the key 
    but an example is the sequence VLSPADKTNVKAAW will output Theoretical pI: 9.88 '''
    
    def pI (self):
        arbitraryCharge =  1000000 #allows for loop to run through infinite positive and update goodPH when any lower charge is found 
        goodPH = 0 
        thepH = 0

        while thepH < 14.01: #loop through all the cases where pH less than 14
            myCharge = self._charge_(thepH)#calculates the absolute value of the charge at currect pH 
            if abs(myCharge) < arbitraryCharge: #if the charge is less than the arbitrary charge 
                #which is set to a high number so the loop can iterate over a lot of values 
                #is this is true this means the current pH will result in a charge closer to zero 
                arbitraryCharge = abs(myCharge) #updates the arbitraryCharge with the current charge value 
                goodPH = thepH 
            thepH += 0.01 #makes sure the best pH value is down to 2 decimals 
        return goodPH
        
    def _charge_ (self, pH):
        posCharge = 0 #initilizes the variables 
        negCharge = 0
        nTerminous = (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)
        cTerminous = (10 ** pH) / (10 ** self.aaCterm + 10 ** pH)
        for aa in self.protein: #loops through the aa that are in the input protein
            if aa in self.aa2chargePos: #if aa are in the dictionary with the positive charges the keys will be searched 
                #countAminoAcids = self.aaComp[aa]
                posCharge += (10 ** self.aa2chargePos[aa]) / (10 ** self.aa2chargePos[aa] + 10 ** pH)
                #((self.aaComp[aa]) * (10 ** self.aa2chargePos[aa])/((10 ** self.aa2chargePos[aa]) + (10 ** pH))) 
                #nTerm = (10 ** self.aaNterm)/((10 ** self.aaNterm) + (10 ** pH))
                
            elif aa in self.aa2chargeNeg: #if aa are in the dictionary with the negative charges the keys will be searched 
                #countAminoAcids = self.aaComp[aa] 
                negCharge += (10 ** pH) / (10 ** self.aa2chargeNeg[aa] + 10 ** pH)
                #((self.aaComp[aa]) * (10 ** pH)/((10 ** self.aa2chargeNeg[aa]) + (10 ** pH)))
                #cTerm = (10 ** pH)/((10 ** self.aaCterm) + (10 ** pH)) 
        negCharge += cTerminous #adds charge to c-terminus
        posCharge += nTerminous #adds charge to n-terminus 

        netCharge = posCharge - negCharge
        return netCharge


    '''the molar Extinction function is used to find the molarExtinction and also the massExtinction 
    from the sequence VLSPADKTNVKAAW the out put for both are as follows molar Extinction coefficient: 5500.00
    mass Extinction coefficient: 3.67'''

    def molarExtinction (self): 
        if "W" in self.protein or "Y" in self.protein or "C" in self.protein: # check if the letters W,Y, or C are in the protein
            y = self.protein.count("Y") #initilize the variables to use in calculations 
            w = self.protein.count("W") #these are the count value for the letters y,w,c respectively
            c = self.protein.count("C")
            aaList = ["Y", "W", "C"] #make a list of the three amino acids that have absorbance at 280 nm 
            
            calculate = 0 #initilize the variable 
            for aa, count in zip(aaList, [y, w, c]): #using the zip function i combined elements from aaList and [y, w, c]
                #this creates a tuple to be able to connect the values in the dictionary to the count of letter int the sequence 
                calculate += count * self.aa2abs280.get(aa, 0) #for each aa found in absorbance dictionary
                #the product of letter count and absorbance for that specific letter are added to the variable calculate  
            return calculate
        else: #if the aa from the sequence can not be found zero will be returned 
            return 0
       
    
    def massExtinction (self):
        pass
        myMW =  self.molecularWeight() 
        return self.molarExtinction() / myMW if myMW else 0.0

    '''using the dictionary of aa weights this function will sum the averages of the aa in the 
    given sequence and also acount for water loss during bonding
    from VLSPADKTNVKAAW the Molecular Weight will be  1499.7''' 
    
    def molecularWeight (self):
        if self.aaNum == 0: #if there is no entry of sequence then the weight will be zero  
            print("0")
        else: 
            self.weight = 0 #initilize the variables 
            loss = self.aaNum - 1 #amount of amino acids minus one will give us the water lost in bonding 
            for aa in self.protein: #for the aa in the sequence 
                self.weight += self.aa2mw.get(aa, 0) #updates and adds the weight if the aa is found in the dictionary 
                #from the dictionary the molecular weights will be added 
                
            return self.weight - (loss * self.mwH2O) #this will store the massExtinction #loss times the molecular weight of water gives the molecular weight of water loss #must subtract the molecular weight of water loss from the weight of the sequence to get total weight 

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
   
    def __init__ (self, inString=''):
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
        
        '''Constructs the NucParams class.
        arguments should include a DNA/RNA sequence- the default is an empty string.'''
            # initilize dictionaries to store sequence information 
        self.sequenceData= {
                'validNuc': {'A', 'C', 'G', 'T', 'U', 'N'},
                'sequence': '',
                'nucComp': {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0, 'N': 0},
                'codonComp': {},
                'aaComp': {}
            }

            # put rnaCodonTable and dna CodonTable in here to be able to call them 
        self.rnaCodonTable = rnaCodonTable
        self.dnaCodonTable = dnaCodonTable

        self.addSequence(inString)
            # call the add sequence in order to handle the sequence 
        
    def addSequence (self, inSeq):
        '''adds new DNA/RNA to the existing data.
        arguments include a sequence that is a string either DNA or RNA  
        
        '''
        inSeq = inSeq.upper() 
        # update the data structure 
        for nucleotide in inSeq:
            if nucleotide in self.sequenceData['validNuc']: # checks if the uppercase version of nucleotides is in the list of valid nucleotides
                self.sequenceData['nucComp'][nucleotide] += 1 # if the nucloetide is in validNuc then the nucComp dictionary will increment the count of nucleotides
        self.sequenceData['sequence'] += inSeq # had to add this line to make sure the input seqeunce are added to the DNA seqeunce 

    def aaComposition(self):
        aaComp = {} # initilize the dictionary 
        for i in range(0, len(self.sequenceData['sequence']), 3): # iterate through the DNA sequence at every 3 positions 
            codon = self.sequenceData['sequence'][i:i+3].upper().replace('T', 'U') # this loop will extract each codon from the DNA sequence 
            if codon in self.rnaCodonTable: # checks if the codon is in the rna codon dictionary 
                aa1 = self.rnaCodonTable[codon] # if the codon is in the table, the matching aa from the table will be retrieved
                aaComp[aa1] = aaComp.get(aa1, 0) + 1 # will update the aaComp count if not in dictionary already 
        return aaComp


    def nucComposition(self):
        return self.sequenceData['nucComp'] # return the 'nucComp' dictionary that hold the valid nuceotides
        
    def codonComposition(self):
        codonComp = {}
        seq = self.sequenceData['sequence']
        for i in range(0, len(seq) - 2, 3): # loops through the sequence and extracts each codon 
            codon = seq[i:i + 3].upper().replace('T', 'U')
            if all(nuc in self.sequenceData['validNuc'] for nuc in codon) and 'N' not in codon:
                # checks if valid nucleotides are in sequence and not 'N'
                codonComp[codon] = codonComp.get(codon, 0) + 1
                # if true then the count of codons in codonComp dictionary will be incremented 
        return codonComp


    def nucCount(self):
        return sum(self.sequenceData['nucComp'].values()) # gets the counts from the nunComp dictionary and sums them 

