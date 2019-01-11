import sys

def count_amino(peptides,aminoacid):
    'Counts the amount of aminoacid in the peptide'
    text = ''.join(peptides)
    n_amino = text.count(aminoacid)
    return n_amino

    
def RK(peptides):
    'Receives a list of peptides and returns the proportion R/K'
    text = ''
    for peptide in peptides:
        text += peptide
    nR = float(text.count('R'))
    nK = float(text.count('K'))
    if nK == 0:
        result = 0
    elif nR == 0:
        result = 0
    else:
        result = nR/nK
    return str(result)


def evaluate(prote, fieldSize,value):
    'Evaluates the net charge and returns the value of the net charge and the stretch'
    Charge = value
    lenght = len(prote)
    last = len(prote)-fieldSize
    find = False
    peptides = []

    if lenght < fieldSize: #prote size < fieldsize
        text = ''
        charge = 0.0
        #read the next lenght caracteres
        for j in range(lenght):
            #hold the sequence of characteres
            text += prote[j]
            #compute the netCharge
            if prote[j] == 'R':
                charge += 1.0
            elif prote[j] == 'K':
                charge += 1.0
            elif prote[j] == 'D':
                charge += -1.0
            elif prote[j] == 'E':
                charge += -1.0
        if charge == Charge:
            result = RK([text]), [text]
        else:
            result = '0',[]
        return result

    #for each character (does not include the last fieldsize characteres)
    for i in range(len(prote)-fieldSize+1):
        text = ''
        charge = 0.0
        #read the next fieldSize caracteres
        for j in range(fieldSize):
            #hold the sequence of characteres
            text += prote[i+j]
            #compute the netCharge
            if text[j] == 'R':
                charge += 1.0
            elif text[j] == 'K':
                charge += 1.0
            elif text[j] == 'D':
                charge += -1.0
            elif text[j] == 'E':
                charge += -1.0
        if charge == Charge:
            find = True
            delta = i + fieldSize
            peptides.append(text)
            
    if find:
        result = RK(peptides), peptides
    else:
        result = '0',[] #does not find the charge
    return result


#main function
def main():
    Charges = range(31)
    names = []
    fieldSize = 30
    blockSize = 1000
    #read the program arguments
    if len(sys.argv) < 3:
        print "Use: python", sys.argv[0], "<input complete file name>", "<output root file name>"
        sys.exit(0)

    filein = sys.argv[1] #input file name
    rootfileout = sys.argv[2] #output root file name

    #read the file's contents
    fin = open(filein, 'r')
    contents = fin.read()
    
    #open the first block of the output file
    block = 1
    fileout  = rootfileout + '.txt'
    fout = open(fileout, 'w')
    #start the block size counter
    count = 0
    
    #split the file's contents by '>'
    proteins = str.split(contents, '>')[1:] #the file header is discarded

    #sequence
    seq_all = []
    for i in Charges:
        seq_all.append([])
        
    #handle each protein...
    for prot in proteins:
        seq = []
        lista = str.split(prot, '\n') #split each protein 
        names.append(lista[0]) #save the names
        
        #concatenate the sequences of amino acids of the protein
        for temp in lista[1:]:
            seq.append(temp[:])    
        seqfinal = ''.join(seq) 
        protName = str.split(lista[0]) #take the first name
        for value in Charges:
            _, peptides = evaluate(seqfinal, fieldSize,value)
            seq_all[value] += peptides
        print protName[0]
        
    for value in Charges:
        peptides = seq_all[value]
        nR = count_amino(peptides,'R')
        nK = count_amino(peptides,'K')
        rk = RK(peptides)
        fout.write('Charge: ' + str(value) + '\t' + 'R: ' + str(nR) + '\t' + 'K: ' + str(nK) + '\t' + 'R/K: ' + str(rk) + '\n')

    fin.close()
    fout.close()
    print 'Done!'        

main()

