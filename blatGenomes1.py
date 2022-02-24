# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 16:12:12 2015

@author: batzer
"""
from __future__ import print_function
import Tkinter,tkFileDialog
import subprocess, os, sys, time, re
from Bio import SeqIO
from subprocess import Popen, PIPE, STDOUT

root = Tkinter.Tk()
root.withdraw()

def revComp(seq1):   #generates reverse complement
        
    revComp = re.sub('(?i)A', '~',seq1)
    revComp = re.sub('(?i)T', 'A', revComp)
    revComp = re.sub('~', 'T', revComp)
    revComp = re.sub('(?i)C', '~', revComp)
    revComp = re.sub('(?i)G', 'C', revComp)
    revComp = re.sub('~', 'G', revComp)
    
    revComp=list(revComp)

    revComp=''.join(revComp)
            
    return revComp[::-1] #reverse string

def prompt4DatabasesandOoc(databasesAndOocFiles,  is_additional_run=False):
    prompt1 = '\nType "add" to add a database and its associated .ooc file OR press ENTER to continue\n\n>>> '
    prompt2 = '\nType "add" AGAIN to select an ADDITIONAL database and its associated .ooc file OR press ENTER to continue\n\n>>> '

    if is_additional_run:
        prompt=prompt2
    
    else:
        prompt=prompt1
        
    choice=raw_input(prompt)
        
    if choice == 'add':
        database = tkFileDialog.askopenfilename(parent=root,title='Choose Database (.2bit)')
        oocfile = tkFileDialog.askopenfilename(parent=root,title='Choose associated oocFile (.ooc)')
        databasesAndOocFiles.append((database,oocfile))
        databasesAndOocFiles = prompt4DatabasesandOoc(databasesAndOocFiles, is_additional_run=True)
    return databasesAndOocFiles

def blat(fasta, blatOutFile, database, ooc):
    subprocess.call(['blat', database, fasta, 'ooc=' + ooc, blatOutFile])
    print('Done with {0}'.format(database.split('/')[-1].split('.')[0]))
    
def pslFilter(blatOutFile, filteredOutFile):
    subprocess.call(['pslFilter', '-noHead', blatOutFile, filteredOutFile])
    with open(filteredOutFile, 'r') as inFile:
        newLines = []
        #fill mateLocs
        previousID = None
        for line in inFile.readlines():
            ID = line.split('\t')[9]
            if previousID != ID:
                newLines.append(line)
            previousID = ID
    with open(filteredOutFile, 'w') as outFile:
        outFile.writelines(newLines)
    
def seqPull(dbName, database, filteredOutFile):
    Queries = []
    Seqs = []
    Masks=[]
    SeqOutput = os.getcwd() + '/temp.fa'
    with open(filteredOutFile, 'r') as inFile:
        for line in inFile.readlines():
            columns = re.split(r'\s+', line)
            query = columns[9]
            chrom = columns[13]
            start = columns[15]
            end = columns[16]
            orientation = columns[8]
            
            Queries.append(query)
            coordinates = [chrom, start, end]
            subprocess.call(['twoBitToFa', '{0}:{1}:{2}-{3}'.format(database, *coordinates), SeqOutput])
            with open(SeqOutput, 'r') as SeqOut:     
                SeqOut=SeqOut.read()
            SeqMask=[]
            for character in (''.join(SeqOut.split('\n')[1:])):
                if character.isupper():
                    SeqMask.append('1')
                else:
                    SeqMask.append('0')
            if orientation=='-':
                reverseComplement = revComp(''.join(SeqOut.split("\n")[1:]))
                HeaderLine=SeqOut.split('\n')[0]
                Seqs.append(HeaderLine+'\n'+reverseComplement+'\n')
                SeqMask=SeqMask[::-1]
                       
            else:
                Seqs.append(SeqOut)
                        
            Masks.append(''.join(SeqMask))                
    return [dbName, Seqs, Masks, Queries]
    
def checkSeqsToAlign(seqsToAlign):
    organismSeqIDs = []
    differences = []
    removalSet = set()
    masterSet = [seq.split('\n')[0][1:] for seq in seqsToAlign[0][1]]
    masterSet = set(masterSet)
    for organism in seqsToAlign[1:]:
        organismSeqIDs.append(set(organism[3]))
    for organismSet in organismSeqIDs:
        differences.append(masterSet ^ organismSet)
    for difference in differences:
        print(difference)
        removalSet = removalSet | difference
    print(removalSet)
    for seq in reversed(seqsToAlign[0][1]):  #reversed to avoid errors editing lists in place
        if seq.split('\n')[0][1:] in removalSet:
            del seqsToAlign[0][2][seqsToAlign[0][1].index(seq)]
            del seqsToAlign[0][1][seqsToAlign[0][1].index(seq)]
            
    for i in range(1, len(seqsToAlign)):
        for seq in reversed(seqsToAlign[i][3]):  #reversed to avoid errors editing lists in place
            if seq in removalSet:
                del seqsToAlign[i][2][seqsToAlign[i][3].index(seq)]
                del seqsToAlign[i][1][seqsToAlign[i][3].index(seq)]
    return seqsToAlign
                
        
    
def alignSeqs(seqsToAlign):
            
    organismCount = len(seqsToAlign)
    print(len(seqsToAlign[0][1]), len(seqsToAlign[1][1]), len(seqsToAlign[2][1]), len(seqsToAlign[3][1]))
    alignments = []
    for i in range(len(seqsToAlign[0][1])):
        fastaToAlign = []
        passed_length_check = True
        sequencelengths = []
        for j in range(organismCount):
            seq = seqsToAlign[j][1][i]
            sequencelengths.append(len(seq))
        for seqLength in sequencelengths:
            if seqLength > 5000:
                passed_length_check = False
        if passed_length_check == True:
            for j in range(organismCount):
                dbName = seqsToAlign[j][0]
                seq = seqsToAlign[j][1][i].replace('>', '>{0}_{1}_'.format(j, dbName)) #Adds database name to seq for identification purposes later on
                fastaToAlign.append(seq.strip())
            p = Popen(['muscle3.8.31_i86linux64'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
            alignment = p.communicate('\n'.join(fastaToAlign))[0]
            alignedseqs = alignment.replace('>','|>').split('|')[1:]
            alignedseqs.sort(key = lambda x: int(x[1]))
            for j in range(organismCount):
                Mask=seqsToAlign[j][2][i]
                seqHead = alignedseqs[j].split('\n')[0]
                seqSequence = list(''.join(alignedseqs[j].split('\n')[1:]))
                sequenceIndex=0
                maskIndex = 0
                for character in Mask:
                    while seqSequence[sequenceIndex] == '-':
                        sequenceIndex += 1
                    if character == '0':
                        seqSequence[sequenceIndex] = seqSequence[sequenceIndex].lower()
                    sequenceIndex += 1
                    maskIndex += 1
                
                maskedSeq = seqHead + '\n' + ''.join(seqSequence)
                alignedseqs[j] = maskedSeq
            alignment = '\n'.join(alignedseqs)
            alignments.append(alignment)
        
    return alignments
        
            
def seqCount(fasta):
    with open(fasta, 'r') as infile:
        seqCount=infile.read().count('>')
        progress_increment=seqCount
    return progress_increment
    
    
def main(fasta, outFile):
    seqsToAlign = []
    progress_increment=seqCount(fasta)
    OOIname = raw_input('\nPlease enter the name of your organism of interest\n\n>>>')
    with open(fasta, 'r') as OOIfasta:
        OOIseqs = OOIfasta.read().replace('>', '|>').split('|')[1:]
        OOImasks = []
        for seq in OOIseqs:
            SeqMask = []
            for character in (''.join(seq.split('\n')[1:])):
                    if character.isupper():
                        SeqMask.append('1')
                    else:
                        SeqMask.append('0')
            OOImasks.append(''.join(SeqMask))
        seqsToAlign.append([OOIname, OOIseqs, OOImasks])
    databasesAndOocFiles = []
    databasesAndOocFiles = prompt4DatabasesandOoc(databasesAndOocFiles)
    print('Currently blatting!')
    print('Sequences to process: ', progress_increment)
    for databasesAndOocFile in databasesAndOocFiles:
        dbName = databasesAndOocFile[0].split('/')[-1].split('.')[0]
        blatOutFile = '{0}/{1}.{2}.psl'.format(os.getcwd(), fasta.split('/')[-1], dbName)
        filteredOutFile = blatOutFile + '.filtered'
        blat(fasta, blatOutFile, *databasesAndOocFile)
        time.sleep(1)
        pslFilter(blatOutFile, filteredOutFile)
        seqsToAlign.append(seqPull(dbName, databasesAndOocFile[0], filteredOutFile))
    checkSeqsToAlign(seqsToAlign)
    alignments = alignSeqs(seqsToAlign)
    with open(outFile, 'w') as outFile:
        outFile.write('\n=\n'.join(alignments))
    
        
    
if __name__ == "__main__":
    main(*sys.argv[1:])