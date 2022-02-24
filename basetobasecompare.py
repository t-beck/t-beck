import Tkinter,tkFileDialog
import re, types
from subprocess import Popen, PIPE, STDOUT

root = Tkinter.Tk()
root.withdraw()

fasta = tkFileDialog.askopenfilename(parent=root,title='Choose comparison file')
try:
    mispriminglibrary = tkFileDialog.askopenfilename(parent=root,title='OPTIONAL: Choose mispriming library')
except:
    mispriminglibrary = None
def compareAndGeneratePrimers(sequences, mispriminglibrary):
    dbNames = []
    seqs=[]
    headers = []
    for sequence in sequences:
        sequence = sequence.strip()
        header = sequence.split('\n')[0]
        headers.append(header)
        dbName = re.search(r'(?<=_)[^_]+(?=_)', header).group()
        dbNames.append(dbName)
        seq = ''.join(sequence.split('\n')[1:]).upper()
        seqs.append(list(seq))
    
    consecMuts=[]
    prevTruthTable=list(range(len(seqs)))
    for i in range(len(seqs[0])):
        truthTable = []
        for j in range(len(seqs)):
            if seqs[j][i] == seqs[0][i]:
                truthTable.append(True)
            else:
                truthTable.append(False)
            if truthTable[j]==False and truthTable[j]==prevTruthTable[j] and seqs[0][i] != '-':
                consecMuts.append(i)
                consecMuts.append(i-1)
        if truthTable.count(False) < 2:
            seqs[0][i] = seqs[0][i].upper()
    
        else:
            seqs[0][i] = seqs[0][i].lower()
    
        prevTruthTable=truthTable
        
    for i in consecMuts:
        seqs[0][i]='N'
    
    compString=seqs[0]
    nonGaps = [character for character in compString if character != '-']
    compStringIndex = 0
    nonGapsToCompStringMap = []
    for i in range(len(nonGaps)):
        while nonGaps[i] != compString[compStringIndex]:
            compStringIndex += 1
        else:
            nonGapsToCompStringMap.append([i, compStringIndex])
            compStringIndex += 1
    TrueMidPoint = int(max([x[0] for x in nonGapsToCompStringMap]) / 2)
    TrueTargetStart, TrueTargetEnd = TrueMidPoint - 325, TrueMidPoint + 325
    compStringTargetStart, compStringTargetEnd = nonGapsToCompStringMap[TrueTargetStart][1], nonGapsToCompStringMap[TrueTargetEnd][1]
    targetLength = compStringTargetEnd - compStringTargetStart
    targetRegion = [compStringTargetStart, targetLength]
    
    primer3Input = ''.join(compString).replace('-', 'N')
    qScores = []
    
    for character in primer3Input:
        if character.isalnum():
            if character.islower():
                qScores.append('-1 ')
            else:
                qScores.append('0 ')
    qScores=''.join(qScores).strip()
    
    Seq_ID=headers[0]
    Seq_Temp=primer3Input
    Seq_Qual=qScores
    if type(mispriminglibrary) == types.StringType:
        library = (mispriminglibrary, '0.0')
    else:
        library = ('','0.0')
    
    
    Primer3Vars='''SEQUENCE_ID={0}
SEQUENCE_TEMPLATE={1}
SEQUENCE_QUALITY={2}
PRIMER_MISPRIMING_LIBRARY={library[0]}
PRIMER_WT_LIBRARY_MISPRIMING={library[1]}
PRIMER_MIN_QUALITY=-2
PRIMER_QUALITY_RANGE_MIN=-1000
PRIMER_LOWERCASE_MASKING=1
PRIMER_GC_CLAMP=1
PRIMER_MAX_POLY_X=4
PRIMER_MIN_END_QUALITY=0
PRIMER_QUALITY_RANGE_MAX=1
PRIMER_WT_SEQ_QUAL=0
SEQUENCE_TARGET={3},{4}
PRIMER_WT_END_QUAL=0
='''.format(Seq_ID, Seq_Temp, Seq_Qual, *targetRegion, library=library)
    
    p = Popen(['/home/batzer/Downloads/primer3-2.3.6/src/primer3_core', '-p3_settings_file=/home/batzer/Downloads/primer3-2.3.6/primer3web_v4_0_0_default_settings.txt'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    
    Primer3Output= p.communicate(Primer3Vars)[0]
   #print(Primer3Output)
    global primerCount
    try:
        seq_ID = re.search(r'(?<=SEQUENCE_ID=)[^\n]+', Primer3Output).group()
        print('_seq_ID\t' + seq_ID)
        left_pos = int(re.search(r'(?<=PRIMER_LEFT_0=)\d+(?=,)', Primer3Output).group())
        right_pos = int(re.search(r'(?<=PRIMER_RIGHT_0=)\d+(?=,)', Primer3Output).group()) + 1
        
        left_primer = re.search(r'(?<=PRIMER_LEFT_0_SEQUENCE=)[actgACTG]+', Primer3Output).group()
        print('Left_primer\t' + left_primer)       
        right_primer = re.search(r'(?<=PRIMER_RIGHT_0_SEQUENCE=)[actgACTG]+', Primer3Output).group()
        print('right_primer\t' + right_primer)
        left_TM = re.search(r'(?<=PRIMER_LEFT_0_TM=)[\d\.]+', Primer3Output).group()
        print('left_TM\t' + left_TM)
        right_TM = re.search(r'(?<=PRIMER_RIGHT_0_TM=)[\d\.]+', Primer3Output).group()
        print('right_TM\t' + right_TM)        
        left_GC = re.search(r'(?<=PRIMER_LEFT_0_GC_PERCENT=)[\d\.]+', Primer3Output).group()
        print('left_GC\t' + left_GC)        
        right_GC = re.search(r'(?<=PRIMER_RIGHT_0_GC_PERCENT=)[\d\.]+', Primer3Output).group()
        print('right_GC\t' + right_GC)        
        product_size = re.search(r'(?<=PRIMER_PAIR_0_PRODUCT_SIZE=)[\d]+', Primer3Output).group()
        for j in range(len(seqs)):
            gapCount = seqs[j][left_pos:right_pos].count('-')
            species_product_size = int(product_size) - gapCount
            print('PRODUCT_SIZE_{0}\t{1}'.format(dbNames[j], species_product_size))
        primerCount += 1

    except:
        print('FAILED')
        pass
    print('=')
    
primerCount = 0
    
with open(fasta, 'r') as fasta:
    alignments = fasta.read().split('\n=\n')
    
for alignment in alignments:
    sequences = alignment.replace('>', '|>').split('|')[1:]
    compareAndGeneratePrimers(sequences, mispriminglibrary)
print(primerCount)




