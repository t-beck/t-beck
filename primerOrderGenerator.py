# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 16:15:57 2015

@author: batzer
"""
#%%
import re, sys

args = sys.argv[1:]

orderTemplate = """{0}\t{1}
{2}\t{3}
"""
#%%
argv1 = args[0]
argv2 = args[1]
argv3= args[2]

with open(argv1, 'r') as inFile:
    string = inFile.read()
    
#%%
primers = string.split('\n=\n')[:-1]
primerNo = 0
primersToWrite =[]
#%%
for primer in primers:
    try:
        chrmnum = re.search(r'(?<=chr)\d+', primer).group()
	startPoint=re.search(r'(?<=:)\d+(?=-)', primer).group()
        left_primer_sequence = re.search(r'(?<=Left_primer\t)[ACTGactg]+', primer).group()
        right_primer_sequence = re.search(r'(?<=right_primer\t)[ACTGactg]+', primer).group()
        left_primer_ID = '{orgID}_{0}F_{chrmnum}_{startPoint}'.format(primerNo, orgID=argv3, chrmnum=chrmnum, startPoint=startPoint)
        right_primerID = '{orgID}_{0}R_{chrmnum}_{startPoint}'.format(primerNo, orgID=argv3, chrmnum=chrmnum, startPoint=startPoint)
        primersToWrite.append([left_primer_ID, left_primer_sequence, right_primerID, right_primer_sequence])
        primerNo += 1

    except:
        pass
#%%
with open(argv2, 'w') as outFile:
    for primerToWrite in primersToWrite:
        outFile.write(orderTemplate.format(*primerToWrite))
#%%
