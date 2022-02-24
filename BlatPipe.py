# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:49:55 2015

@author: batzer
"""
###pulls coords from blat files 


import subprocess
import os
import sys


def pullcords(filename):
    coordinates=[]
    with open(filename, 'r') as infile:
        lines=infile.readlines()
        for line in lines:
            chrm=line.split(':')[0]
            startCoord=int(line.split(':')[1].split('-')[0])
            endCoord=int(line.split(':')[1].split('-')[1].strip())
            coordinates.append([chrm,str(startCoord),str(endCoord)])
    return coordinates
    
def blatSeqs(twobitGenome,coordinates):
    results=[]
    for coordinate in coordinates:
        seqAndRange=':{0}:{1}-{2}'.format(*coordinate)
        outpath=os.getcwd()+'temp.fa'
        subprocess.call(['twoBitToFa', twobitGenome+seqAndRange, outpath])
        with open(outpath, 'r') as resultfile:
            result=resultfile.read()
        os.remove(outpath)
        results.append(result)
    return results
    
def main(filenameCoords,twobitPath):
    
    coordinates=pullcords(filenameCoords)
    results=blatSeqs(twobitPath, coordinates)
    for result in results:
        print result.strip()


if __name__ == "__main__":
    main(*sys.argv[1:])
    
    
    
