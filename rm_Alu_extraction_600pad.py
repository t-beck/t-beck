#!/usr/bin/env python3

import sys

with open(sys.argv[1],'r') as rpmsk:
    with open(sys.argv[2],'w') as bed:
        for line in rpmsk:
            spl=line.split()
            if len(spl)>11:
                spl_replace11a=spl[11].replace('(','-')
                spl_replace11b=spl_replace11a.replace(')','')
                spl_replace13a=spl[13].replace('(','-')
                spl_replace13b=spl_replace13a.replace(')','')
                if line.split()[10]=='SINE/Alu':
                    orient=spl[8]
                    if orient=='C':
                        orient='-'
                    if int(spl_replace11b)>=1 and not int(spl_replace11b)>=5 and int(spl[12])>=267:
                        start,end=str(int(spl[5])-600),str(int(spl[6])+600)
                        if int(start)>=1:
                            bed.write(spl[4]+'\t'+start+'\t'+end+'\t'+spl[9]+'_'+'rmsk_Alurange='+spl[4]+':'+start+'-'+end+"/5'pad=600/3'pad=600/strand="+orient+'\t'+'.'+'\t'+orient+'\n')
                        if int(start)<1:
                            start=1
                            bed.write(spl[4]+'\t'+str(start)+'\t'+end+'\t'+spl[9]+'_'+'rmsk_Alurange='+spl[4]+':'+str(start)+'-'+end+"5'pad=600/3'pad=600/strand="+orient+'\t'+'.'+'\t'+orient+'\n')
                    if int(spl_replace13b)>=1 and not int(spl_replace13b)>=5 and int(spl[12])>=267:
                        start,end=str(int(spl[5])-600),str(int(spl[6])+600)
                        if int(start)>=1:
                            bed.write(spl[4]+'\t'+start+'\t'+end+'\t'+spl[9]+'_'+'rmsk_Alurange='+spl[4]+':'+start+'-'+end+"/5'pad=600/3'pad=600/strand="+orient+'\t'+'.'+'\t'+orient+'\n')
                        if int(start)<1:
                            start=1
                            bed.write(spl[4]+'\t'+str(start)+'\t'+end+'\t'+spl[9]+'_'+'rmsk_AluRange='+spl[4]+':'+str(start)+'-'+end+"/5'pad=600/3'pad=600/strand="+orient+'\t'+'.'+'\t'+orient+'\n')
