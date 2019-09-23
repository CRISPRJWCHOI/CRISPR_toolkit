#!/home/hkimlab/anaconda2/bin/python2.7



import sys
import subprocess as sp


sFile_path = sys.argv[1]

iSplit_line = int(sys.argv[2]) #400000
iSplit_num  = int(sys.argv[3]) #11

def Split():

    with open(sFile_path) as fq:

        for num in range(1, iSplit_num+1):
            with open('%s_%s.fq' % (sFile_path, num), 'w') as out:
                iCount = 0
                for sRow in fq:
                    iCount += 1
                    out.write(sRow)
                    if iCount == iSplit_line:
                        break         


def Make_filelist():
    
    with open('./LongGuide_Synthetic_2nd.txt', 'w') as filelist:
        
        for sFilename in sp.check_output('ls', shell=True).split('\n'):

            lFilename = sFilename.split('.')
            #print(lFilename)  
            if lFilename[-1] == 'fq':
                filelist.write(sFilename+'\n')


#Split()
Make_filelist()
