#!/home/hkimlab/anaconda2/bin/python2.7

import os
import sys
import pdb
from datetime import datetime
from collections import namedtuple as nt
from collections import OrderedDict

sOutput_dir = sys.argv[1]

def Calculate_indel_freq():

    if not os.path.isdir('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir))

    for sFile in os.listdir('{outdir}/result/freq'.format(outdir=sOutput_dir)):
        #print sFile
        if os.path.isfile(os.path.join('{outdir}/result/freq'.format(outdir=sOutput_dir), sFile)):
            with open(os.path.join('{outdir}/result/freq'.format(outdir=sOutput_dir), sFile)) as Input_freq,\
                open(os.path.join('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir), sFile), 'w') as Output_freq:

                sRef       = Input_freq.readline()  # first row is ref.
                sDelemiter = Input_freq.readline()  # second row is '-------' delemiter.
                Output_freq.write(sRef+sDelemiter)

                lSeq_indel = []   # [namedtuple1(['TGCA', '30M3I']) namedtuple2 ...
                dFreq_count = {}   # {'30M3I':2 ... }

                for sRow in Input_freq:
                    Seq_indel = nt('Seq_indel', ['seq', 'indel', 'freq', 'ref_needle', 'query_needle'])

                    if sRow == sRef: continue
                    if sRow[0] == '-': continue
                    
                    try:
                        lCol   = sRow.replace('\n', '').split('\t')
                        Seq_indel.seq          = lCol[0]
                        Seq_indel.indel        = lCol[1]
                        Seq_indel.ref_needle   = lCol[3]
                        Seq_indel.query_needle = lCol[4]
                        lSeq_indel.append(Seq_indel)
                    except IndexError:
                        print sFile, lCol                        
                        continue
                        
                    try:
                        dFreq_count[Seq_indel.indel] += 1
                    except KeyError:
                        dFreq_count[Seq_indel.indel] = 1
                #end: for sRow

                # Add freq infomation pre-result data.
                lResult = []
                iTotal = len(lSeq_indel)

                #print 'dFreq_count', dFreq_count
                #print 'lSeq_indel', lSeq_indel

                for Seq_indel in lSeq_indel:
                    iCount = dFreq_count[Seq_indel.indel]
                    Seq_indel.freq  = float(iCount) / iTotal
                    lResult.append(Seq_indel)

                lResult.sort(key=lambda x: x.indel)
                lResult.sort(key=lambda x: x.freq, reverse=True)

                #print 'lResult', lResult

                for Seq_indel in lResult:
                    #print Seq_indel.__dict__
                    Output_freq.write('\t'.join(map(str, [Seq_indel.seq, Seq_indel.indel, Seq_indel.freq, Seq_indel.ref_needle, Seq_indel.query_needle]))+'\n')
            #end: with open
        #end: if os.path
    #end: sFile


def Make_indel_summary():

    lOutput = []

    for sFile in os.listdir('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir)):
        if os.path.isfile(os.path.join('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir), sFile)):
            with open(os.path.join('{outdir}/result/freq/freq_result'.format(outdir=sOutput_dir), sFile)) as Input_freq:

                sRef       = Input_freq.readline()  # first row is ref.
                sDelemiter = Input_freq.readline()  # second row is '-------' delemiter.

                dINDEL = OrderedDict()

                lTable = [sRow.replace('\n', '').split('\t') for sRow in Input_freq]
                iTotal = len(lTable)

                for lCol in lTable:
                    sINDEL = lCol[1]
                    try:
                        dINDEL[sINDEL] += 1
                    except KeyError:
                        dINDEL[sINDEL] = 1

                dINDEL = OrderedDict(sorted(dINDEL.items(), key=lambda t: t[1], reverse=True))
                
                llINDEL = [[sKey, iValue, round(iValue/float(iTotal),3)*100] for sKey, iValue in dINDEL.items()]
                sINDEL_result = ''.join([':'.join(map(str, lINDEL))+', ' for lINDEL in llINDEL])[:-2]

                lOutput.append([sFile, iTotal, sINDEL_result])
                #Output_freq.write('\t'.join([sFile, sTotal, sINDEL_result]) + '\n')

    lOutput = sorted(lOutput, key=lambda x: x[1], reverse=True)

    with open('{outdir}/result/freq/freq_result/Indel_summary.txt'.format(outdir=sOutput_dir), 'w') as Output_freq:
        for lCol in lOutput:
            Output_freq.write('\t'.join(map(str, lCol)) + '\n')


if __name__ == '__main__':
    print 'Indel frequency calculator start: ', datetime.now()
    Calculate_indel_freq()
    Make_indel_summary()
    print 'Indel frequency calculator end: ', datetime.now()
