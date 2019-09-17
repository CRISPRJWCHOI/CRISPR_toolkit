import os, sys, logging

from pdb import set_trace
from datetime import datetime
from collections import OrderedDict
from collections import namedtuple as nt

strOutputDir = sys.argv[1]
strSample    = sys.argv[2]
strLogPath   = sys.argv[3]

logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                    level=logging.DEBUG,
                    filename=strLogPath,
                    filemode='a')
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def MakeIndelSummary():

    """
    Input
    TTTGCAGAGTATATCACACCATATCA AGTCAGACAAGGAGCACCACACGGTGGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA 17M1I 0.134 AGTCAGACAAGGAGCAC-ACACGGTGGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTC------- AGTCAGACAAGGAGCACCACACGGTGGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA

    0: barcode
    1: target region
    2: indel pos
    3: total freq
    4: ref seq
    5: query seq

    Output
    TTTGTCTCGTACACTCGTATGCTGCA 2 18M2D:1:50.0, 24M1I:1:50.0
    TTTGACATCTACAGTGTCTCTCCACA 2 22M1I:2:100.0
    """

    listOutput = []

    with open('{outdir}/Tmp/{sample}_Indel_freq.txt'.format(sample=strSample, outdir=strOutputDir)) as InputFreq,\
        open('{outdir}/Tmp/{sample}_Indel_summary.txt'.format(sample=strSample, outdir=strOutputDir), 'w') as OutputFreq:

        listTable = [strRow.replace('\n', '').split('\t') for strRow in InputFreq]
        intTotal  = len(listTable)

        #strBarcode = listCol[0]
        dictINDEL = OrderedDict({listCol[0]:OrderedDict({'Total':0}) for listCol in listTable}) ## {'TTTGACATCTACAGTGTCTCTCCACA': {22M1I : 2, ...}}

        for listCol in listTable:
            strBarcode = listCol[0]
            strIndel   = listCol[2]

            dictINDEL[strBarcode]['Total'] += 1

            try:
                dictINDEL[strBarcode][strIndel] += 1
            except KeyError:
                dictINDEL[strBarcode][strIndel] = 1

        #dictINDEL = OrderedDict(sorted(dictINDEL.items(), key=lambda t: t[1], reverse=True)) ## sort value count.

        list2Result = []
        for strBarcode in dictINDEL:
            intTotal       = dictINDEL[strBarcode]['Total']
            list2INDEL     = [[strIndel, intCount, round(intCount/float(intTotal),3)*100] for strIndel, intCount in dictINDEL[strBarcode].items()]
            list2INDEL     = sorted(list2INDEL, key=lambda x: x[1], reverse=True)
            strIndelResult = ''.join([':'.join(map(str, listINDEL))+', ' for listINDEL in list2INDEL if listINDEL[0] != 'Total'])
            list2Result.append([strBarcode, intTotal, strIndelResult])

        for listResult in sorted(list2Result, key=lambda x: x[1], reverse=True):
            OutputFreq.write('\t'.join(map(str, listResult)) + '\n')


if __name__ == '__main__':
    logging.info('Indel frequency calculator start: %s' % str(datetime.now()))
    MakeIndelSummary()
    logging.info('Indel frequency calculator end: %s' % str(datetime.now()))
