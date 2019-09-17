import os, sys, logging

from pdb import set_trace

import pandas as pd

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import SplitSampleInfo, AttachSeqToIndel, Helper

logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                    level=logging.INFO)


strProjectFile    = sys.argv[1]
strUserName       = sys.argv[2]
strProjectName    = sys.argv[3]

"""
/media/hkim/Pipeline/Indel_searcher_2/miniconda2/bin/python ./Indel_normalization.py User/JaeWoo/JaeWoo_test_samples.txt JaeWoo JaeWoo_test_samples
"""

def MakeIndelSeqDict():

    """
    dictD0Indel:
    {'sample_*': {'ACGATCGAT': {'Total': 300, {'ACGATCGAT_30M2I_AG': {'IndelCount': 3}}}}}}}

    validation
    ./Output/JaeWoo/JaeWoo_test_samples/190819_Nahye_12K_D7_2_D0_1-Cas9D7/Tmp
    grep TTTGGATCGTCTATCGTCG 190819_Nahye_12K_D7_2_D0_1-Cas9D7_Indel_freq.txt | grep 18M16D | wc -l
    -> Indel count
    """

    dictD0Indel  = {}
    dictExpIndel = {}

    with open(strProjectFile) as SampleList:

        for strSample in SampleList:
            print(strSample)

            tupSampleInfo = SplitSampleInfo(strSample)
            if not tupSampleInfo: continue
            strSample, strRef, strExpCtrl = tupSampleInfo

            if strExpCtrl == 'CTRL':
                dictD0Indel[strSample] = {}
            elif strExpCtrl == 'EXP':
                dictExpIndel[strSample] = {}

            with open('./Output/{user}/{project}/{sample}/Tmp/{sample}_Indel_freq.txt'.format(
                      user=strUserName, project=strProjectName, sample=strSample)) as IndelFreq,\
                open('./Output/{user}/{project}/{sample}/Result/{sample}_Summary_result.tsv'.format(
                      user=strUserName, project=strProjectName, sample=strSample)) as TotalResult:

                for strRow in IndelFreq:
                    listCol      = strRow.replace('\n','').split('\t')
                    strBarcode   = listCol[0]
                    strIndelPos  = listCol[2]
                    strRefseq    = listCol[4]
                    strQueryseq  = listCol[5]

                    if strExpCtrl == 'CTRL':
                        AttachSeqToIndel(strSample, strBarcode, strIndelPos, strRefseq, strQueryseq, dictD0Indel)
                    elif strExpCtrl == 'EXP':
                        AttachSeqToIndel(strSample, strBarcode, strIndelPos, strRefseq, strQueryseq, dictExpIndel)

                TotalResult.readline() ## skip header
                for strRow in TotalResult:
                    listCol    = strRow.replace('\n', '').split('\t')
                    strBarcode = listCol[0]
                    intTotal   = int(listCol[1])

                    try:
                        dictD0Indel[strSample][strBarcode]['Total'] = intTotal
                    except KeyError:
                        pass

                    try:
                        dictExpIndel[strSample][strBarcode]['Total'] = intTotal
                    except KeyError:
                        pass

        #set_trace()
        #print(dictSub.items())#
    return (dictD0Indel, dictExpIndel)


def MakeTmp(dictD0Indel, dictExpIndel):

    for dictIndel in [dictD0Indel, dictExpIndel]:
        for strSample, dictBarcode in dictIndel.items():
            strTmpDir = './Output/{user}/{project}/{sample}/Tmp'.format(user=strUserName,
                                                                        project=strProjectName,
                                                                        sample=strSample)
            with open(os.path.join(strTmpDir, strSample+'_indel_seq_count.txt'), 'w') as Output:
                for strBarcode, dictCountTotalAndIndel in dictBarcode.items():
                    for strIndelSeq, dictCount in dictCountTotalAndIndel.items():
                        if strIndelSeq == 'Total': continue
                        Output.write('\t'.join([strIndelSeq, str(dictCount['IndelCount'])])+'\n')


def MergeD0SampleResults(dictD0Indel):

    """
    dictD0Indel:
    {'sample_*': {'ACGATCGAT': {'Total': 300, {'ACGATCGAT_30M2I_AG': {'IndelCount': 3}}}}}}}

    -> sum total, sum indelcount

    dictD0IndelMerge:
    {'ACGATCGAT': {'Total': 600, {'ACGATCGAT_30M2I_AG': {'IndelCount': 5}}}}}}}
    """

    dictD0IndelMerge = {}

    for strD0SampleName in dictD0Indel:
        for strBarcode, dictCountTotalAndIndel in dictD0Indel[strD0SampleName].items():

            try:
                dictD0IndelMerge[strBarcode]['Total'] += dictCountTotalAndIndel['Total']
            except KeyError:
                dictD0IndelMerge[strBarcode] = {}
                dictD0IndelMerge[strBarcode]['Total'] = dictCountTotalAndIndel['Total']

            for strIndelSeq, dictCount in dictCountTotalAndIndel.items():  ## dcitCount : {'TTTGAGCATATCACACGAT:33M1D_T': {'IndelCount': 0}}
                if strIndelSeq == 'Total': continue

                try:
                    dictD0IndelMerge[strBarcode][strIndelSeq]['IndelCount'] += dictCount['IndelCount']
                except KeyError:
                    dictD0IndelMerge[strBarcode][strIndelSeq] = {}
                    dictD0IndelMerge[strBarcode][strIndelSeq]['IndelCount'] = dictCount['IndelCount']

    return dictD0IndelMerge


def SubtractIndelWithD0(dictD0IndelMerge, dictExpIndel):

    """
    dictD0IndelMerge: indel proportion - dictExpIndel: indel proportion
    """
    strD0SubResultDir = './Output/{user}/{project}/All_results/D0SubResult'.format(user=strUserName, project=strProjectName)
    Helper.MakeFolderIfNot(strD0SubResultDir)

    for strSample, dictBarcode in dictExpIndel.items():
        with open(os.path.join(strD0SubResultDir, '{sample}_D0SubResult.txt').format(sample=strSample), 'w') as Output:
            Output.write('Barcode_indel_seq\tD0_total\tD0_indel_prop\tExp_total\tExp_indel_prop\tD0_sub_indel_prop\n')

            for strBarcode, dictCountTotalAndIndel in dictBarcode.items():

                intExpTotal = dictCountTotalAndIndel['Total']

                for strIndelSeq, dictCount in dictCountTotalAndIndel.items():
                    if strIndelSeq == 'Total': continue

                    try:
                        intD0Total = dictD0IndelMerge[strBarcode]['Total']
                        intD0Count = dictD0IndelMerge[strBarcode][strIndelSeq]['IndelCount']

                        floD0Prop  = round(intD0Count / float(intD0Total), 6)

                        intExpCount = dictCount['IndelCount']
                        floExpProp  = round(intExpCount / float(intExpTotal), 6)

                        floSubExpIndel = floExpProp - floD0Prop
                        if floSubExpIndel < 0:
                            floSubExpIndel = 0

                        Output.write('\t'.join(map(str, [strIndelSeq,intD0Total, floD0Prop,
                                                         intExpTotal, floExpProp, floSubExpIndel]))+'\n')
                    except KeyError:
                        intExpCount = dictCount['IndelCount']
                        floExpProp  = round(intExpCount / float(intExpTotal), 6)

                        Output.write('\t'.join(map(str, [strIndelSeq, 'None', 'None',
                                                         intExpTotal, floExpProp, floExpProp]))+'\n')


def Main():
    logging.info("Indel normalization Start")
    logging.info("MakeIndelSeqDict")
    dictD0Indel, dictExpIndel = MakeIndelSeqDict()
    logging.info("MakeTmp")
    MakeTmp(dictD0Indel, dictExpIndel)
    logging.info("MergeD0SampleResults")
    dictD0IndelMerge = MergeD0SampleResults(dictD0Indel)
    logging.info("SubtractIndelWithD0")
    SubtractIndelWithD0(dictD0IndelMerge, dictExpIndel)
    logging.info("Indel normalization End")


Main()

