#!/extdata1/JaeWoo/Tools/Python/miniconda2/bin/python2.7

import os, sys
from pdb import set_trace

import cPickle

#strSampleFolder = sys.argv[1]


def MakeIndelContrDict():

    for strSampleFolder in ['18K_D0_1','18K_D0_2','18K_D0_3']:
        with open('./Output/%s/%s_IndelSubtarction.txt' % (strSampleFolder, strSampleFolder), 'w') as Output:

            dictSub = {}

            for strFile in os.listdir('./Output/%s/result' % strSampleFolder):
                if 'filtered' in strFile:
                    with open('./Output/%s/result/%s' % (strSampleFolder, strFile)) as Input:

                        strBarcodeName = strFile.replace('_filtered_indel.txt','')
                        for strRow in Input:
                            listCol      = strRow.replace('\n','').split('\t')
                            #set_trace()
                            strIndelPos  = listCol[2].replace("['",'').replace("']",'')
                            listIndelPos = strIndelPos.split('M')
                            intMatch     = int(listIndelPos[0])
                            strRefseq    = listCol[4]
                            strQueryseq  = listCol[5]

                            if 'I' in strIndelPos: ## insertion
                                intInsertion    = int(listIndelPos[1].replace('I', ''))
                                strInsertseq    = strQueryseq[intMatch:intMatch+intInsertion]
                                #set_trace()
                                strInsertPosSeq = strIndelPos+'_'+strInsertseq

                                try:
                                    dictSub[strBarcodeName+':'+strInsertPosSeq].append([strInsertPosSeq, strRefseq, strQueryseq])
                                except KeyError:
                                    dictSub[strBarcodeName+':'+strInsertPosSeq] = [[strInsertPosSeq, strRefseq, strQueryseq]]

                            elif 'D' in strIndelPos:
                                intDeletion     = int(listIndelPos[1].replace('D', ''))
                                strDeleteSeq    = strRefseq[intMatch:intMatch+intDeletion]
                                strDeletePosSeq = strIndelPos+'_'+strDeleteSeq

                                try:
                                    dictSub[strBarcodeName+':'+strDeletePosSeq].append([strDeletePosSeq, strRefseq, strQueryseq])
                                except KeyError:
                                    dictSub[strBarcodeName+':'+strDeletePosSeq] = [[strDeletePosSeq, strRefseq, strQueryseq]]

            for strBarcodeName, list2IndelPosSeq in dictSub.items():
                for listIndelPosSeq in list2IndelPosSeq:
                    Output.write('\t'.join([strBarcodeName] + listIndelPosSeq) + '\n')


def ConcatContrDict():

    DictSubNoDup = {}

    for strSampleFolder in ['18K_D0_1', '18K_D0_2', '18K_D0_3']:
        with open('./Output/%s/%s_IndelSubtarction.txt' % (strSampleFolder, strSampleFolder)) as Input:

            for strRow in Input:
                listCol = strRow.replace('\n', '').split('\t')
                try:
                    DictSubNoDup[listCol[0]] += 1
                except KeyError:
                    DictSubNoDup[listCol[0]] = 1

    #print(DictSubNoDup)
    with open('./Output/DictSubNoDup.pickle', 'wb') as PickleObj:
        cPickle.dump(DictSubNoDup, PickleObj)


def Main():
    #MakeIndelContrDict()
    ConcatContrDict()


Main()