#!/media/hkim/Pipeline/Indel_searcher_2/miniconda2/bin/python2.7

import os, re, sys, pickle
import subprocess as sp
from Bio import AlignIO
from pdb import set_trace


strFastq='/media/hkim/Pipeline/CRISPR_Indel_searcher/Input/FASTQ/190807_Nahye_24k_NG_rep1-24kLib/NG_rep1.extendedFrags.fastq'
strBarcode='TTTGGTGATCTCACTCTCGACAACTC'

sRef_fa = './Input/Reference/190807_Nahye_24k_NG_rep1-24kLib/Reference.fa'
sBarcode_PAM_pos='Foward'


def CountBar():
    with open(strFastq) as Input:
        
        intCnt=0

        for strRow in Input:
            if strBarcode in strRow:
                intCnt+=1

        print(intCnt)

def ExtractFastq():

    with open(strFastq) as Input,\
        open('./Input/FASTQ/Test1/Test1.fastq_target', 'w') as Output:

        listFastq  = []

        for i, strRow in enumerate(Input):
            listFastq.append(strRow.replace('\n', ''))
            if i % 4 == 3:
                #print(listFastq)
                if strBarcode in listFastq[1]:
                    Output.write('\n'.join(listFastq)+'\n')
                listFastq = []

            
def LoadPickle():
    
    with open('Output/Test1/Pickle/Test1.fastq_1.fq.pickle', 'rb') as Input:
        obj = pickle.load(Input)
        set_trace()


def CheckSearch():
    dRef = {}
    dResult = {}

    with open(sRef_fa) as Ref:

        iCount         = 0
        sBarcode       = ""
        sTarget_region = ""

        for sRow in Ref:
            iCount += 1

            if iCount % 2 != 0:
                #      barcode               target region
                # >CGCTCTACGTAGACA:CTCTATTACTCGCCCCACCTCCCCCAGCCC
                sBarcode_indel_seq = sRow.strip().replace('\n', '').replace('\r', '').split(':')
                sBarcode = sBarcode_indel_seq[0].replace('>', '')
                sTarget_region = sBarcode_indel_seq[1]

                ## Reverse the sentence. If it is done, all methods are same before work.
                if sBarcode_PAM_pos == 'Reverse':
                    sBarcode = sBarcode[::-1]
                    sTarget_region = sTarget_region[::-1]

            elif iCount % 2 == 0:
                ## Reverse
                sRef_seq = sRow.strip().replace('\n', '').replace('\r', '')

                if sBarcode_PAM_pos == 'Reverse':
                    sRef_seq = sRef_seq[::-1]

                Seq_matcher = re.compile(r'(?=(%s))' % sTarget_region)
                # iIndel_start_pos       = sRef_seq.index(sTarget_region)               # There is possible to exist two indel.
                iIndel_start_pos = Seq_matcher.finditer(sRef_seq)

                for i, match in enumerate(iIndel_start_pos):
                    iIndel_start_pos = match.start()
                # print iIndel_start_pos
                # print len(sTarget_region)
                # print sRef_seq
                iIndel_end_pos = iIndel_start_pos + len(sTarget_region) - 1

                try:
                    iBarcode_start_pos = sRef_seq.index(sBarcode)

                    #if iIndel_start_pos <= iBarcode_start_pos:
                    #    print(iIndel_start_pos, iBarcode_start_pos)
                    #    raise IndexError('indel is before barcode')

                    iBarcode_end_pos = iBarcode_start_pos + len(sBarcode) - 1
                    sRef_seq_after_barcode = sRef_seq[iBarcode_end_pos + 1:]

                    # modified. to -1
                    iIndel_end_next_pos_from_barcode_end = iIndel_end_pos - iBarcode_end_pos - 1

                    iIndel_start_next_pos_from_barcode_end = iIndel_start_pos - iBarcode_end_pos - 1

                    #  "barcode"-------------*(N) that distance.
                    #          ^  ^            ^
                    #   *NNNN*NNNN
                    #    ^    ^     indel pos, the sequence matcher selects indel event pos front of it.

                    dRef[sBarcode] = (sRef_seq, sTarget_region, sRef_seq_after_barcode, iIndel_start_next_pos_from_barcode_end,
                                      iIndel_end_next_pos_from_barcode_end, iIndel_start_pos,iIndel_end_pos)  # total matched reads, insertion, deletion, complex
                    dResult[sBarcode] = [0, 0, 0, 0, [], [], [], [], []]
                except ValueError:
                    continue

    with open('test.seq') as Input:
        
        iBarcode_matched = 0

        for sSeq in Input:
            sSeq = sSeq.replace('\n','') 

            listSeqWindow = [sSeq[i:i + 26] for i in range(len(sSeq))[:-25]]

            iNeedle_matched  = 0
            iInsert_count    = 0
            iDelete_count    = 0
            iComplex_count   = 0

            intFirstBarcode  = 0 ## check whether a barcode is one in a sequence.

            for strSeqWindow in listSeqWindow:

                if intFirstBarcode == 1: break ## A second barcode in a sequence is not considerable.

                try:
                    lCol_ref = dRef[strSeqWindow]
                    sBarcode = strSeqWindow
                    intFirstBarcode = 1
                except KeyError:
                    continue

                iBarcode_matched += 1

        print(iBarcode_matched)


def CheckNeedle():
    
    sBarcode = 'TTTGACTAGTCATCACTATAGCATAA'
    sRef_seq_after_barcode = 'TACAGTGTTTTTTTTTTTTCAGAGGAAGCTTGGCGTAACTAGATCT'
    sQuery_seq_after_barcode = 'TACAGTGTTTTTTTTTTTCAGAGGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAATA'

    sRef_seq    = r'<(echo -e ">{name}\n{seq}")'.format(name='Ref', seq=sRef_seq_after_barcode)
    sQuery_seq  = r'<(echo -e ">{name}\n{seq}")'.format(name='Query', seq=sQuery_seq_after_barcode)

    sNeedle_cmd = r"/bin/bash -c 'needle -filter {0} {1} -outfile stdout -gapopen {2} -gapextend {3} -endweight Y -endopen {4} -endextend {5}'".format(sRef_seq, sQuery_seq, '20', '1', '20', '1')
        
    Needle_result = sp.Popen(sNeedle_cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True,shell=True)
    lResult = [Instance.seq._data for Instance in AlignIO.read(Needle_result.stdout, "emboss")]
    print(lResult)


def LoggingTest():
    import logging
    logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                    level=logging.DEBUG, filename='test.log', filemode='a'
                    )
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info('test')
    a = a * 10


def Main():
    #CountBar()
    #ExtractFastq()
    #LoadPickle()
    #CheckSearch()
    #CheckNeedle()
    LoggingTest()

Main()
