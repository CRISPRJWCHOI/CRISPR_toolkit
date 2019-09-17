import os, re, sys, pickle, logging

from time import time
import subprocess as sp
from pdb import set_trace

import numpy as np

from Bio import pairwise2
from CRISPResso2 import CRISPResso2Align


#EDNAFULL='/data/Indel_searcher_2/miniconda2/pkgs/crispresso2-2.0.30-py27h14c3975_0/lib/python2.7/site-packages/CRISPResso2/EDNAFULL'


"""
variable prefix

s: string
l: list
i: int
f: float
d: dictionary
"""

sManual = """
Usage:

python2.7 ./indel_search_ver1.0.py splitted_input_1.fq splitted_input_2.fq reference.fa

splitted_input_1.fq : forward
splitted_input_2.fq : reverse

Total FASTQ(fq) lines / 4 = remainder 0.
"""
#iQual_cutoff = 20

if len(sys.argv) > 1:
    sForward_fq_path = sys.argv[1]
    sReverse_fq_path = sys.argv[2]
    sRef_fa          = sys.argv[3]
    sPair            = sys.argv[4]
    sOG              = int(sys.argv[5])
    sOE              = int(sys.argv[6])
    iInsertion_win   = int(sys.argv[7])
    iDeletion_win    = int(sys.argv[8])
    sPAM_type        = sys.argv[9]       # Cpf1, Cas9
    sBarcode_PAM_pos = sys.argv[10]      # PAM - BARCODE type (reverse) or BARCODE - PAM type (forward)
    iQual_cutoff     = int(sys.argv[11])
    sOutput_dir      = sys.argv[12]
    sEDNAFULL        = sys.argv[13]
    sLogPath         = sys.argv[14]

else:
    print(sManual)
    sys.exit()

logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                    level=logging.DEBUG,
                    filename=sLogPath,
                    filemode='a')

logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

# index name, constant variable.
gNum_of_total = 0
gNum_of_ins   = 1
gNum_of_del   = 2
gNum_of_com   = 3
gTotal_FASTQ  = 4
gIns_FASTQ    = 5
gDel_FASTQ    = 6
gCom_FASTQ    = 7
gINDEL_info   = 8

aln_matrix = CRISPResso2Align.read_matrix(sEDNAFULL)

def Open_FASTQ_files():

    lFASTQ_forward = []
    lFASTQ_reverse = []
    lStore         = []

    # performance 20% up
    #lStore_append  = lStore.append
    #lFASTQ_forward_append = lFASTQ_forward.append
    #lFASTQ_reverse_append = lFASTQ_reverse.append

    # forward
    with open(sForward_fq_path) as fa_1:

        for i, sRow in enumerate(fa_1):

            i = i + 1
            sRow = sRow.replace('\n', '').upper()

            if i % 4 == 1 or i % 4 == 2:
                lStore.append(sRow)
            elif i % 4 == 0:
                lQual = [ord(i) - 33 for i in sRow]
                lStore.append(lQual)
                lFASTQ_forward.append(tuple(lStore))
                lStore = []

    # reverse
    dRev = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    if sPair == 'True':
        #with open('./6_AsD0_2_small_test.fq') as fa_2:
        with open(sReverse_fq_path) as fa_2:

            for i, sRow in enumerate(fa_2):
                i = i + 1
                sRow = sRow.replace('\n', '').upper()

                if i % 4 == 1:
                    lStore.append(sRow)
                elif i % 4 == 2:
                    lStore.append(''.join([dRev[sNucle] for sNucle in sRow[::-1]]))
                elif i % 4 == 0:
                    lQual = [ord(i) - 33 for i in sRow][::-1]
                    lStore.append(lQual)
                    lFASTQ_reverse.append(tuple(lStore))
                    lStore = []

    return lFASTQ_forward, lFASTQ_reverse
    #end1: return
#end: def

def Search_barcode_indel_postion():

    dRef = {}
    dResult = {}

    with open(sRef_fa) as Ref:

        iCount     = 0
        sBarcode   = ""
        sTarget_region = ""

        for sRow in Ref:
            iCount += 1

            if iCount % 2 != 0:
                #      barcode               target region
                # >CGCTCTACGTAGACA:CTCTATTACTCGCCCCACCTCCCCCAGCCC
                sBarcode_indel_seq = sRow.strip().replace('\n', '').replace('\r', '').split(':')
                sBarcode           = sBarcode_indel_seq[0].replace('>', '')
                sTarget_region     = sBarcode_indel_seq[1]

            elif iCount % 2 == 0:
                sRef_seq = sRow.strip().replace('\n', '').replace('\r', '')

                Seq_matcher = re.compile(r'(?=(%s))' % sTarget_region)

                #iIndel_start_pos       = sRef_seq.index(sTarget_region)               # There is possible to exist two indel.
                iIndel_start_pos        = Seq_matcher.finditer(sRef_seq)

                for i, match in enumerate(iIndel_start_pos):
                    iIndel_start_pos = match.start()
                #print iIndel_start_pos
                #print len(sTarget_region)
                #print sRef_seq
                iIndel_end_pos = iIndel_start_pos + len(sTarget_region) - 1

                try:
                    iBarcode_start_pos = sRef_seq.index(sBarcode)

                    if iIndel_start_pos <= iBarcode_start_pos:
                        raise IndexError('indel is before barcode')

                    iBarcode_end_pos        = iBarcode_start_pos + len(sBarcode) - 1
                    sRef_seq_after_barcode  = sRef_seq[iBarcode_end_pos + 1:]

                    # modified. to -1
                    iIndel_end_next_pos_from_barcode_end = iIndel_end_pos - iBarcode_end_pos -1

                    iIndel_start_next_pos_from_barcode_end = iIndel_start_pos - iBarcode_end_pos -1

                    #  "barcode"-------------*(N) that distance.
                    #          ^  ^            ^
                    #   *NNNN*NNNN
                    #    ^    ^     indel pos, the sequence matcher selects indel event pos front of it.

                    dRef[sBarcode] = (sRef_seq, sTarget_region, sRef_seq_after_barcode, iIndel_start_next_pos_from_barcode_end,
                                      iIndel_end_next_pos_from_barcode_end, iIndel_start_pos, iIndel_end_pos)  # total matched reads, insertion, deletion, complex
                    dResult[sBarcode] = [0, 0, 0, 0, [], [], [], [], []]
                except ValueError: continue

    assert len(dRef.keys()) == len(dResult.keys())
    return dRef, dResult
    #end1: return

#end: def

def Make_intermediate_file(lFASTQ_forward, lFASTQ_reverse, dRef, dResult):

    if not os.path.isdir('./dump'): os.mkdir('./dump')

    with open('./dump/lFASTQ_forward', 'wb') as lFASTQ_forward_dump,\
        open('./dump/lFASTQ_reverse', 'wb') as lFASTQ_reverse_dump, \
        open('./dump/dRef', 'wb') as dRef_dump,\
        open('./dump/Pre_dResult', 'wb') as dResult_dump:

        cPickle.dump(lFASTQ_forward, lFASTQ_forward_dump)
        cPickle.dump(lFASTQ_reverse, lFASTQ_reverse_dump)
        cPickle.dump(dRef, dRef_dump)
        cPickle.dump(dResult, dResult_dump)


def Search_indel(lFASTQ=[], dRef = {}, dResult={}, bRev=False):

    #lFASTQ : [(seq, qual),(seq, qual)]
    #lRef   : [(ref_seq, ref_seq_after_barcode, barcode, barcode end pos, indel end pos, indel from barcode),(...)]
    #dResult = [# of total, # of ins, # of del, # of com, [total FASTQ], [ins FASTQ], [del FASTQ], [com FASTQ]]
    iCount         = 0

    for lCol_FASTQ in lFASTQ:

        sName = lCol_FASTQ[0]
        sSeq  = lCol_FASTQ[1]
        lQual = lCol_FASTQ[2]

        assert isinstance(sName, str) and isinstance(sSeq, str) and isinstance(lQual, list)

        listSeqWindow = [sSeq[i:i + 26] for i in range(len(sSeq))[:-25]]

        iBarcode_matched = 0
        iNeedle_matched  = 0
        iInsert_count    = 0
        iDelete_count    = 0
        iComplex_count   = 0

        intFirstBarcode = 0  ## check whether a barcode is one in a sequence.

        for strSeqWindow in listSeqWindow:

            if intFirstBarcode == 1: break  ## A second barcode in a sequence is not considerable.

            try:
                lCol_ref = dRef[strSeqWindow]
                sBarcode = strSeqWindow
                intFirstBarcode = 1
            except KeyError:
                continue

            sRef_seq                      = lCol_ref[0]
            sTarget_region                = lCol_ref[1]
            iIndel_seq_len                = len(sTarget_region)
            sRef_seq_after_barcode        = lCol_ref[2]
            iIndel_start_from_barcode_pos = lCol_ref[3]
            iIndel_end_from_barcode_pos   = lCol_ref[4]
            i5bp_front_Indel_end          = iIndel_end_from_barcode_pos - 4    # NN(N)*NNN(N)*NNNN
            iIndel_start_pos              = lCol_ref[4]
            iIndel_end_pos                = lCol_ref[5]

            # Check the barcode pos and remove it.
            iBarcode_start_pos_FASTQ = sSeq.index(sBarcode)
            iBarcode_matched += 1
            print(iBarcode_matched)
            iBarcode_end_pos_FASTQ   = iBarcode_start_pos_FASTQ + len(sBarcode) - 1
            """
                junk seq  target region
            ref: AGGAG    AGAGAGAGAGA
            que: AGGAG    AGAGAGAGAGA
            But, It doesnt know where is the target region because of existed indels.
            So, There is no way not to include it.
            """
            # Use this.
            sQuery_seq_after_barcode  = sSeq[iBarcode_end_pos_FASTQ + 1:]
            lQuery_qual_after_barcode = lQual[iBarcode_end_pos_FASTQ:]

            sRef_seq    = r'<(echo -e ">{name}\n{seq}")'.format(name=sBarcode, seq=sRef_seq_after_barcode)
            sQuery_seq  = r'<(echo -e ">{name}\n{seq}")'.format(name=sName.split()[0].replace(':','_'), seq=sQuery_seq_after_barcode)

            sNeedle_cmd = r"/bin/bash -c 'needle -filter {0} {1} -outfile stdout -gapopen {2} -gapextend {3} -endweight Y -endopen {4} -endextend {5}'".format(sRef_seq, sQuery_seq, sOG, sOE, sEG, sEE)

            if iCount == 0:
               print(sNeedle_cmd)
               iCount = 1

            Needle_result = sp.Popen(sNeedle_cmd, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True,shell=True)
            lResult = [Instance.seq._data for Instance in AlignIO.read(Needle_result.stdout, "emboss")]
            sRef_needle_ori   = lResult[0]
            sQuery_needle_ori = lResult[1]

            Needle_result.stdout.close()

            # detach forward ---, backward ---
            # e.g.    ref   ------AAAGGCTACGATCTGCG------
            #         query AAAAAAAAATCGCTCTCGCTCTCCGATCT
            # trimmed ref         AAAGGCTACGATCTGCG
            # trimmed qeury       AAATCGCTCTCGCTCTC
            iReal_ref_needle_start = 0
            iReal_ref_needle_end   = len(sRef_needle_ori)
            iRef_needle_len        = len(sRef_needle_ori)

            for i, sRef_nucle in enumerate(sRef_needle_ori):
                if sRef_nucle in ['A', 'C', 'G', 'T']:
                    iReal_ref_needle_start = i
                    break

            for i, sRef_nucle in enumerate(sRef_needle_ori[::-1]):
                if sRef_nucle in ['A', 'C', 'G', 'T']:
                    iReal_ref_needle_end = iRef_needle_len - (i + 1)
                    # forward 0 1 2  len : 3
                    # reverse 2 1 0,  len - (2 + 1) = 0
                    break

            sRef_needle = sRef_needle_ori[iReal_ref_needle_start:iReal_ref_needle_end + 1]
            if iReal_ref_needle_start:
                sQuery_needle = sQuery_needle_ori[:iReal_ref_needle_end]
            sQuery_needle = sQuery_needle_ori[:len(sRef_needle)]
            # detaching completion

            # indel info making.
            iNeedle_match_pos_ref = 0
            iNeedle_match_pos_query = 0
            iNeedle_insertion = 0
            iNeedle_deletion = 0

            lInsertion_in_read = [] # insertion result [[100, 1], [119, 13]]
            lDeletion_in_read  = []  # deletion result  [[97, 1], [102, 3]]

            #print 'sRef_needle', sRef_needle
            #print 'sQuery_needle', sQuery_needle
            for i, (sRef_nucle, sQuery_nucle) in enumerate(zip(sRef_needle, sQuery_needle)):

                if sRef_nucle == '-':
                    iNeedle_insertion += 1

                if sQuery_nucle == '-':
                    iNeedle_deletion += 1

                if sRef_nucle in ['A', 'C', 'G', 'T']:
                    if iNeedle_insertion:
                        lInsertion_in_read.append([iNeedle_match_pos_ref, iNeedle_insertion])
                        iNeedle_insertion = 0
                    iNeedle_match_pos_ref += 1

                if sQuery_nucle in ['A', 'C', 'G', 'T']:
                    if iNeedle_deletion:
                        lDeletion_in_read.append([iNeedle_match_pos_query, iNeedle_deletion])
                        iNeedle_match_pos_query += iNeedle_deletion
                        iNeedle_deletion = 0
                    iNeedle_match_pos_query += 1
            #print 'sRef_needle', sRef_needle
            #print 'sQuery_needle', sQuery_needle
            #print 'lInsertion_in_read: onebase', lInsertion_in_read
            #print 'lDeletion_in_read: onebase', lDeletion_in_read
            #print 'i5bp_front_Indel_end', i5bp_front_Indel_end
            #print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos

            lTarget_indel_result = []  # ['20M2I', '23M3D' ...]
            """
            ins case
            ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNN*NNNN*AGCTT
            """
            for iMatch_pos, iInsertion_pos in lInsertion_in_read:
                #if i5bp_front_Indel_end == iMatch_pos -1 or iIndel_end_from_barcode_pos == iMatch_pos -1: # iMatch_pos is one base # original ver
                if i5bp_front_Indel_end - iInsertion_win <= iMatch_pos -1  <= i5bp_front_Indel_end + iInsertion_win or \
                iIndel_end_from_barcode_pos - iInsertion_win <= iMatch_pos -1 <= iIndel_end_from_barcode_pos + iInsertion_win: # iMatch_pos is one base
                    iInsert_count = 1
                    lTarget_indel_result.append(str(iMatch_pos)+'M'+str(iInsertion_pos)+'I')
            """
            del case 1
            ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNN(*)NN**AGCTT
            del case 2
            ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNN(*)NNNN**CTT
            """
            for iMatch_pos, iDeletion_pos in lDeletion_in_read:
                #print 'iMatch_pos, iDeletion_pos', iMatch_pos, iDeletion_pos
                #print 'i5bp_front_Indel_end, iIndel_end_from_barcode_pos', i5bp_front_Indel_end, iIndel_end_from_barcode_pos
                #print 'iMatch_pos -1, iMatch_pos + iDeletion_pos -1', iMatch_pos - 1, iMatch_pos + iDeletion_pos - 1

                """
                Insertion: 30M3I
                       ^
                ACGT---ACGT
                ACGTTTTACGT -> check this seq
                Insertion just check two position

                Deletion: 30M3D
                     ^
                ACGTTTTACGT
                ACGT---ACGT -> check this seq
                But deletion has to includes overlap deletion.
                """
                #if (i5bp_front_Indel_end <= iMatch_pos -1 <= iIndel_end_from_barcode_pos) or \
                #    (i5bp_front_Indel_end <= iMatch_pos + iDeletion_pos -1 <= iIndel_end_from_barcode_pos): # iMatch_pos is zero base
                if (iMatch_pos -iDeletion_win -1 <= i5bp_front_Indel_end and i5bp_front_Indel_end < (iMatch_pos + iDeletion_pos + iDeletion_win-1)) or \
                    (iMatch_pos -iDeletion_win -1 <= iIndel_end_from_barcode_pos and iIndel_end_from_barcode_pos < (iMatch_pos + iDeletion_pos + iDeletion_win-1)):
                    iDelete_count = 1
                    lTarget_indel_result.append(str(iMatch_pos)+'M'+str(iDeletion_pos)+'D')

            if iInsert_count == 1 and iDelete_count == 1:
                iComplex_count = 1
                iInsert_count  = 0
                iDelete_count  = 0

            #""" test set
            #print 'sBarcode', sBarcode
            #print 'sTarget_region', sTarget_region
            #print 'sRef_seq_after_barcode', sRef_seq_after_barcode
            #print 'sSeq_after_barcode', sQuery_seq
            #print 'iIndel_start_from_barcode_pos', iIndel_start_from_barcode_pos
            #print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos
            #"""

            lResult_FASTQ = [sName, sSeq, '+', ''.join(chr(i+33) for i in lQual)]
            dResult[sBarcode][gTotal_FASTQ].append(lResult_FASTQ)

            """
            iQual_end_pos + 1 is not correct, because the position is like this.
            *NNNN*(N)
            So, '+ 1' is removed.
            Howerver, seqeunce inspects until (N) position. indel is detected front of *(N).
            """

            # len(sQuery_seq_after_barcode) == len(lQuery_qual_after_barcode)
            if np.mean(lQuery_qual_after_barcode[iIndel_start_from_barcode_pos : iIndel_end_from_barcode_pos + 1]) >= iQual_cutoff:

                """
                23M3I
                23M is included junk_seq after barcode,

                barcorde  junk   targetseq   others
                *********ACCCT-------------ACACACACC
                so should select target region.
                If junk seq is removed by target region seq index pos.
                """
                # filter start,
                iTarget_start_from_barcode = sRef_seq_after_barcode.index(sTarget_region)
                sRef_seq_after_barcode     = sRef_seq_after_barcode[iTarget_start_from_barcode:]
                sQuery_seq_after_barcode   = sQuery_seq_after_barcode[iTarget_start_from_barcode:]
                lTrimmed_target_indel_result = []
                sContinue = 0

                for sINDEL in lTarget_indel_result:
                    # B - A is not included B position, so +1
                    iMatch_target_start = int(sINDEL.split('M')[0]) - iTarget_start_from_barcode
                    if iMatch_target_start < 0:
                        sContinue = 1
                    lTrimmed_target_indel_result.append(str(iMatch_target_start)+'M' + sINDEL.split('M')[1])
                # filter end

                if sContinue == 1:
                    continue
                #print 'Check'
                #print sRef_seq_after_barcode
                #print sQuery_seq_after_barcode
                #print lTrimmed_target_indel_result

                dResult[sBarcode][gINDEL_info].append([sRef_seq_after_barcode, sQuery_seq_after_barcode, lTrimmed_target_indel_result, sTarget_region, sRef_needle_ori, sQuery_needle_ori])
                if iInsert_count:
                    dResult[sBarcode][gIns_FASTQ].append(lResult_FASTQ)
                elif iDelete_count:
                    dResult[sBarcode][gDel_FASTQ].append(lResult_FASTQ)
                elif iComplex_count:
                    dResult[sBarcode][gCom_FASTQ].append(lResult_FASTQ)
            else:
                iInsert_count  = 0
                iDelete_count  = 0
                iComplex_count = 0

            # total matched reads, insertion, deletion, complex
            dResult[sBarcode][gNum_of_total] += iBarcode_matched
            dResult[sBarcode][gNum_of_ins]   += iInsert_count
            dResult[sBarcode][gNum_of_del]   += iDelete_count
            dResult[sBarcode][gNum_of_com]   += iComplex_count

            iBarcode_matched = 0
            iInsert_count = 0
            iDelete_count = 0
            iComplex_count = 0
        #end: for lCol_FASTQ
    #END:for
    return dResult


def Calculate_INDEL_frequency(dResult):
    dResult_INDEL_freq = {}

    for sBarcode, lValue in dResult.items():  #lValue[gINDEL_info] : [[sRef_seq_after_barcode, sQuery_seq_after_barcode, lTarget_indel_result, sTarget_region], ..])
        sRef_seq_loop = ''
        llINDEL_store = []    # ['ACAGACAGA', ['20M2I', '23M3D']]
        dINDEL_freq = {}

        if lValue[gINDEL_info]:
            for sRef_seq_loop, sQuery_seq, lINDEL, sTarget_region, sRef_needle, sQuery_needle in lValue[gINDEL_info]: # llINDEL : [['20M2I', '23M3D'], ...]
                #print 'lINDEL', lINDEL
                for sINDEL in lINDEL:
                    llINDEL_store.append([sQuery_seq, sINDEL, sRef_needle, sQuery_needle])

            iTotal = len([lINDEL for sQuery_seq, lINDEL, sRef_needle, sQuery_needle in llINDEL_store])

            for sQuery_seq, sINDEL, sRef_needle, sQuery_needle in llINDEL_store:
                dINDEL_freq[sINDEL] = [[], 0, [], []]

            for sQuery_seq, sINDEL, sRef_needle, sQuery_needle in llINDEL_store:
                dINDEL_freq[sINDEL][1] += 1
                dINDEL_freq[sINDEL][0].append(sQuery_seq)
                dINDEL_freq[sINDEL][2].append(sRef_needle)
                dINDEL_freq[sINDEL][3].append(sQuery_needle)

            for sINDEL in dINDEL_freq:
                lQuery        = dINDEL_freq[sINDEL][0]
                iFreq         = dINDEL_freq[sINDEL][1]
                lRef_needle   = dINDEL_freq[sINDEL][2]
                lQuery_needle = dINDEL_freq[sINDEL][3]

                try:
                    dResult_INDEL_freq[sBarcode].append([sRef_seq_loop, lQuery, sINDEL, float(iFreq)/iTotal,
                                                         sTarget_region, lRef_needle, lQuery_needle])
                except (KeyError, TypeError, AttributeError) as e:
                    dResult_INDEL_freq[sBarcode] = []
                    dResult_INDEL_freq[sBarcode].append([sRef_seq_loop, lQuery, sINDEL, float(iFreq)/iTotal,
                                                         sTarget_region, lRef_needle, lQuery_needle])
        #end: if lValue[gINDEL_info]
    #end: for sBarcode, lValue
    return dResult_INDEL_freq
    #end1: return
#end: def


class CAS9_parser():

    def __init__(self):
        return

    def Search_barcode_indel_position(self, sBarcode_PAM_pos):

        dRef = {}
        dResult = {}

        with open(sRef_fa) as Ref:

            iCount         = 0
            sBarcode       = ""
            sTarget_region = ""
            intBarcodeLen  = 0

            for sRow in Ref:
                iCount += 1

                if iCount % 2 != 0:
                    #      barcode               target region
                    # >CGCTCTACGTAGACA:CTCTATTACTCGCCCCACCTCCCCCAGCCC
                    sBarcode_indel_seq = sRow.strip().replace('\n', '').replace('\r', '').split(':')
                    sBarcode = sBarcode_indel_seq[0].replace('>', '')

                    if intBarcodeLen > 0:
                        assert intBarcodeLen == len(sBarcode), 'All of the barcode lengths must be same.'
                    intBarcodeLen = len(sBarcode)

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

        assert len(dRef.keys()) == len(dResult.keys())

        return dRef, dResult
        # end1: return

    def Search_barcode_indel_position_forward(self):
        return


    def Search_indel(self, lFASTQ=[], dRef = {}, dResult={}, sBarcode_PAM_pos=""):

        # lFASTQ : [(seq, qual),(seq, qual)]
        # lRef   : [(ref_seq, ref_seq_after_barcode, barcode, barcode end pos, indel end pos, indel from barcode),(...)]
        # dResult = [# of total, # of ins, # of del, # of com, [total FASTQ], [ins FASTQ], [del FASTQ], [com FASTQ]]
        iCount = 0

        intBarcodeLen = len(dRef.keys()[0])
        #print('intBarcodeLen', intBarcodeLen)

        for lCol_FASTQ in lFASTQ:
            sName = lCol_FASTQ[0]
            if sBarcode_PAM_pos == 'Reverse':
                sSeq  = lCol_FASTQ[1][::-1]
                lQual = lCol_FASTQ[2][::-1]
            else:
                sSeq  = lCol_FASTQ[1]
                lQual = lCol_FASTQ[2]

            assert isinstance(sName, str) and isinstance(sSeq, str) and isinstance(lQual, list)

            listSeqWindow = [sSeq[i:i + intBarcodeLen] for i in range(len(sSeq))[:-intBarcodeLen - 1]]

            iBarcode_matched = 0
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

                sRef_seq = lCol_ref[0]
                sTarget_region = lCol_ref[1]
                iIndel_seq_len = len(sTarget_region)
                sRef_seq_after_barcode = lCol_ref[2]
                iIndel_start_from_barcode_pos = lCol_ref[3]
                iIndel_end_from_barcode_pos = lCol_ref[4]
                try:
                    i6bp_front_Indel_end = iIndel_end_from_barcode_pos - 6  # NN(N)*NNN(N)*NNNN
                except Exception:
                    set_trace()

                """
                                                     *     ^ : iIndel_end_from_barcode_pos
                                  GGCG   TCGCTCATGTACCTCCCGT
                TATAGTCTGTCATGCGATGGCG---TCGCTCATGTACCTCCCGTTACAGCCACAAAGCAGGA
                     *
                GGCGTC GCTCATGTACCTCCCGT
                  6          17 
                """

                iIndel_start_pos = lCol_ref[4]
                iIndel_end_pos = lCol_ref[5]

                ## bug fix
                if sBarcode == "": continue

                # Check the barcode pos and remove it.
                iBarcode_start_pos_FASTQ = sSeq.index(sBarcode)
                iBarcode_matched += 1
                iBarcode_end_pos_FASTQ = iBarcode_start_pos_FASTQ + len(sBarcode) - 1
                """
                    junk seq  target region
                ref: AGGAG    AGAGAGAGAGA
                que: AGGAG    AGAGAGAGAGA
                But, It doesnt know where is the target region because of existed indels.
                So, There is no way not to include it.
                """
                # Use this.
                sQuery_seq_after_barcode = sSeq[iBarcode_end_pos_FASTQ + 1:]
                lQuery_qual_after_barcode = lQual[iBarcode_end_pos_FASTQ:]

                #if iCount == 0: ## For printing the cmd
                #    print(sNeedle_cmd)
                #    iCount = 1

                #""" cripsress no incentive == gotoh
                amp_len = len(sRef_seq_after_barcode)
                gap_incentive = np.zeros(amp_len + 1, dtype=np.int)

                try:
                    lResult = CRISPResso2Align.global_align(sQuery_seq_after_barcode.upper(), sRef_seq_after_barcode.upper(),
                                                            matrix=aln_matrix, gap_open=sOG, gap_extend=sOE, gap_incentive=gap_incentive)
                except Exception as e:
                    logging.error(e, exc_info=True)
                    continue

                sQuery_needle_ori = lResult[0]
                sRef_needle_ori   = lResult[1]

                # detach forward ---, backward ---
                # e.g.    ref   ------AAAGGCTACGATCTGCG------
                #         query AAAAAAAAATCGCTCTCGCTCTCCGATCT
                # trimmed ref         AAAGGCTACGATCTGCG
                # trimmed qeury       AAATCGCTCTCGCTCTC
                iReal_ref_needle_start = 0
                iReal_ref_needle_end = len(sRef_needle_ori)
                iRef_needle_len = len(sRef_needle_ori)

                for i, sRef_nucle in enumerate(sRef_needle_ori):
                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        iReal_ref_needle_start = i
                        break

                for i, sRef_nucle in enumerate(sRef_needle_ori[::-1]):
                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        iReal_ref_needle_end = iRef_needle_len - (i + 1)
                        # forward 0 1 2  len : 3
                        # reverse 2 1 0,  len - (2 + 1) = 0
                        break

                sRef_needle = sRef_needle_ori[iReal_ref_needle_start:iReal_ref_needle_end + 1]
                if iReal_ref_needle_start:
                    sQuery_needle = sQuery_needle_ori[:iReal_ref_needle_end]
                sQuery_needle = sQuery_needle_ori[:len(sRef_needle)]
                # detaching completion

                # indel info making.
                iNeedle_match_pos_ref = 0
                iNeedle_match_pos_query = 0
                iNeedle_insertion = 0
                iNeedle_deletion = 0

                lInsertion_in_read = []  # insertion result [[100, 1], [119, 13]]
                lDeletion_in_read = []  # deletion result  [[97, 1], [102, 3]]

                # print 'sRef_needle', sRef_needle
                # print 'sQuery_needle', sQuery_needle
                for i, (sRef_nucle, sQuery_nucle) in enumerate(zip(sRef_needle, sQuery_needle)):

                    if sRef_nucle == '-':
                        iNeedle_insertion += 1

                    if sQuery_nucle == '-':
                        iNeedle_deletion += 1

                    if sRef_nucle in ['A', 'C', 'G', 'T']:
                        if iNeedle_insertion:
                            lInsertion_in_read.append([iNeedle_match_pos_ref, iNeedle_insertion])
                            iNeedle_insertion = 0
                        iNeedle_match_pos_ref += 1

                    if sQuery_nucle in ['A', 'C', 'G', 'T']:
                        if iNeedle_deletion:
                            lDeletion_in_read.append([iNeedle_match_pos_query, iNeedle_deletion])
                            iNeedle_match_pos_query += iNeedle_deletion
                            iNeedle_deletion = 0
                        iNeedle_match_pos_query += 1
                        # print 'sRef_needle', sRef_needle

                # print 'sQuery_needle', sQuery_needle
                # print 'lInsertion_in_read: onebase', lInsertion_in_read
                # print 'lDeletion_in_read: onebase', lDeletion_in_read
                # print 'i5bp_front_Indel_end', i5bp_front_Indel_end
                # print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos

                lTarget_indel_result = []  # ['20M2I', '23M3D' ...]

                """
                ins case
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNN*NNNNNAGCTT
                """
                for iMatch_pos, iInsertion_pos in lInsertion_in_read:
                    # if i5bp_front_Indel_end == iMatch_pos -1 or iIndel_end_from_barcode_pos == iMatch_pos -1: # iMatch_pos is one base # original ver
                    if i6bp_front_Indel_end - iInsertion_win <= iMatch_pos - 1 <= i6bp_front_Indel_end + iInsertion_win:  # iMatch_pos is one base
                        iInsert_count = 1
                        lTarget_indel_result.append(str(iMatch_pos) + 'M' + str(iInsertion_pos) + 'I')
                """
                del case 1
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNAGCTT
                del case 2
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNNNCTT
                """
                for iMatch_pos, iDeletion_pos in lDeletion_in_read:
                    # print 'iMatch_pos, iDeletion_pos', iMatch_pos, iDeletion_pos
                    # print 'i5bp_front_Indel_end, iIndel_end_from_barcode_pos', i5bp_front_Indel_end, iIndel_end_from_barcode_pos
                    # print 'iMatch_pos -1, iMatch_pos + iDeletion_pos -1', iMatch_pos - 1, iMatch_pos + iDeletion_pos - 1
                    """
                    Insertion: 30M3I
                           ^
                    ACGT---ACGT
                    ACGTTTTACGT -> check this seq
                    Insertion just check two position

                    Deletion: 30M3D
                         ^
                    ACGTTTTACGT
                    ACGT---ACGT -> check this seq
                    But deletion has to includes overlap deletion.
                    """
                    # if (i5bp_front_Indel_end <= iMatch_pos -1 <= iIndel_end_from_barcode_pos) or \
                    #    (i5bp_front_Indel_end <= iMatch_pos + iDeletion_pos -1 <= iIndel_end_from_barcode_pos): # iMatch_pos is zero base
                    if (iMatch_pos - iDeletion_win - 1 <= i6bp_front_Indel_end and i6bp_front_Indel_end < (
                                iMatch_pos + iDeletion_pos + iDeletion_win - 1)):
                        iDelete_count = 1
                        lTarget_indel_result.append(str(iMatch_pos) + 'M' + str(iDeletion_pos) + 'D')

                if iInsert_count == 1 and iDelete_count == 1:
                    iComplex_count = 1
                    iInsert_count = 0
                    iDelete_count = 0
                # """ test set
                # print 'sBarcode', sBarcode
                # print 'sTarget_region', sTarget_region
                # print 'sRef_seq_after_barcode', sRef_seq_after_barcode
                # print 'sSeq_after_barcode', sQuery_seq
                # print 'iIndel_start_from_barcode_pos', iIndel_start_from_barcode_pos
                # print 'iIndel_end_from_barcode_pos', iIndel_end_from_barcode_pos
                # """
                lResult_FASTQ = [sName, sSeq, '+', ''.join(chr(i + 33) for i in lQual)]
                dResult[sBarcode][gTotal_FASTQ].append(lResult_FASTQ)
                """
                iQual_end_pos + 1 is not correct, because the position is like this.
                *NNNN*(N)
                So, '+ 1' is removed.
                Howerver, seqeunce inspects until (N) position. indel is detected front of *(N).
                """
                ################################################################
                #print(lTarget_indel_result)
                #set_trace()
                # len(sQuery_seq_after_barcode) == len(lQuery_qual_after_barcode)
                if np.mean(lQuery_qual_after_barcode[
                           iIndel_start_from_barcode_pos: iIndel_end_from_barcode_pos + 1]) >= iQual_cutoff:
                    """
                    23M3I
                    23M is included junk_seq after barcode,

                    barcorde  junk   targetseq   others
                    *********ACCCT-------------ACACACACC
                    so should select target region.
                    If junk seq is removed by target region seq index pos.
                    """
                    # filter start,
                    iTarget_start_from_barcode = sRef_seq_after_barcode.index(sTarget_region)
                    sRef_seq_after_barcode = sRef_seq_after_barcode[iTarget_start_from_barcode:]
                    sQuery_seq_after_barcode = sQuery_seq_after_barcode[iTarget_start_from_barcode:]
                    lTrimmed_target_indel_result = []
                    sContinue = 0
                    for sINDEL in lTarget_indel_result:
                        # B - A is not included B position, so +1
                        iMatch_target_start = int(sINDEL.split('M')[0]) - iTarget_start_from_barcode
                        """ This part determines a deletion range.
                                                  ^ current match pos                                           
                        AGCTACGATCAGCATCTGACTTACTTC[barcode]
                        

                                       ^ fix the match start at here. (target region)                                           
                        AGCTACGATCAGCATC TGACTTACTTC[barcode]
                        
                        if iMatch_target_start < 0:
                            sContinue = 1
                        
                        But, this method has some problems.
                        
                                       ^ barcode start
                        AGCTACGATCAGCAT*********C[barcode]
                        Like this pattern doesn't seleted. because, deletion checking is begun the target region start position. 
                        Thus, I have fixed this problem.
                        """

                        if iMatch_target_start <= -(iTarget_start_from_barcode):
                            #print(iMatch_target_start, iTarget_start_from_barcode)
                            continue
                            #sContinue = 1

                        lTrimmed_target_indel_result.append(str(iMatch_target_start) + 'M' + sINDEL.split('M')[1])
                    # filter end

                    if sContinue == 1:
                        continue
                    # print 'Check'
                    # print sRef_seq_after_barcode
                    # print sQuery_seq_after_barcode
                    # print lTrimmed_target_indel_result
                    # print('Trimmed', lTrimmed_target_indel_result)

                    dResult[sBarcode][gINDEL_info].append(
                        [sRef_seq_after_barcode, sQuery_seq_after_barcode, lTrimmed_target_indel_result,
                         sTarget_region, sRef_needle_ori, sQuery_needle_ori])
                    if iInsert_count:
                        dResult[sBarcode][gIns_FASTQ].append(lResult_FASTQ)
                    elif iDelete_count:
                        dResult[sBarcode][gDel_FASTQ].append(lResult_FASTQ)
                    elif iComplex_count:
                        dResult[sBarcode][gCom_FASTQ].append(lResult_FASTQ)
                else:
                    iInsert_count  = 0
                    iDelete_count  = 0
                    iComplex_count = 0

                # total matched reads, insertion, deletion, complex
                dResult[sBarcode][gNum_of_total] += iBarcode_matched
                dResult[sBarcode][gNum_of_ins] += iInsert_count
                dResult[sBarcode][gNum_of_del] += iDelete_count
                dResult[sBarcode][gNum_of_com] += iComplex_count

                iBarcode_matched = 0
                iInsert_count    = 0
                iDelete_count    = 0
                iComplex_count   = 0

            #End:for
        #END:for
        return dResult

    def Calculate_INDEL_frequency(self):
        return

#END:class


def MakePickleOutput(dictResult, dictResultIndelFreq, strBarcodePamPos=''):

    dictOutput = {'dictResult' : dictResult,
                  'dictResultIndelFreq' : dictResultIndelFreq,
                  'strBarcodePamPos' : strBarcodePamPos}

    with open('{outdir}/Tmp/Pickle/{fq}.pickle'.format(outdir=sOutput_dir, fq=os.path.basename(sForward_fq_path)), 'wb') as Pickle:
        pickle.dump(dictOutput, Pickle)


def Main():
    
    # Output: 1. Count information of matched barcode e.g. TACGATCTA\t# total\tins\t# del\t# com
    # Output: 2. classify FASTQ.    e.g. TAGAATATACACG.insertion.fastq
    #if not os.path.isdir('{outdir}'.format(outdir=sOutput_dir)): os.mkdir('{outdir}'.format(outdir=sOutput_dir))
    #if not os.path.isdir('{outdir}/Summary/'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/Summary/'.format(outdir=sOutput_dir))
    #if not os.path.isdir('{outdir}/Result/'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/Result/'.format(outdir=sOutput_dir))
    #if not os.path.isdir('{outdir}/Pickle/'.format(outdir=sOutput_dir)): os.mkdir('{outdir}/Pickle/'.format(outdir=sOutput_dir))

    logging.info('Program start : %s' % sForward_fq_path)

    logging.info('File Open')
    lFASTQ_forward, lFASTQ_reverse = Open_FASTQ_files()

    if sPAM_type == 'Cpf1':
        logging.info('Search barcode INDEL pos')
        dRef, dResult = Search_barcode_indel_postion()  # ref check.

        logging.info('Search INDEL forward')
        dResult_forward = Search_indel(lFASTQ_forward, dRef, dResult)

        if sPair == 'True':
            logging.info('Search INDEL reverse')
            dResult_reverse = Search_indel(lFASTQ_reverse, dRef, dResult_forward)

            logging.info('Calculate INDEL frequency')
            dResult_INDEL_freq = Calculate_INDEL_frequency(dResult_reverse)

            logging.info('Make pickle output reverse')
            MakePickleOutput(dResult_reverse, dResult_INDEL_freq)

        else:
            logging.info('Calculate INDEL frequency')
            dResult_INDEL_freq = Calculate_INDEL_frequency(dResult_forward)

            logging.info('Make pickle output forward')
            MakePickleOutput(dResult_forward, dResult_INDEL_freq)

    elif sPAM_type == 'Cas9':
        Ins_CAS9        = CAS9_parser()
        logging.info('Search barcode INDEL pos')
        dRef, dResult   = Ins_CAS9.Search_barcode_indel_position(sBarcode_PAM_pos)
        logging.info('Search INDEL')
        dResult_forward = Ins_CAS9.Search_indel(lFASTQ_forward, dRef, dResult, sBarcode_PAM_pos)
        logging.info('Calculate INDEL frequency')
        dResult_INDEL_freq = Calculate_INDEL_frequency(dResult_forward)

        logging.info('Make pickle output forward')
        MakePickleOutput(dResult_forward, dResult_INDEL_freq, sBarcode_PAM_pos)

    logging.info('Program end : %s' % sForward_fq_path)
#END:def


if __name__ == '__main__':
    Main()


