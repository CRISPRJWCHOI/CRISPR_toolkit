import os, re, sys, logging

import numpy as np
import subprocess as sp
import cPickle as pickle

from pdb import set_trace
from datetime import datetime
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import CoreGotoh


class clsParameter(object):

    """
    ./BaseEdit_freq_crispresso.py {forw} {GapO} {GapE} {barcode} {ref} {target_window} {indel_check_pos}
     {target_ref_alt} {outdir} {file_name} {PAM_seq} {PAM_pos} {Guide_pos} {ednafull} {log}
    """
    def __init__(self):

        if len(sys.argv) > 1:
            self.strForwPath        = sys.argv[1]
            self.floOg              = float(sys.argv[2])
            self.floOe              = float(sys.argv[3])
            self.strBarcode         = sys.argv[4]
            strRef                  = sys.argv[5]
            self.strRef             = strRef[strRef.index(self.strBarcode):]  ## 'ACTG'<barcode>ACGACACACGCAT, leftside bases are redundant.
            self.listTargetWindow   = sys.argv[6].split('-')
            self.listIndelCheckPos  = sys.argv[7].split('-')
            self.listTargetRefAlt   = sys.argv[8].split(',')
            self.strOutputDir       = sys.argv[9]
            self.strFileName        = sys.argv[10]
            self.strPamSeq          = sys.argv[11]
            self.listPamPos         = sys.argv[12].split('-')
            self.listGuidePos       = sys.argv[13].split('-')
            self.strEDNAFULL        = os.path.abspath('../EDNAFULL')
            self.strLogPath         = sys.argv[14]

        else:
            sManual = """
            Usage:

            python2.7 ./indel_search_ver1.0.py splitted_input_1.fq splitted_input_2.fq reference.fa

            splitted_input_1.fq : forward
            splitted_input_2.fq : reverse

            Total FASTQ(fq) lines / 4 = remainder 0.
            """
            print sManual
            sys.exit()


class clsBaseEditParser():

    def __init__(self, InstParameter):
        self.strForwPath       = InstParameter.strForwPath
        self.strRef            = InstParameter.strRef
        self.strBarcode        = InstParameter.strBarcode
        self.strEDNAFULL       = InstParameter.strEDNAFULL
        self.floOg             = InstParameter.floOg
        self.floOe             = InstParameter.floOe
        self.listIndelCheckPos = InstParameter.listIndelCheckPos
        self.listTargetWindow  = InstParameter.listTargetWindow

    def OpenSequenceFiles(self):
        lSequence_forward = []
        with open(self.strForwPath) as fa_1:
            lSequence_forward = [sRow.replace('\n', '').upper() for sRow in fa_1]
        return lSequence_forward

    def CalculateBaseEditFreq(self, lQuery_seq=[]):

        dRef    = {}
        dResult = {}

        dRef[self.strBarcode]    = (self.strRef)  # total matched reads, insertion, deletion, complex
        dResult[self.strBarcode] = [0, 0, 0, 0, [], [], [], [], [], [], []]

        # lRef   : [(ref_seq, ref_seq_after_barcode, barcode, barcode end pos, indel end pos, indel from barcode),(...)]
        # dResult = [# of total, # of ins, # of del, # of com, [total FASTQ], [ins FASTQ], [del FASTQ], [com FASTQ], info]
        iCount = 0

        InstGotoh = CoreGotoh(strEDNAFULL=self.strEDNAFULL, floOg=self.floOg, floOe=self.floOe)

        for sQuery_seq_raw in lQuery_seq:

            iBarcode_matched = 0
            iNeedle_matched  = 0
            iInsert_count    = 0
            iDelete_count    = 0
            iComplex_count   = 0

            try:
                # Check the barcode pos and remove it.
                sQuery_seq_raw = sQuery_seq_raw.replace('\r', '')
                iBarcode_start_pos = sQuery_seq_raw.index(self.strBarcode)
                iBarcode_matched += 1

                sQuery_seq_with_barcode = sQuery_seq_raw[iBarcode_start_pos:]  ## this is not after barcode seq. including barcode

                npGapIncentive = InstGotoh.GapIncentive(self.strRef)

                try:
                    lResult = InstGotoh.RunCRISPResso2(sQuery_seq_with_barcode.upper(), self.strRef.upper(), npGapIncentive)
                except Exception as e:
                    logging.error(e, exc_info=True)
                    continue

                sQuery_needle_ori = lResult[0]
                sRef_needle_ori   = lResult[1]

                # if _check == 1:
                #     print(sRef_needle_ori)
                #     print(sQuery_needle_ori)
                #     set_trace()

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
                iNeedle_match_pos_ref   = 0
                iNeedle_match_pos_query = 0
                iNeedle_insertion       = 0
                iNeedle_deletion        = 0

                lInsertion_in_read = []  # insertion result [[100, 1], [119, 13]]
                lDeletion_in_read  = []  # deletion result  [[97, 1], [102, 3]]

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

                iCleavage_window_start = int(self.listIndelCheckPos[0])
                iCleavage_window_end = int(self.listIndelCheckPos[1]) - 1

                for iMatch_pos, iInsertion_pos in lInsertion_in_read:
                    if iCleavage_window_start <= iMatch_pos <= iCleavage_window_end:  # iMatch_pos is one base
                        iInsert_count = 1
                        lTarget_indel_result.append(str(iMatch_pos) + 'M' + str(iInsertion_pos) + 'I')
                """
                del case 1
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNAGCTT
                del case 2
                ...............................NNNNNNNNNNNNNN....NNNNNNNNNNNNNNNNNNNNN**NNNNNCTT
                """
                for iMatch_pos, iDeletion_pos in lDeletion_in_read:

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

                    if iMatch_pos <= iCleavage_window_end and iCleavage_window_start <= (iMatch_pos + iDeletion_pos):
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

                    """
                    23M3I
                    23M is included junk_seq after barcode,

                    barcorde  junk   targetseq   others
                    *********ACCCT-------------ACACACACC
                    so should select target region.
                    If junk seq is removed by target region seq index pos.
                    """

                ## 8: indel info
                dResult[self.strBarcode][8].append(
                    [self.strRef, sQuery_seq_raw, lTarget_indel_result,
                     "", sRef_needle_ori, sQuery_needle_ori])  ## "" -> target seq, but this is not used this project.

            # end: try
            except ValueError as e:
                print(e)
                continue

            # total matched reads, insertion, deletion, complex
            dResult[self.strBarcode][0] += iBarcode_matched
            dResult[self.strBarcode][1] += iInsert_count
            dResult[self.strBarcode][2] += iDelete_count
            dResult[self.strBarcode][3] += iComplex_count

            ## base editing frequency
            """
                   BaseEditPos : 0                                                    1                                  2
            [OrderedDict([('A',0),('C',0),('G',0),('T',0)]), OrderedDict([('A',0),('C',0),('G',0),('T',0)]), ...

            and sum the counts each position
            """

            ## No indel reads only
            if iInsert_count == 0 and iDelete_count == 0 and iComplex_count == 0:

                lBaseEdit = []
                iTarget_len = int(self.listTargetWindow[1]) - int(self.listTargetWindow[0]) + 1

                for i in range(iTarget_len):
                    lBaseEdit.append(OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]))

                iTarget_start = int(self.listTargetWindow[0]) - 1
                iTarget_end = int(self.listTargetWindow[1])

                """
                                       cleavage window start
                                        ^
                [barcode]ACGACGTACGACGT[cleavage]
                [barcode]ACGACGTACGACGT[cleavage]
                """

                iBase_edit_event = 0

                for i, tRef_Query_base in enumerate(zip(sRef_needle[iTarget_start: iTarget_end], sQuery_needle[iTarget_start: iTarget_end])):
                    sRef_base   = tRef_Query_base[0]
                    sQuery_base = tRef_Query_base[1]

                    if sRef_base == '-' or sQuery_base == '-': continue

                    if sRef_base != sQuery_base and sQuery_base != 'N':
                        iBase_edit_event = 1
                        lBaseEdit[i][sQuery_base] += 1
                        # print(sQuery_needle)

                dResult[self.strBarcode][9].append(lBaseEdit)
                ## Processed indel filtering and store aligned alt mut read.
                if iBase_edit_event == 1:
                    dResult[self.strBarcode][10].append([self.strRef, sQuery_seq_raw, lTarget_indel_result, [list(orderedDict.values()) for orderedDict in lBaseEdit], sRef_needle_ori, sQuery_needle_ori])
                # dResult[sBarcode] = [0, 0, 0, 0, [], [], [], [], [], [BaseEdit_freq_data]]

            iBarcode_matched = 0
            iInsert_count = 0
            iDelete_count = 0
            iComplex_count = 0
            # end: for sBarcode, lCol_ref
        # end: for lCol_FASTQ
        return dResult


class clsOutputMaker():

    def __init__(self, InstParameter):

        self.strForwPath      = InstParameter.strForwPath
        self.strRef           = InstParameter.strRef
        self.strFileName      = InstParameter.strFileName
        self.strOutputDir     = InstParameter.strOutputDir
        self.listTargetRefAlt = InstParameter.listTargetRefAlt
        self.listTargetWindow = InstParameter.listTargetWindow
        self.strPamSeq        = InstParameter.strPamSeq
        self.listPamPos       = InstParameter.listPamPos
        self.listGuidePos     = InstParameter.listGuidePos

        # index name, constant variable.
        self.intNumOfTotal = 0
        self.intNumOfIns   = 1
        self.intNumOfDel   = 2
        self.intNumOfCom   = 3
        self.intTotalFastq = 4
        self.intInsFastq   = 5
        self.intDelFastq   = 6
        self.intComFastq   = 7
        self.intIndelInfo  = 8

    def MakeOutput(self, dResult):
        """
       {'TTTGGTGCACACACATATA': [6, 2, 2, 0, [], [], [], [], [['TATCTCTA..ref', 'GAGTCGGTG...query', [13M5D], '',
       'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG------------------------------------------------------------',
       'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGA-TGCGGGGCGTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA'],
       ['TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG', '', '', '',
       'TTTGGTGCACACACATATAACTGGAACACAAAGCATAGACTGCGGGGCG------------------------------------------------------------', ...
       [[OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 1)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)])],
       [OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 1)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)]), OrderedDict([('A', 0), ('C', 0), ('G', 0), ('T', 0)])]]]}
        """

        with open('{outdir}/Tmp/Alignment/{file_name}_filtered_indel.txt'.format(outdir=self.strOutputDir, file_name=self.strFileName), 'w') as Filtered,\
            open('{outdir}/Tmp/Alignment/{file_name}_aligned_BaseEdit.txt'.format(outdir=self.strOutputDir, file_name=self.strFileName), 'w') as Ref_Alt_edit:

            for sBarcode in dResult:
                for lAligned_indel_result in dResult[sBarcode][8]:  # 8 : indel list
                    if lAligned_indel_result[2]:
                        Filtered.write('\t'.join(map(str, lAligned_indel_result)) + '\n')

                for lAligned_alt_result in dResult[sBarcode][10]:  # 10 : alt base list
                    if lAligned_alt_result:
                        lAligned_alt_result[2] = str(lAligned_alt_result[2])
                        try:
                            Ref_Alt_edit.write('\t'.join(map(str, lAligned_alt_result)) + '\n')
                        except Exception:
                            set_trace()

            """
            lAligned_result
            ['TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC', 'TCTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCAAAGGTCCAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAAT\r',
            ['38M1D'], '', 'TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC-----------------------------------------------------------------', 'TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGT-CAAAGGTCCAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAAT']
            """

        dSelect_base = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        sTarget_ref = self.listTargetRefAlt[0]
        sTarget_alt = self.listTargetRefAlt[1]

        iTarget_base = dSelect_base[sTarget_alt]

        try:
            if not os.path.isdir('{outdir}/Tmp/All'.format(outdir=self.strOutputDir)):
                os.mkdir('{outdir}/Tmp/All'.format(outdir=self.strOutputDir))
            if not os.path.isdir('{outdir}/Tmp/Target'.format(outdir=self.strOutputDir)):
                os.mkdir('{outdir}/Tmp/Target'.format(outdir=self.strOutputDir))
        except OSError:
            pass

        for sBarcode, lValue in dResult.items():

            iBarcode_start_pos       = self.strRef.index(sBarcode)
            sRef_seq_without_barcode = self.strRef[iBarcode_start_pos+len(sBarcode):]

            llBaseEdit = lValue[9]
            lSum = []

            for i, lBaseEdit in enumerate(llBaseEdit):

                if not lSum:
                    lSum = [[0, 0, 0, 0] for iQuery in range(len(lBaseEdit))]

                for j in range(len(lBaseEdit)):
                    for k, iCount in enumerate(list(llBaseEdit[i][j].values())):
                        lSum[j][k] += iCount

            with open('{outdir}/Tmp/All/{file_name}_Summary.txt'.format(outdir=self.strOutputDir, file_name=self.strFileName), 'w') as Summary, \
                open('{outdir}/Tmp/Target/{file_name}_{target}_Summary.txt'.format(outdir=self.strOutputDir, file_name=self.strFileName, target=sTarget_ref + 'to' + sTarget_alt), 'w') as Target_summary:

                ## This Ref has barcode.
                sRef_target = self.strRef[int(self.listTargetWindow[0]) - 1:int(self.listTargetWindow[1])]

                iPAM_start    = int(self.listPamPos[0]) - 1
                iPAM_end      = int(self.listPamPos[1])
                iGuide_start  = int(self.listGuidePos[0]) - 1
                iGuide_end    = int(self.listGuidePos[1])
                iGuide_len    = iGuide_end - iGuide_start
                iBarcode_len  = len(sBarcode)

                """
                barcode Guide st,ed 
                <----><----------> NGG
                ACGTACGTACGTACGTACGTGGACG
                """

                #sRef_target[iPAM_start:iPAM_end] = sPAM_seq
                ## iWithout_target_len = len(sRef_target[iBarcode_len:iGuide_start]) -> weird part.
                ## So I corrected it.
                iWithout_target_len = iGuide_start - iBarcode_len
                lWithout_target_pos = [-(i+1) for i in range(iWithout_target_len)][::-1]

                lWith_target_pos = [i + 1 for i in range(iGuide_len)]
                lAfter_PAM_pos   = [i + 1 for i in range(len(self.strRef) - iPAM_end + 1)]

                lPos_num           = lWithout_target_pos + lWith_target_pos + list(self.strPamSeq) + lAfter_PAM_pos
                lPos_annotated_ref = [str(i)+'.'+str(j) for i,j in zip(sRef_target, lPos_num)]
                ## ['A.-7', 'C.-6', 'A.-5', 'A.-4', 'G.-3', 'C.-2', 'A.-1', 'T.1', 'G.2', 'C.3', 'A.4', 'A.5', 'T.6', 'C.7', 'A.8', 'C.9', 'C.10', 'T.11', 'T.12', 'G.13', 'G.14',

                lMasked_pos_annotated_ref_target = []   ## '' '' '' A '' '' '' A A '' ''

                for sBase_pos in lPos_annotated_ref:
                    sBase_only = sBase_pos.split('.')[0]
                    if sBase_only != sTarget_ref:
                        lMasked_pos_annotated_ref_target.append(' ')
                    else:
                        lMasked_pos_annotated_ref_target.append(sBase_pos)

                #set_trace()

                strFormat = "{sample}\t{bar}\t{ref}\t{NumTot}\t{NumIns}\t{NumDel}\t{NumCom}\t{BaseEditCount}\n"
                ## Making a header
                Summary.write("Sample\tBarcode\tRef\t# of Total\t# of Insertion\t# of Deletion\t# of Combination\t{refseq}\n".format(refseq='\t'.join(lPos_annotated_ref)))
                Target_summary.write("Sample\tBarcode\tRef\t# of Total\t# of Insertion\t# of Deletion\t# of Combination\t{refseq}\n".format(refseq='\t'.join(lMasked_pos_annotated_ref_target)))

                for i, lBase_count in enumerate(zip(*lSum)):  ## lBase_count [(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0)]

                    if i == 0:
                        Summary.write(strFormat.format(sample=self.strFileName, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[self.intNumOfTotal], NumIns=lValue[self.intNumOfIns], NumDel=lValue[self.intNumOfDel], NumCom=lValue[self.intNumOfCom],
                            BaseEditCount='\t'.join(map(str, lBase_count))))
                    else:
                        Summary.write("\t\t\t\t\t\t\t{BaseEditCount}\n".format(BaseEditCount='\t'.join(map(str, lBase_count))))

                try:
                    lTarget_base_count = zip(*lSum)[iTarget_base]
                    lMasked_target_base_count = []  ## '' 20 '' 30 '' '' '' '' 20 ''

                    for sMasked_ref, fCount in zip(lMasked_pos_annotated_ref_target, lTarget_base_count):

                        if sMasked_ref == ' ':
                            lMasked_target_base_count.append(' ')
                        else:
                            lMasked_target_base_count.append(fCount)

                    Target_summary.write((strFormat.format(sample=self.strFileName, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[self.intNumOfTotal],
                                                           NumIns=lValue[self.intNumOfIns], NumDel=lValue[self.intNumOfDel], NumCom=lValue[self.intNumOfCom],
                        BaseEditCount='\t'.join(map(str, lMasked_target_base_count)))))

                except IndexError:
                    print('Null query: ', self.strForwPath)
                    ## Null query base count is all zero.
                    Target_summary.write(
                        (strFormat.format(sample=self.strFileName, bar=sBarcode, ref=sRef_seq_without_barcode, NumTot=lValue[self.intNumOfTotal],
                            NumIns=lValue[self.intNumOfIns], NumDel=lValue[self.intNumOfDel], NumCom=lValue[self.intNumOfCom],
                            BaseEditCount='\t'.join(['0'] * len(lPos_annotated_ref)))))


def Main():

    InstParameter = clsParameter()
    logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                        level=logging.DEBUG,
                        filename=InstParameter.strLogPath,
                        filemode='a')

    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    # Output: 1. Count information of matched barcode e.g. TACGATCTA\t# total\tins\t# del\t# com
    # Output: 2. classify FASTQ.    e.g. TAGAATATACACG.insertion.fastq

    logging.info('Program start : %s' % InstParameter.strFileName)

    InstParser = clsBaseEditParser(InstParameter)
    logging.info('File Open : %s' % InstParameter.strFileName)
    listSequenceForward = InstParser.OpenSequenceFiles()

    logging.info('Calculate base edit frequency : %s' % InstParameter.strFileName)
    dictResultForward  = InstParser.CalculateBaseEditFreq(listSequenceForward)

    logging.info('Make output forward : %s' % InstParameter.strFileName)
    InstOutput = clsOutputMaker(InstParameter)
    InstOutput.MakeOutput(dictResultForward)

    logging.info('Program end : %s' % InstParameter.strFileName)
# end: def Main


if __name__ == '__main__':
    Main()


