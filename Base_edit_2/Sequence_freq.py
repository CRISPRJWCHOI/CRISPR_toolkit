#!/home/hkim/anaconda2/bin/python2.7

import os,sys
import numpy as np
from collections import Counter
from collections import OrderedDict
import multiprocessing as mp

import logging
from pdb import set_trace
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S', level=logging.DEBUG)

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import Helper


try:
    strUser    = sys.argv[1]
    strProject = sys.argv[2]
    lWindow    = sys.argv[3].split('-')
    iWinStart  = int(lWindow[0])
    iWinEnd    = int(lWindow[1])
    iCore      = int(sys.argv[4])


except IndexError:
    print('\nUsage: ./Sequence_freq.py SH 24K_screening 25-33 10\n'
          '         ./Sequence_freq.py user_name project_name window_range thread\n')
    sys.exit()


def Count_seq_freq(lPara):

    """ aligned_BaseEdit.txt
    ACTAGCTATCGCTCACTCTGGGGTCAGGGACAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA        CGCTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGACTAGCTATCGCTCACTCTGGGGTCAGGGGCAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAATA  []
    [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
    ACTAGCTATCGCTCACTCTGGGGTCAGGGACAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTA--A      ACTAGCTATCGCTCACTCTGGGGTCAGGGGCAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAATA
    """

    strSample           = lPara[0]
    sFile_path          = lPara[1]
    sTotal_readcnt_path = lPara[2]
    dInput_fa           = lPara[3]
    print ("Count_seq_freq: ", strSample, sFile_path, sTotal_readcnt_path)

    try:
        with open('./Output/{user}/{project}/{sample}/Result/Seq_freq.txt'.format(user=strUser,
                                                                                   project=strProject,
                                                                                   sample=strSample), 'w') as Output:

            Output.write('Filename\tSeq\tMotif\tCount\tTotal_cnt\tProportion\tSubstitution\n')

            ## A project has many file. The total read count is summation in the each file.
            for iFile_num, sFile in enumerate(os.listdir(sFile_path)):
                #set_trace()
                if 'aligned' in sFile:
                    # print(iFile_num)
                    sFilename = sFile.replace('_aligned_BaseEdit', '').split('.')[:-1][0] ## To search filename in other folder.
                    sTotal_readcnt_file = sFilename + '_Summary.txt'                      ## extract totral read count for the sequence frequency.

                    with open(sFile_path + '/' + sFile) as aligned_BaseEdit,\
                        open(sTotal_readcnt_path + '/' + sTotal_readcnt_file) as Total_readcnt:
                        #print(sFile_path + '/' + sFile)

                        iTotal_wo_indel = 0

                        for i, sRow in enumerate(Total_readcnt):
                            if i == 0: continue
                            lCol = sRow.replace('\n', '').split('\t')
                            iTotal_read = int(lCol[3])                   ## This is read counts of a matched barcode.
                            iIndel_read = int(lCol[4]) + int(lCol[5]) + int(lCol[6])
                            iTotal_wo_indel = iTotal_read - iIndel_read  ## Total read is without indel reads
                            break                                        ## 2 row is target, over 3 is none

                        lTarget_seq        = []
                        sRef_seq           = ''
                        dSeq_wt_extend     = {}  ## WT + motif(target sequence) + WT

                        for i, sRow in enumerate(aligned_BaseEdit):

                            lCol       = sRow.replace('\n', '').split('\t')
                            sQuery_seq = lCol[5]

                            if sRef_seq == '':   ## Reference is same in the file, so store once.
                                sRef_seq = lCol[0]
                                dSeq_wt_extend[sRef_seq[iWinStart - 1: iWinEnd]] = sRef_seq

                            lRef_seq_with_motif = list(sRef_seq)
                            lRef_seq_with_motif[ iWinStart-1 : iWinEnd ] = list(sQuery_seq[ iWinStart-1 : iWinEnd ])
                            sRef_seq_with_motif = ''.join(lRef_seq_with_motif)

                            dSeq_wt_extend[sQuery_seq[ iWinStart-1 : iWinEnd ]] = sRef_seq_with_motif
                            lTarget_seq.append(sQuery_seq[ iWinStart-1 : iWinEnd ])

                        iNormal  = iTotal_wo_indel - len(lTarget_seq)
                        sRef_seq = sRef_seq[ iWinStart-1 : iWinEnd ]
                        dSeq_cnt = Counter(lTarget_seq)

                        try:
                            iRef_cnt_in_aligned = dSeq_cnt[sRef_seq]  ## check normal sequence because substitution exists outside of window size.
                            iNormal = iNormal + iRef_cnt_in_aligned
                            del dSeq_cnt[sRef_seq]
                        except KeyError:
                            pass

                        if iNormal > 0:
                            if sRef_seq == '':                     ## aligned result file can be none result file. So extract from input file.
                                sRef_seq = dInput_fa[sFilename][1] ## dInput_fa[0] : full ref, dInput_fa[1] : target ref
                                dSeq_wt_extend[sRef_seq] = dInput_fa[sFilename][0]
                            try:
                                Output.write('\t'.join(map(str, [sFilename, dSeq_wt_extend[sRef_seq], sRef_seq, iNormal, iTotal_wo_indel, round(iNormal/float(iTotal_wo_indel),4), 'ref_from_result']))+'\n')
                            except Exception as e:
                                print(e, 'line150')
                                set_trace()

                        elif iNormal == 0:  ## if iNormal = 0, that means no result generation. because aligned_BaseEdit file is not contained non-read file.
                            sRef_seq = dInput_fa[sFilename][1]
                            dSeq_wt_extend[sRef_seq] = dInput_fa[sFilename][0]
                            try:
                                Output.write('\t'.join(map(str, [sFilename, dSeq_wt_extend[sRef_seq], sRef_seq, iNormal, iTotal_wo_indel, iNormal, 'ref_from_input'])) + '\n')
                            except Exception as e:
                                print(e, 'line158')
                                set_trace()

                        for sSeq, iCnt in dSeq_cnt.most_common():
                            try:
                                Output.write('\t'.join(map(str, [sFilename, dSeq_wt_extend[sSeq], sSeq, iCnt, iTotal_wo_indel, round(iCnt/float(iTotal_wo_indel),4), 'alt']))+'\n')
                            except Exception as e:
                                print(lPara[0], sFilename)
                                print(iCnt, iTotal_wo_indel)
                                print(e, 'line175')
                                #pass
                                set_trace()
                    #END: for
                #END: with
            #END: for
        #END: with
    except Exception as e:
        print(e)
        print("Error in the input: ", strSample, sFilename, sTotal_readcnt_file)
        pass
#END: def


def Make_ref_dict(strRef):

    dInput_fa = {}

    with open('./Input/{user}/Reference/{project}/{ref}/Reference.fa'.format(user=strUser,
                                                                             project=strProject,
                                                                             ref=strRef)) as Input_ref:
        """
        YSKim_0525+01614_98_repeat1     TATACACGCATGTAT TTTGTATACACGCATGTATGCATCCTGCAGGTCTCGCTCTGACATGTGGGAAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA

        1 file has 1 barcode.
        This should be done.
        """
        for sRow in Input_ref:
            lCol = sRow.replace('\n', '').split('\t')

            sInputFile = lCol[0]
            sBarcode   = lCol[1]
            sInputRef  = lCol[2]

            iBarcode_start        = sInputRef.index(sBarcode)
            sBarcode_start_ref    = sInputRef[iBarcode_start:]
            dInput_fa[sInputFile] = [sBarcode_start_ref, sBarcode_start_ref[iWinStart - 1: iWinEnd]]

    return dInput_fa


def Count_group():

    """
    Filename        Seq     Count   Total_cnt       Proportion      Substitution
    Doench2014_1000 AGGGACA 13      14      0.9286  ref_from_result
    Doench2014_1000 AG----- 1       14      0.0714  alt
    Doench2014_1001 GGCGCCA 17      26      0.6538  ref_from_result
    Doench2014_1001 GGTGCCA 5       26      0.1923  alt
    Doench2014_1001 GGAGCCA 2       26      0.0769  alt
    Doench2014_1001 GGCGCTA 1       26      0.0385  alt
    """

    sHeader    = ''
    dTotal_cnt = {}

    ## Make dictionary to sum the total reads count of the group. The total reads count is always same in their group.
    with open('Group_list.txt') as Group_list:
        for sGroupname in Group_list:
            if sGroupname[0] == "#": continue
            sGroupname = sGroupname.replace('\n', '').strip()
            if not os.path.isdir('./Output/Group_result'): os.mkdir('./Output/Group_result')

            for sFile in os.listdir('./Output'):
                if sGroupname in sFile:  ## matched group names -> Sum the counts
                    with open('./Output/%s/Summary/Seq_freq.txt' % sFile) as SeqFreq:

                        sHeader = SeqFreq.readline()
                        dSelect_one_total_cnt = {}
                        for sRow in SeqFreq:
                            lCol            = sRow.replace('\n', '').split('\t')
                            sFilename       = lCol[0]
                            try:
                                iTotal_read_cnt = int(lCol[4])
                            except IndexError:
                                set_trace()
                            dSelect_one_total_cnt[sFilename] = iTotal_read_cnt

                        for sFilename, iTotal_read_cnt in dSelect_one_total_cnt.items():
                            try:
                                dTotal_cnt[sGroupname + '_' + sFilename] += iTotal_read_cnt
                            except KeyError:
                                dTotal_cnt[sGroupname + '_' + sFilename] = iTotal_read_cnt

    with open('Group_list.txt') as Group_list:
        for sGroupname in Group_list:
            if sGroupname[0] == "#": continue
            sGroupname = sGroupname.replace('\n', '').strip()
            dSeq_freq  = OrderedDict()   ## ('GECKO_6367_GATCTGCTC', ['GECKO_6367', 'GATCTGCTC', 2, 156, '0.0128']),
                                         ## Unique key, only one list.
            if not os.path.isdir('./Output/Group_result'): os.mkdir('./Output/Group_result')

            for sFile in os.listdir('./Output'):
                if sGroupname in sFile:  ## matched group names -> Sum the counts
                    with open('./Output/%s/Summary/Seq_freq.txt' % sFile) as SeqFreq:

                        sHeader = SeqFreq.readline()

                        for sRow in SeqFreq:
                            lCol            = sRow.replace('\n', '').split('\t')
                            sFilename       = lCol[0]
                            sSeq_wt_extend  = lCol[1]
                            sFile_seq       = lCol[0] + '_' + lCol[2]  ## Unique name : Doench2014_1000_CTCTGGGGT
                            iCount          = int(lCol[3])
                            iTotal_read_cnt = dTotal_cnt[sGroupname + '_' + sFilename]

                            lCol[3] = iCount
                            lCol[4] = iTotal_read_cnt

                            try:
                                _ = dSeq_freq[sFile_seq]

                                dSeq_freq[sFile_seq][3] += iCount
                                #dSeq_freq[sFile_seq][4] = iTotal_read_cnt

                            except KeyError:
                                dSeq_freq[sFile_seq] = lCol  ## initial assignment

            ## x[0] : key, x[1] : value, int(x[1][5]) : proportion, x[1][6]: alt, wt category, x[1][0]: filename,
            llSeq_freq = sorted(sorted(dSeq_freq.items(), key=lambda x:x[1][6], reverse=True), key=lambda x:x[1][0])
            if not os.path.isdir('./Output/Group_result/%s' % sGroupname): os.mkdir('./Output/Group_result/%s' % sGroupname)
            with open('./Output/Group_result/%s/Seq_freq.txt' % sGroupname, 'w') as Output:

                Output.write(sHeader)

                for sFile_seq, lCol in llSeq_freq:
                    try:
                        try:
                            lCol[5] = round(float(lCol[3])/lCol[4], 4)  ## proportion calculation, previous proportion is not correct.
                        except ZeroDivisionError:
                            lCol[5] = 0
                    except Exception:
                        set_trace()
                    Output.write('\t'.join(map(str, lCol)).replace('ref_from_result', 'wt').replace('ref_from_input', 'wt')+'\n')
        #END: for
    #END: with


def Trim_data():

    """
    Remove gap seqs (e.g. AC---)
    """
    with open('Group_list.txt') as Group_list:
        for sGroupname in Group_list:
            if sGroupname[0] == "#": continue
            sGroupname = sGroupname.replace('\n', '').strip()
            dSeq_freq  = OrderedDict()

            with open('./Output/Group_result/%s/Seq_freq.txt' % sGroupname) as Group_result,\
                open('./Output/Group_result/%s/Trimmed_seq_freq.txt' % sGroupname, 'w') as Trimmed_result:

                sHeader = ''

                for i, sRow in enumerate(Group_result):

                    if i == 0:
                        sHeader = sRow
                        continue

                    lCol            = sRow.replace('\n', '').split('\t')
                    sFilename       = lCol[0]  ## Doench2014_1000

                    try:
                        dSeq_freq[sFilename].append(lCol)
                    except KeyError:
                        dSeq_freq[sFilename] = [lCol]

                for sFilename in dSeq_freq:
                    llFilename = dSeq_freq[sFilename]  ## [[Doench2014_1000,ACAGCAGCGAAC...,ACGCATC, 12,30,0.4][],[]...
                                                           ## A Same file name chunk in the group file.
                    iRecal_total      = 0  ## sub the gap seq cnt
                    #lDele_key     = []
                    llPre_recal_total = []
                    llRecal_total     = []

                    for i, lFilename in enumerate(llFilename):
                        sMotif          = lFilename[2]
                        iMotif_cnt      = int(lFilename[3])
                        iTotal_read_cnt = int(lFilename[4])

                        if lFilename[6] == 'wt':
                            iRecal_total = iTotal_read_cnt
                            llPre_recal_total.append(lFilename)

                        elif '-' in sMotif:
                            iRecal_total -= iMotif_cnt
                            continue
                        else:
                            llPre_recal_total.append(lFilename) ## store AC----- row key

                    for lPre_recal_total in llPre_recal_total:
                        lPre_recal_total[4] = iRecal_total
                        try:
                            lPre_recal_total[5] = round(float(lPre_recal_total[3])/iRecal_total,4)  ## recal proportion because of sub.
                        except ZeroDivisionError:
                            pass
                        llRecal_total.append(lPre_recal_total)

                    #llRecal_total[1:]    = sorted(llRecal_total[1:], key=lambda x: float(x[5]), reverse=True)
                    dSeq_freq[sFilename] = llRecal_total ## reassign the total cnt
                #END for

                llFilename_chunk = sorted(dSeq_freq) ## key is a filename
                for sKey in llFilename_chunk:
                    llCol = dSeq_freq[sKey]
                    llCol = sorted(llCol, key=lambda x: x[6], reverse=True) ## wild type category first

                    if llCol[0][6] != 'wt':
                        logging.critical('error, wildtype must be fisrt row. If you see this error message, please contact the developer.')
                        logging.critical('This program will be terminated.')
                        sys.exit()

                    if len(llCol) > 1:  ## It has alt. only a wt file does not necessary.
                        llCol[1:] = sorted(llCol[1:], key=lambda x: float(x[5]), reverse=True)
                        dSeq_freq[sKey] = llCol

                Trimmed_result.write(sHeader)
                for llRecal_total_final in dSeq_freq.values():
                    for lRecal_total_final in llRecal_total_final:
                        Trimmed_result.write('\t'.join(map(str,lRecal_total_final))+'\n')
            #END with
        #END for
    #END with


def Main():

    logging.info('Program Start')

    logging.info('Make commands for a multiple processing')
    lPara = []
    with open('./User/{user}/{project}.txt'.format(user=strUser, project=strProject)) as Project_list:

        for strSample in Project_list:
            if strSample[0] == '#': continue

            tupSampleInfo = Helper.SplitSampleInfo(strSample)
            if not tupSampleInfo: continue
            strSample, strRef, strExpCtrl = tupSampleInfo

            strSample           = strSample.replace('\n', '').replace('\r', '')
            sFile_path          = './Output/{user}/{project}/{sample}/Tmp/Alignment'.format(user=strUser, project=strProject, sample=strSample)
            sTotal_readcnt_path = './Output/{user}/{project}/{sample}/Tmp/All'.format(user=strUser, project=strProject, sample=strSample)
            dInput_fa           = Make_ref_dict(strRef)

            lPara.append([strSample, sFile_path, sTotal_readcnt_path, dInput_fa])

    logging.info('Multiple processing Start')
    p = mp.Pool(iCore)
    p.map_async(Count_seq_freq, lPara).get()
    logging.info('Multiple processing End')

    #logging.info('Count group Start')
    #Count_group()
    #logging.info('Count group End')

    #logging.info('Trim data Start')
    #Trim_data()
    #logging.info('Trim data End')

    logging.info('Program End')


Main()
