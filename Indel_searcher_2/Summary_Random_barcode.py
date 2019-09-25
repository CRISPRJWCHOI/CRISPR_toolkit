import os, sys
import logging
import multiprocessing as mp

from argparse import ArgumentParser
from collections import OrderedDict

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import  Helper


class clsParameters():

    def __init__(self, options):
        self.strUser    = options.user_name
        self.strProject = options.project_name.replace('.txt', '') ## A user can be confused the input. So I prevented from it using 'replace'.
        self.strGroup   = options.group
        self.intCore    = options.thread

        self.strSampleList = 'User/{user}/{project}.txt'.format(user=options.user_name, project=options.project_name)


def SummaryRandomBarcode(sFile_path):

    """
    /Tmp
    190819_Nahye_24k_2_D0_2-24kLib_Classified_Indel_barcode.fastq* -> process target
    190819_Nahye_24k_2_D0_2-24kLib_Indel_freq.txt*
    190819_Nahye_24k_2_D0_2-24kLib_Indel_summary.txt*
    190819_Nahye_24k_2_D0_2-24kLib_Summary.txt*
    Pickle

    dBarcode_cnt = {'ACGTACTC_sorting_barcode': {'ACATACAC_random': 5, 'CGTGTTGA_random': 3, ...}
    """
    dictBarcodeCnt = {}
    strClassCheck  = ''

    strSample = sFile_path.split('/')[-1]
    logging.info('Summary_random_barcode start : %s, %s' % (sFile_path, strSample))

    for sFile in os.listdir(sFile_path+'/Tmp/'):
        if '.fastq' in sFile:
            with open(sFile_path+'/Tmp/'+sFile) as Input:
                for i, strRow in enumerate(Input):

                    # @D00235:683:CE1P6ANXX:6:1114:2135:5231 1:N:0:CTGAAGCT+CCTATCCT:Barcode_TTTGCTATCTCGACGTATGGACAGTG:total
                    if i % 4 == 0:
                        listBarClass = strRow.replace('\n','').split('Barcode_')[1].split(':')
                        strBarcode   = listBarClass[0]
                        strClass     = listBarClass[1]

                        if strClass == 'total':
                            strClassCheck = 'total'

                    if i % 4 == 1 and strClassCheck == 'total':
                        strRow = strRow.replace('\n','').upper()
                        intBarcodeStart   = strRow.find(strBarcode)
                        strRandom_barcode = strRow[intBarcodeStart-8:intBarcodeStart]

                        try:
                            _ = dictBarcodeCnt[strBarcode]
                        except KeyError:
                            dictBarcodeCnt[strBarcode] = {}
                        try:
                            dictBarcodeCnt[strBarcode][strRandom_barcode] += 1
                        except KeyError:
                            dictBarcodeCnt[strBarcode][strRandom_barcode] = 1
                        #print(sBarcode, sRandom_barcode, iBarcode_start, sRow)

                        strClassCheck = ''

    if not os.path.isdir(sFile_path + '/Summary_Random_barcode'): os.mkdir(sFile_path + '/Summary_Random_barcode')
    with open(sFile_path + '/Summary_Random_barcode/%s_all_random_barcode.txt' % strSample, 'w') as All_random,\
        open(sFile_path + '/Summary_Random_barcode/%s_Unique_RandomBarcodeNumber_In_SortingBarcode.txt' % strSample, 'w') as Random_sorting:

        All_random.write('Sorting_barcode\tUnique_RandomBarcodeNumber_In_SortingBarcode\tRandomBarcode\tEach_RandomBarcode_read_count\n')
        Random_sorting.write('Sorting_barcode\tUnique_RandomBarcodeNumber_In_SortingBarcode\n')

        for sBarcode, dRandom_barcode_cnt in dictBarcodeCnt.items():
            iRandom_barcode_num = len(dRandom_barcode_cnt.keys())
            Random_sorting.write('\t'.join(map(str, [sBarcode, iRandom_barcode_num]))+'\n')

            for sRandom_barcode, iCnt in dRandom_barcode_cnt.items():
                All_random.write('\t'.join(map(str, [sBarcode, iRandom_barcode_num, sRandom_barcode, iCnt]))+'\n')

    logging.info('Summary_random_barcode end: %s' % sFile_path)

## on going
def CountGroup(InstParameters):
    """
    Sorting_barcode Unique_RandomBarcodeNumber_In_SortingBarcode    RandomBarcode   Each_RandomBarcode_read_count
    TATATCATAGCGTACTCATC    8       TGCGTTTG        3
    TATATCATAGCGTACTCATC    8       CGCGTTTG        3
    TATATCATAGCGTACTCATC    8       TAGTTTTG        1
    TATATCATAGCGTACTCATC    8       ATAGTTTG        1
    """

    sHeader = ''

    with open(InstParameters.strSampleList) as Sample: ## tmp input

        listSample = Sample.readlines()

        setGroup = set([strRow.replace('\n', '').split('\t')[2].upper() for strRow in listSample])

        for strGroup in setGroup:
            if strGroup == 'CTRL': continue

            for strRow in listSample:
                if strGroup == strGroupOfSample:  ## matched group names -> Sum the counts
                    listCol          = strRow.replace('\n', '').split('\t')
                    strSample        = listCol[0]
                    strRef           = listCol[1]
                    strGroupOfSample = listCol[2]

                    strProjectDir = './Output/{user}/{project}'.format(user=InstParameters.strUser,
                                                                       project=InstParameters.strProject)
                    strGroupDir = os.path.join(strProjectDir, 'Group_result')
                    Helper.MakeFolderIfNot(strGroupDir)

                    dTotal_RandomBarcode_cnt_in_SortingBarcode = OrderedDict() ## ('GECKO_6367_GATCTGCTC', ['GECKO_6367', 'GATCTGCTC', 2, 156, '0.0128']),
                                                                               ## Unique key, only one list.

                    with open('{project_dir}/{sample}_all_random_barcode.txt'.format(project_dir=strProjectDir,
                                                                                     sample=strSample)) as RandomBarcode_SeqFreq:
                        sHeader = RandomBarcode_SeqFreq.readline()

                        for sRow in RandomBarcode_SeqFreq:
                            lCol = sRow.replace('\n', '').split('\t')

                            sSortingBarcode                             = lCol[0]
                            #iTotal_RandomBarcode_cnt_in_SortingBarcode  = int(lCol[1])
                            sSorting_and_Random_barcode_seq             = lCol[0] + '_' + lCol[2]  ## Unique name : Doench2014_1000_CTCTGGGGT
                            iRandomBarcode_count                        = int(lCol[3])

                            lCol[3] = iRandomBarcode_count

                            try:
                                _ = dTotal_RandomBarcode_cnt_in_SortingBarcode[sSorting_and_Random_barcode_seq]

                                dTotal_RandomBarcode_cnt_in_SortingBarcode[sSorting_and_Random_barcode_seq][3] += iRandomBarcode_count

                            except KeyError:
                                dTotal_RandomBarcode_cnt_in_SortingBarcode[sSorting_and_Random_barcode_seq] = lCol  ## initial assignment
                    #END for
                    dRecal_total_kind_of_RandomBarcode = OrderedDict()
                    for sSort_Rand_seq in dTotal_RandomBarcode_cnt_in_SortingBarcode:  ## sSorting_and_Random_barcode_seq
                        sSortBarcode = sSort_Rand_seq.split('_')[0]
                        try:
                            dRecal_total_kind_of_RandomBarcode[sSortBarcode].append(dTotal_RandomBarcode_cnt_in_SortingBarcode[sSort_Rand_seq])
                        except KeyError:
                            dRecal_total_kind_of_RandomBarcode[sSortBarcode] = [dTotal_RandomBarcode_cnt_in_SortingBarcode[sSort_Rand_seq]]

                    for sKey, llValue in dRecal_total_kind_of_RandomBarcode.items():
                        ## sKey: TATATCATAGCGTACTCATC, llValue : [[TATATCATAGCGTACTCATC, 8, TGCGTTTG, 3],[],[] ...
                        iKind_of_RandomBarcode = len(llValue)  ################## why do I make like this ?????
                        for lValue in llValue:
                            lValue[1] = iKind_of_RandomBarcode ## Recal using group total cnt.

                        llValue = sorted(llValue, key=lambda x:x[3], reverse=True)
                        dRecal_total_kind_of_RandomBarcode[sKey] = llValue

                    strEachGroup = './Output/Group_result/%s' % strGroup
                    Helper.MakeFolderIfNot(strEachGroup)

                    with open(os.path.join(strEachGroup, 'Summary_all_random_barcode_in_group.txt'), 'w') as Sort_Random_cnt,\
                        open(os.path.join(strEachGroup, 'Summary_Unique_RandomBarcodeNumber_in_group.txt'), 'w') as Uniq_random_cnt:

                        Sort_Random_cnt.write(sHeader)
                        Uniq_random_cnt.write('Sorting_barcode\tUnique_RandomBarcodeNumber_In_SortingBarcode\n')

                        for sSortBarcode, llCol in dRecal_total_kind_of_RandomBarcode.items():
                            Uniq_random_cnt.write('\t'.join(map(str, [sSortBarcode, len(llCol)]))+'\n')
                            for lCol in llCol:
                                Sort_Random_cnt.write('\t'.join(map(str, lCol))+'\n')
        #END: for
    #END: with


def Main():

    logging.info('Program Start')
    logging.info('Make commands for a multiple processing')

    parser = ArgumentParser(description='Script for counting the random barcodes')

    parser.add_argument('-u', '--user_name', type=str, dest='user_name', help='The user name in the /user subdir')
    parser.add_argument('-p', '--project_name', type=str, dest='project_name', help='The project name in the /user/user_name/ subdir')
    parser.add_argument('-g', '--group',  type=str, dest='group', default='false', help='The group sum run of the barcodes, default: false')
    parser.add_argument('-t', '--thread', type=int, dest='thread', default='15', help='The multicore number 1~15')
    options = parser.parse_args()

    InstParameters = clsParameters(options)

    lPara = []

    with open(InstParameters.strSampleList) as SampleList:

        for strSample in SampleList:
            if strSample[0] == '#' or strSample[0] in ['', ' ', '\r', '\n', '\r\n']: continue
            strSample  = strSample.replace('\n', '').replace('\r', '').split('\t')[0]
            sFile_path = './Output/{user}/{project}/{sample}'.format(user=options.user_name,
                                                                     project=options.project_name,
                                                                     sample=strSample)
            #print('sFile_path', sFile_path)
            lPara.append(sFile_path)

    ## single_test
    #Summary_random_barcode(lPara[0])

    logging.info('Multiple processing Start')
    p = mp.Pool(options.thread)
    p.map_async(SummaryRandomBarcode, lPara).get()
    logging.info('Multiple processing End')

    #logging.info('Count group Start')
    #CountGroup(InstParameters)
    #logging.info('Count group End')

    #logging.info('Program End')

Main()
