#!/home/hkim/anaconda2/bin/python2.7

from pdb import set_trace

## D0 Sub list
"""
Euchromatin_206_repeat5 TCTATCGTACATCGC Euchromatin_206_repeat5:39M2D_AC        1
ExtremeGC_811   CTACATCGTCATACA ExtremeGC_811:39M1D_G   1
"""
#strSubHiseq1 = './Sub_indel_result/Summation_Project_list_sub_indel.txt'  ## total indel cnt : 8929
#strSubHiseq2 = './Sub_indel_result/Summation_Project_list2_sub_indel.txt' ## total indel cnt : 8367
#strSubNeon1  = './Sub_indel_result/Summation_Project_list3_sub_indel.txt' ## total indel cnt : 9396
#3strSubNeon2  = './Sub_indel_result/Summation_Project_list4_sub_indel.txt' ## total indel cnt : 8321
strSubHiseq1 = './Output/Summation_Project_list_sub_indel.txt'  ## total indel cnt : 8929
strSubHiseq2 = './Output/Summation_Project_list2_sub_indel.txt' ## total indel cnt : 8367
strSubNeon1  = './Output/Summation_Project_list3_sub_indel.txt' ## total indel cnt : 9396
strSubNeon2  = './Output/Summation_Project_list4_sub_indel.txt' ## total indel cnt : 8321

## Total sum file
"""
       Sample  Barcode Ref     # of Total      # of Insertion  # of Deletion   # of Combination        A.-7  
0       Doench2014_1    CGCATATCATCATCA TAGATTGAAGAGAGACAGTACATGCCCTGGGAGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA       322     0       0       0       0      
"""
strTotalHiseq1 = './Output/Summation_Project_list.txt'
strTotalHiseq2 = './Output/Summation_Project_list2.txt'
strTotalNeon1  = './Output/Summation_Project_list3.txt'
strTotalNeon2  = './Output/Summation_Project_list4.txt'

## Freq result file
"""
Filename        Seq     Motif   Count   Total_cnt       Proportion      Substitution
Doench2014_1    CGCATATCATCATCATAGATTGAAGAGAGACAGTACATGCCCTGGGAGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA
        AGAGAGACA       246     257     0.9572  wt
"""
strFreqHiseq1 = './Output/Group_result/180903_split_hiseq_R1/Seq_freq.txt'
strFreqHiseq2 = './Output/Group_result/180903_split_hiseq_R2/Seq_freq.txt'
strFreqNeon1  = './Output/Group_result/190311_Neon_splitBE4_R1/Seq_freq.txt'
strFreqNeon2  = './Output/Group_result/190311_Neon_splitBE4_R2/Seq_freq.txt'

## Result
"""
Filename        Seq     Motif   Count   Total_cnt       Proportion      Substitution
Doench2014_1    CGCATATCATCATCATAGATTGAAGAGAGACAGTACATGCCCTGGGAGAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA
        AGAGAGACA       246     257-(D0 indel count)     0.9572  wt
-> next line
  + indelcompelex sum count        
        
"""
strResultHiseq1 = './Output/Seq_freq_add_info_Hiseq1.txt'
strResultHiseq2 = './Output/Seq_freq_add_info_Hiseq2.txt'
strResultNeon1  = './Output/Seq_freq_add_info_Neon1.txt'
strResultNeon2  = './Output/Seq_freq_add_info_Neon2.txt'


def Add_info_result(strSub, strTotal, strFreq, strResult):

    with open(strSub) as Sub,\
        open(strTotal) as Total,\
        open(strFreq) as Freq,\
        open(strResult, 'w') as Result:

        dictSubCnt = {}  ## Doench2016_1948:39M1D_G, Doench2016_1948:39M1I_G -> Neon: two file name is same but the pattern is different.
                         ## I should merge these pattern based on the file name.
        dictTotalAppliedSub = {}
        dictIndelSum        = {}

        for strRow in Sub:
            listCol  = strRow.replace('\n', '').split('\t')
            strFile  = listCol[0]
            intCount = int(listCol[3])

            try:
                dictSubCnt[strFile] += intCount
            except KeyError:
                dictSubCnt[strFile] = intCount

        # intSubIndelAllCnt = sum([v for k,v in dictSubCnt.items()])
        # print('%s sub indel all count : %s' % (strSub, intSubIndelAllCnt)) ## checked all count is correct.

        #"""
        for i, strRow in enumerate(Total):
            if i == 0: continue ## header skip.
            listCol  = strRow.replace('\n', '').split('\t')
            strFile  = listCol[1]
            intTotal = int(listCol[4])
            intIns   = int(listCol[5])
            intDel   = int(listCol[6])
            intCom   = int(listCol[7])

            try:
                intSub  = dictSubCnt[strFile]
            except KeyError:
                intSub = 0

            intTotalAppliedSub = intTotal - intSub            ## The total count is not subtracted by DOindel, so apply it.

            dictTotalAppliedSub[strFile] = intTotalAppliedSub
            dictIndelSum[strFile] = intIns + intDel + intCom  ## each file row indel complex count sum

        dictFreq = {}  ## {'GECKO_346': [[Filename	Seq	Motif	Count	Total_cnt	Proportion	Substitution],[],[],[]]}
        strHeader = ''

        for i, strRow in enumerate(Freq):                        ## Freq total was removed by crispr indel.
            if i == 0:
                strHeader = strRow
                continue             ## header skip.
            listCol = strRow.replace('\n', '').split('\t')
            strFile  = listCol[0]
            intCount = int(listCol[3])
            intTotal = int(listCol[4])
            floProp  = float(listCol[5])
            listCol[3] = intCount
            listCol[4] = intTotal
            listCol[5] = floProp

            try:
                dictFreq[strFile].append(listCol)
            except KeyError:
                dictFreq[strFile] = [listCol]

        Result.write(strHeader.replace('\n','')+'\tTotal(D0)\tD0_indel\n')

        for strFile, list2Col in dictFreq.items():

            list2Col[1:] = sorted(list2Col[1:], key=lambda x: x[5], reverse=True) ## sort by proportion

            intAltAllCnt = sum([listAlt[3] for listAlt in list2Col[1:]])

            ## for validation
            listCountCheck = []
            intTotalCheck  = 0

            for i, listCol in enumerate(list2Col):
                strSubstitution = listCol[6]
                intTotal        = listCol[4]

                intTotalAppliedSub = dictTotalAppliedSub[strFile]

                intIndelSum = dictIndelSum[strFile]   ## intIns + intDel + intCom
                intTotalD0  = intTotal + intIndelSum  ## freq total are substrated by indel sum, so add it again

                intD0IndelCount = intTotalD0 - intTotalAppliedSub

                if strSubstitution == 'wt':  ## modify WT count. Total - alt count = wt count
                    intModiCount = intTotalAppliedSub - intAltAllCnt - intIndelSum
                    if intModiCount < 0:
                        print('minus value error, this integer is positive.')
                        set_trace()

                    listCol[3]   = intModiCount
                    #if listCol[0] == 'GECKO_7232': Neon1, 2761
                    #    set_trace()
                listCountCheck.append(listCol[3])  ## for validation

                listCol[4] = intTotalAppliedSub
                try:
                    listCol[5] = round(float(listCol[3]) / listCol[4], 4)
                except Exception:
                    listCol[5] = 0

                Result.write('\t'.join(map(str, listCol + [intTotalD0, intD0IndelCount]))+'\n')

                if i == 0:
                    listCountCheck.append(intIndelSum)
                    intTotalCheck = listCol[4] ## for validation

                    listResultCol = len(listCol) * ['~'] + ['~', '~']
                    listResultCol[0] = strFile
                    listResultCol[6] = 'Indel'
                    listResultCol[3] = intIndelSum

                    Result.write('\t'.join(map(str, listResultCol))+'\n')

            #if strFile == 'GECKO_7232':
            #    set_trace()
            intCountCheckTotal = sum(listCountCheck)
            if intCountCheckTotal != intTotalCheck:
                print('Count total is diffrent. result:%s, file:%s, CountCheckTotal:%s, TotalCheck:%s' % (strResult, strFile, intCountCheckTotal, intTotalCheck))
        #"""


def Main():

    for strSub, strTotal, strFreq, strResult in [[strSubHiseq1, strTotalHiseq1, strFreqHiseq1, strResultHiseq1],
                                      [strSubHiseq2, strTotalHiseq2, strFreqHiseq2, strResultHiseq2],
                                      [strSubNeon1, strTotalNeon1, strFreqNeon1, strResultNeon1],
                                      [strSubNeon2, strTotalNeon2, strFreqNeon2, strResultNeon2]]:

        Add_info_result(strSub, strTotal, strFreq, strResult)

Main()

"""
## deprecated
def Merge_sub_indel_and_dict(strInput1, strInput2):
    dictSubIndel = {}

    with open(strInput1) as Input1, \
            open(strInput2) as Input2:

        for strRow in Input1:
            listCol = strRow.replace('\n', '').split('\t')
            strFile = listCol[0]
            strBarcode = listCol[1]
            strPattern = listCol[2]
            intCount = int(listCol[3])

            dictSubIndel[strFile] = [strBarcode, strPattern, intCount]

        for strRow in Input2:
            listCol = strRow.replace('\n', '').split('\t')
            strFile = listCol[0]
            intCount = int(listCol[3])

            dictSubIndel[strFile][2] += intCount
"""