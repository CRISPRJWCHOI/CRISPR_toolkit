#!/home/hkimlab/anaconda2/bin/python2.7

import os, sys
from pdb import set_trace

try:
    strUser        = sys.argv[1]
    strProject     = sys.argv[2]

except IndexError:
    print('\n')
    print('usage   : ./Each_base_summary.py user_name project_name\n')
    print('example : ./Each_base_summary.py SH p53_screening\n')
    sys.exit()


def Make_target_ref_alt_summary(strSample='', strRef='', strAlt='', strFirstOutput=''):

    """ row 0: header, 1: A and info, 2: C, 3: G, 4: T
    Sample          Barcode         Ref                                                                                          # of Total # of Insertion # of Deletion  # of Combination  C.-7    T.-6    C.-5    T.-4    G.-3    G.-2    G.-1    G.1     T.2     C.3     A.4     G.5     G.6     G.7      A.8     C.9     A.10    G.11    T.12    G.13    G.14    A.15    C.16    T.17    C.18    G.19    A.20    A.N     G.G     G.G     A.1     G.2     A.3
    Doench2014_1000 ACTAGCTATCGCTCA CTCTGGGGTCAGGGACAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA       5       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       00       0
                                                                                                                                            0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       00       0       0       0       0       0       0       0       0       0       0       0       0       0
                                                                                                                                            0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       00       0       0       0       0       0       0       0       0       0       0       0       0       0
                                                                                                                                            0       0       0       0       0       0       0       0       0       0       0       1       0       0       0       0       0       0       00       0       0       0       0       0       0       0       0       0       0       0       0       0
    """
    dAlt = {'A' : 1, 'C' : 2, 'G' : 3, 'T' : 4}
    lHeader = []
    llResult = []

    strSampleDir  = './Output/{user}/{project}/{sample}'.format(user=strUser,
                                                                project=strProject,
                                                                sample=strSample)
    strSummaryDir = os.path.join(strSampleDir, 'Result')
    strMergeTargetDir = os.path.join(strSummaryDir, 'Merge_target_result')

    with open(os.path.join(strMergeTargetDir,strFirstOutput)) as Fisrt_output,\
        open(os.path.join(strMergeTargetDir, '{sample}_{ref}to{alt}_Summary_addition.txt'.format(sample=strSample,
                                                                                                 ref=strRef,
                                                                                                 alt=strAlt)), 'w') as Output:

        strSummaryAllDir = os.path.join(strSampleDir,'Tmp/All')
        for iFile_cnt, sFile in enumerate(os.listdir(strSummaryAllDir)):

            with open(os.path.join(strSummaryAllDir, sFile)) as Input:
                lNone_alt_col  = []
                lBaseEdit_Info = []

                for i, sRow in enumerate(Input):
                    lCol = sRow.replace('\n', '').split('\t')

                    if i == 0:
                        for j, sCol_name in enumerate(lCol[7:]):
                            if strRef not in sCol_name:
                                lNone_alt_col.append(7+j)
                                lCol[7+j] = ' '

                        if lHeader == []:
                            lHeader = lCol
                        elif lHeader:
                            for iHeader_col, tHeader in enumerate(zip(lHeader[7:], lCol[7:])):
                                sHeader_current, sHeader_update = tHeader

                                if sHeader_update == ' ': continue

                                if sHeader_current == ' ':
                                    lHeader[iHeader_col+7] = sHeader_update

                                else:
                                    assert  sHeader_current == sHeader_update, 'Check header %s %s' % (repr(sHeader_current), repr(sHeader_update))

                    elif i == 1:
                        lBaseEdit_Info = lCol[:7]

                    elif i == dAlt[strAlt]:
                        for iNon_col in lNone_alt_col:
                            lCol[iNon_col] = ' '
                        lCol[:7] = lBaseEdit_Info
                        #print(i, lCol)
                        #(3, ['Doench2014_1000', 'ACTAGCTATCGCTCA', 'CTCTGGGGTCAGGGACAGTGGACTCGAAGGAGAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA', '5', '0', '0', '0', '', '', '', '', '', '', '', '', '', '', '0', '', '', '', '0', '', '0', '', '', '', '', '0', '', '', '', '', '0', '0', '', '', '0', '', '0'])
                        llResult.append(lCol)


        print('Total_files: ', iFile_cnt + 1)
        Output.write('\t'.join(lHeader) + '\n')

        """
        All folder doesn't able to have any indel information if it hasn't any counts of alterantive alleles.
        That file has only a header.
        Hence, I check the first merged summary output data, then extract it doesn't have current additional output.   
        """

        dAdditional_output = {} ## dictionary to check for only header files in the 'all' folder.

        for lResult in llResult:
            sSample = lResult[0]
            dAdditional_output[sSample] = '\t'.join(lResult) + '\n'

        for i, sRow in enumerate(Fisrt_output):
            if i == 0: continue ## header skip
            lCol    = sRow.replace('\n', '').split('\t')
            sSample = lCol[0]

            try:
                Output.write(dAdditional_output[sSample])
            except KeyError:     ## Exclusive possession
                Output.write(sRow)


def Main():

    with open('./User/{user}/Additional_BaseEdit_process_list.tsv'.format(user=strUser)) as Input:
        for sRow in Input:
            if sRow[0] == '#': continue
            lCol = sRow.replace('\n', '').replace('\r', '').split('\t')
            if len(lCol) == 1:
                lCol = lCol[0].split()
            print(lCol)

            strSample      = lCol[0]
            listRefAlt     = lCol[1].split(',')
            strRef         = listRefAlt[0]
            strAlt         = listRefAlt[1]
            strFirstOutput = lCol[2]

            Make_target_ref_alt_summary(strSample=strSample,
                                        strRef=strRef,
                                        strAlt=strAlt,
                                        strFirstOutput=strFirstOutput)

Main()
