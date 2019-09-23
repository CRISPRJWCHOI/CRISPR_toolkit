#!/usr/bin/env python

import os, sys
import pandas as pd

from pdb import set_trace

strProjectList = sys.argv[1]
#strProjectList = 'Project_list2.txt'


def Summation_all_final_result():

    with open(strProjectList) as Input:

        listdfResult = []
        for i, strSample in enumerate(Input):
            #print(strSample)
            #if i == 2: break
            strSample = strSample.replace('\n','').replace('\r','').strip()
            strFinalResultDir = './Output/%s/Summary/Merge_target_result/' % strSample

            for j, strFinalResultFile in enumerate(os.listdir(strFinalResultDir)):
                if j > 0:
                    print('I expected one file, but there are more. check the target base change file')
                    sys.exit(1)
                    
                print(strFinalResultFile)
                strFinalResultPath = './Output/%s/Summary/Merge_target_result/%s' % (strSample, strFinalResultFile)

                listdfResult.append(pd.read_table(strFinalResultPath, low_memory=False))

        dfAll        = pd.concat(listdfResult)
        dfForw       = dfAll.iloc[:,0:3]
        dfReve       = dfAll.iloc[:,3:].replace(' ', '0').astype('int64')
        dfAllResult  = pd.concat([dfForw, dfReve], axis=1).groupby(['Sample','Barcode','Ref']).sum()
        dfAllResult.reset_index(inplace=True)

        dfAllResult.to_csv('./Output/Summation_'+strProjectList, sep='\t')

        #with open('./Output/%s/Summary/Merge_target_result/%s' % (strSample, strFinalResultFile)) as FinalResult:
        """
            for strRow in FinalResult:
                listCol       = strRow.replace('\n','').split('\t')
                listSamBarRef = listCol[:3]
                = listCol[3:]
        """


def SummationSubIndel():

    with open(strProjectList) as Input,\
        open('./Output/Summation_' + strProjectList.replace('.txt','') + '_sub_indel.txt', 'w') as Output:

        dictResult = {}

        for i, strSample in enumerate(Input):
            print(strSample)
            #if i == 2: break
            strSample = strSample.replace('\n','').replace('\r','').strip()
            strSubIndelDir = './Output/%s/result' % strSample

            for strSubIndelFile in os.listdir(strSubIndelDir):
                if 'sub' in strSubIndelFile:
                    with open(strSubIndelDir + '/' + strSubIndelFile) as SubIndel:
                        for strRow in SubIndel:
                            listCol         = strRow.replace('\n','').split('\t')
                            setIndelPattern = set(listCol[3].split(','))
                            intCount        = int(listCol[2])
                            strNameBarcodePattern  = '-'.join(listCol[0:2])+'-'+''.join(setIndelPattern)

                            try:
                                dictResult[strNameBarcodePattern] += intCount
                            except KeyError:
                                dictResult[strNameBarcodePattern] = intCount

        for strNameBarcodePattern, intCount in dictResult.items():
            Output.write('\t'.join(strNameBarcodePattern.split('-')) + '\t' + str(intCount) + '\n')


def ConfirmValidation():

    with open(strProjectList) as Input:

        listdfResult = []
        for i, strSample in enumerate(Input):
            if i == 2: break
            print(strSample)
            strSample = strSample.replace('\n','').replace('\r','').strip()
            strFinalResultDir = './Output/%s/Summary/Merge_target_result/' % strSample

            for strFinalResultFile in os.listdir(strFinalResultDir):
                print(strFinalResultFile)
                strFinalResultPath = './Output/%s/Summary/Merge_target_result/%s' % (strSample, strFinalResultFile)

                listdfResult.append(pd.read_table(strFinalResultPath, low_memory=False))

        dfAll        = pd.concat(listdfResult)
        dfForw       = dfAll.iloc[:,0:3]
        dfReve       = dfAll.iloc[:,3:].replace(' ', '0').astype('int64')
        dfAllResult  = pd.concat([dfForw, dfReve], axis=1).groupby(['Sample','Barcode','Ref']).sum()
        dfAllResult.reset_index(inplace=True)
        print(dfAllResult.iloc[:, 3:].sum().values.tolist())


def Main():
    Summation_all_final_result()
    SummationSubIndel()
    #ConfirmValidation()


Main()
