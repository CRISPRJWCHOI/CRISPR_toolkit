import os,sys
from pdb import set_trace

import multiprocessing as mp

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import Helper

strUser    = sys.argv[1]
strProject = sys.argv[2]

print('Usage : python ./BaseEdit_input_converter.py user_name project_name')
print('Usage : python ./BaseEdit_input_converter.py JaeWoo Test_samples')


"""
--> Conversion format
Barcode.txt
ACACACACACACAGCTCATA:ACACACACACACAGCTCATA
Reference.txt
ACACACACACACAGCTCATA:TTTGTATACACGCATGTATGCATCCTGCAGGTCTCGCTCTGACATGTGGGAAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA
Query reads
ACACACACACACAGCTCATA.txt

BaseEdit output
Barcode.txt
YSKim_0525+01614_98_repeat1:TATACACGCATGTAT
...
Reference.txt
YSKim_0525+01614_98_repeat1:TTTGTATACACGCATGTAT GCATCCTGCAGGTCTCGCTCTGACATGTGGGAAAGCTTGGCGTAACTAGATCTCTACTCTACCACTTGTACTTCAGCGGTCAGCTTACTCGACTTAA
...
Read
YSKim_0525+01614_98_repeat1.txt
"""

def Convert_Indelsearcher_output(strSampleRefGroup):

    listSampleRefGroup = strSampleRefGroup.replace('\n', '').replace('\r', '').split('\t')

    strSample = listSampleRefGroup[0]
    strRef    = listSampleRefGroup[1]

    print('Processing: %s, %s' % (strSample, strRef))

    strBaseEditRefFolder   = '../Base_edit_2/Input/{user}/Reference/{project}/{ref}'.format(user=strUser,
                                                                                            project=strProject,
                                                                                            ref=strRef)
    strBaseEditQueryFolder = '../Base_edit_2/Input/{user}/Query/{project}/{sample}'.format(user=strUser,
                                                                                           project=strProject,
                                                                                           sample=strSample)
    try:
        Helper.MakeFolderIfNot(strBaseEditRefFolder)
        Helper.MakeFolderIfNot(strBaseEditQueryFolder)
    except OSError as e:
        print(e)
        pass

    ## BaseEdit refer format : filename, barcode, reference
    ReferenceFile_in_IndelSearcher = open('./Input/{user}/Reference/{project}/{ref}/Reference_sequence.txt'.format(user=strUser,
                                                                                                                   project=strProject,
                                                                                                                   ref=strRef))
    BarcodeFile_in_IndelSearcher   = open('./Input/{user}/Reference/{project}/{ref}/Barcode.txt'.format(user=strUser,
                                                                                                        project=strProject,
                                                                                                        ref=strRef))
    BarcodeFile_for_BaseEdit       = open('../Base_edit_2/Input/{user}/Reference/{project}/{ref}/Barcode.txt'.format(user=strUser,
                                                                                                                     project=strProject,
                                                                                                                     ref=strRef), 'w')
    Reference_for_BaseEdit         = open('../Base_edit_2/Input/{user}/Reference/{project}/{ref}/Reference.txt'.format(user=strUser,
                                                                                                                       ref=strRef,
                                                                                                                       project=strProject), 'w') ## conversion target to barcode:refseq

    dictBarcodeSeq = {}

    for strBarcodeIndelSearcher, strReferenceIndelSearcher in zip(BarcodeFile_in_IndelSearcher, ReferenceFile_in_IndelSearcher):

        strBarcodeIndelSearcher   = strBarcodeIndelSearcher.replace('\n', '').strip()
        strReferenceIndelSearcher = strReferenceIndelSearcher.replace('\n', '').strip()

        dictBarcodeSeq[strBarcodeIndelSearcher] = []
        BarcodeFile_for_BaseEdit.write(strBarcodeIndelSearcher + ':' + strBarcodeIndelSearcher + '\n') ## first is filename, second is barcode. BaseEdit barcode format
        Reference_for_BaseEdit.write(strBarcodeIndelSearcher + ':' + strReferenceIndelSearcher + '\n')

    ReferenceFile_in_IndelSearcher.close()
    BarcodeFile_in_IndelSearcher.close()
    Reference_for_BaseEdit.close()

    Total_result_file = open('./Output/{user}/{project}/{sample}/Tmp/{sample}_Classified_Indel_barcode.fastq'.format(user=strUser,
                                                                                                                     project=strProject,
                                                                                                                     sample=strSample))

    intCheckTotLine = 0
    intOneLineMore  = 0

    for i, strRow in enumerate(Total_result_file):  ## for query reads

        if intOneLineMore == 1:
            intCheckTotLine = 0
            intOneLineMore  = 0

        if i % 4 == 0: ## Classified_Indel_barcode has all total sequence.
            strBarcode = strRow.split('Barcode_')[1].split(':')[0]
            intCheckTotLine = 1

        elif intCheckTotLine == 1:
            dictBarcodeSeq[strBarcode].append(strRow)
            intOneLineMore = 1

    for strBarcode, listSeq in dictBarcodeSeq.items():
        with open('../Base_edit_2/Input/{user}/Query/{project}/{sample}/{barcode}.txt'.format(
                user=strUser, project=strProject, sample=strSample, barcode=strBarcode), 'w') as Output:
            Output.write(''.join(listSeq))

    Total_result_file.close()


def Main():
    print('Program Start')
    p = mp.Pool(2)

    with open('./User/{user}/{project}.txt'.format(user=strUser, project=strProject)) as SampleList:
        listSampleRefGroup = [strSampleRefGroup for strSampleRefGroup in SampleList if strSampleRefGroup[0] != '#']
        p.map_async(Convert_Indelsearcher_output, listSampleRefGroup).get()

    p.close()

    print('Program End')


Main()
