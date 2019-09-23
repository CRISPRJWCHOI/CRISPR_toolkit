#!/media/hkim/7434A5B334A5792E/bin/Python/Python2/bin/python2

import os,sys
import numpy as np
from pdb import set_trace


sProject = sys.argv[1]


def Sum_all_freq():

    sFile_path = './Output/%s/Summary/All' % sProject
    sHeader    = ''

    """
    Sample	         Barcode	     Ref	                 # of Total	# of Insertion	# of Deletion	# of Combination	T.-7	A.-6  G.-5
    Doench2014_1001	ATACATAGCTACATG	CAGCGGTCAGCTTACTCGACTTAA... 	60	         0	          0	               0	          0	     0     	0
    																														  0	     0	    0 
    																														  0	     0	    0
    																														  0	     0	    0
    """

    lSum_total_and_indel_data = []
    lSum_target_data = []

    for iFile_num, sFile in enumerate(os.listdir(sFile_path)):
        #print(iFile_num)
        with open(sFile_path + '/' + sFile) as Input:
            lSum_target = []

            for i, sRow in enumerate(Input):
                if i == 0:
                    sHeader = sRow
                    continue

                lCol = sRow.replace('\n','').split('\t')

                if i == 1: ## This data is in the second row
                    lTotal_and_indel_col = map(int, lCol[3:7])
                    if lSum_total_and_indel_data == []:
                        lSum_total_and_indel_data = np.zeros((len(lTotal_and_indel_col)), int)
                    lSum_total_and_indel_data += lTotal_and_indel_col

                lTarget_col = map(int, lCol[7:])
                if lSum_target_data == []:
                    lSum_target_data = np.zeros((4, len(lTarget_col)), int)

                lSum_target.append(lTarget_col)

            if lSum_target:
                lSum_target_data += lSum_target

    print(lSum_target_data)

    with open('./Output/%s/Summary/Alt_freq.txt' % sProject, 'w') as Output:

        lHeader = sHeader.split('\t')
        lHeader[7:] = [sCol.split('.')[1] for sCol in lHeader[7:]]
        Output.write('Alt_base\t' + '\t'.join(lHeader[3:]))
        
        cnt = -1

        for sBase, lSum in zip(['A','C','G','T'], lSum_target_data):
            cnt += 1
            if cnt == 0:
                Output.write(sBase + '\t' + '\t'.join(map(str, lSum_total_and_indel_data)) + '\t' + '\t'.join(map(str, lSum)) + '\n')
            else:
                Output.write(sBase + '\t\t\t\t\t' + '\t'.join(map(str, lSum)) + '\n')


def Main():
    Sum_all_freq()


Main()
