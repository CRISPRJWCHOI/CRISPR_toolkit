#!/home/hkimlab/anaconda2/bin/python2.7

import os, sys
import subprocess as sp
from pdb import set_trace

print('Additional_BaseEdit_process_list.tsv')
print('format example: ABE_Rep1_1_2_TF4\tA,T\tABE_Rep1_1_2_TF4_CtoA_Summary.txt')
print('table format  :    project      ref,alt        first merged result      ')


with open('Additional_BaseEdit_process_list.tsv') as Input:
	for sRow in Input:
		if sRow[0] == '#': continue
		lCol = sRow.replace('\n', '').replace('\r','').split('\t')
		if len(lCol) == 1:
			lCol = lCol[0].split()
		print(lCol)
		sCmd = './Each_base_summary.py {project} {ref_alt} {first}'.format(project=lCol[0], ref_alt=lCol[1].replace('"',''), first=lCol[2])
		print(repr(sCmd))
		sp.call((sCmd), shell=True)
