#!/home/hkimlab/anaconda2/bin/python2.7



import os

"""
    with open(r"E:\ABMRC_document\Kim_hui_kwan\Indel_searcher\IndelSearcher_ver2\Indel_20161027(must_detect)\Make_ref_with_barcode_and_target\Barcode.txt") as Barcode,\
        open(r"E:\ABMRC_document\Kim_hui_kwan\Indel_searcher\IndelSearcher_ver2\Indel_20161027(must_detect)\Make_ref_with_barcode_and_target\target_region.txt") as Target,\
        open(r"E:\ABMRC_document\Kim_hui_kwan\Indel_searcher\IndelSearcher_ver2\Indel_20161027(must_detect)\Make_ref_with_barcode_and_target\Reference sequence.txt") as Ref,\
        open('./Consensus.fa', 'w') as Output:
"""

def Make_ref():

    with open(r"./Input/Reference/Barcode.txt") as Barcode,\
        open(r"./Input/Reference/Target_region.txt") as Target,\
        open(r"./Input/Reference/Reference_sequence.txt") as Ref,\
        open('./Input/Reference/Consensus.fa', 'w') as Output:

        lName = [sBar.replace('\n', '') + ':' + sTar for sBar, sTar in zip(Barcode, Target)]

        for i, sRow in enumerate(Ref):
            Output.write('>'+lName[i]+sRow)

Make_ref()
