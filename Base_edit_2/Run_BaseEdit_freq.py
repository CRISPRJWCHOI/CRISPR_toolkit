import os, re, sys, pdb, math, logging

import subprocess as sp

from pdb import set_trace
from datetime import datetime
from optparse import OptionParser

sys.path.insert(0, os.path.dirname(os.getcwd()))
from Core.CoreSystem import InitialFolder, UserFolderAdmin, Helper, RunMulticore, CheckProcessedFiles


class clsBaseEditRunner(UserFolderAdmin):

    def __init__(self, strSample, strRef, options, InstInitFolder):
        UserFolderAdmin.__init__(self, strSample, strRef, options, InstInitFolder.strLogPath)
        self.MakeSampleFolder()

        self.strSample        = strSample
        self.strRef           = strRef
        self.intCore          = options.multicore
        self.strGapOpen       = options.gap_open
        self.strGapExtend     = options.gap_extend
        self.strTargetWindow  = options.target_window
        self.strIndelCheckPos = options.indel_check_pos
        self.strTargetRefAlt  = options.target_ref_alt

        self.strBarcodeFile      = os.path.join(self.strRefDir, 'Barcode.txt')
        self.strReferenceSeqFile = os.path.join(self.strRefDir, 'Reference.txt')
        self.strRefFile          = os.path.join(self.strRefDir, 'Reference.fa')

        self.strPamSeq   = options.PAM_seq
        self.strPamPos   = options.PAM_pos
        self.strGuidePos = options.Guide_pos

        Helper.MakeFolderIfNot('./Output/{user}/{project}/{sample}/Tmp/Alignment'.format(user=self.strUser,
                                                                                         project=self.strProject,
                                                                                         sample=self.strSample))

    def MakeReference(self):

        with open(self.strBarcodeFile) as Barcode, \
            open(self.strReferenceSeqFile) as Ref, \
            open(self.strRefFile, 'w') as Output:

            listBarcode = Helper.RemoveNullAndBadKeyword(Barcode)
            listRef     = Helper.RemoveNullAndBadKeyword(Ref)

            ## defensive
            assert len(listBarcode) == len(listRef), 'Barcode and Reference must be a same row number.'

            dictBarcode = {}

            for strBarcode in listBarcode:
                strBarcode = strBarcode.replace('\n','').replace('\r','')
                Helper.CheckIntegrity(self.strBarcodeFile, strBarcode) ## defensive
                listBarcode   = strBarcode.split(':')
                strBarSample  = listBarcode[0]
                strBarcode    = listBarcode[1]
                dictBarcode[strBarSample] = strBarcode

            for strRef in listRef:
                strRef = strRef.replace('\n','').replace('\r','')
                Helper.CheckIntegrity(self.strBarcodeFile, strRef) ## defensive
                listRef      = strRef.split(':')
                strRefSample = listRef[0]
                strRef       = listRef[1]

                try:
                    sBarcode = dictBarcode[strRefSample]
                    Output.write('%s\t%s\t%s\n' % (strRefSample, sBarcode, strRef))
                except KeyError:
                    logging.error('no matching')
                    logging.error(strRefSample,strRef)

    def MakeIndelSearcherCmd(self):

        listCmd = []

        with open(self.strRefFile) as BarcodeRef:

            for strBarcodeRef in BarcodeRef:
                listBarcodeRef = strBarcodeRef.replace('\n', '').replace('\r','').split('\t')
                strFileName    = listBarcodeRef[0]
                strBarcode     = listBarcodeRef[1]
                strRef         = listBarcodeRef[2]

                strForwardQueryFile = './Input/{user}/Query/{project}/{sample}/{file_name}.txt'.format (user=self.strUser,
                                                                                                       project=self.strProject,
                                                                                                       sample=self.strSample,
                                                                                                       file_name=strFileName)
                strCmd = ('{python} ./BaseEdit_freq_crispresso.py {forw} {GapO} {GapE} {barcode} {ref} {target_window} {indel_check_pos}'
                          ' {target_ref_alt} {outdir} {file_name} {PAM_seq} {PAM_pos} {guide_pos} {ednafull} {log}').format(
                        python=self.strPython, forw=strForwardQueryFile, GapO=self.strGapOpen, GapE=self.strGapExtend,
                        barcode=strBarcode, ref=strRef, target_window=self.strTargetWindow, indel_check_pos=self.strIndelCheckPos,
                        target_ref_alt=self.strTargetRefAlt, outdir=self.strOutSampleDir, file_name=strFileName,
                        PAM_seq=self.strPamSeq, PAM_pos=self.strPamPos, guide_pos=self.strGuidePos, ednafull=self.strEdnafull, log=self.strLogPath)
                listCmd.append(strCmd)

        return listCmd

    def MakeMergeTarget(self):
        strCmd = '{python} ./Summary_all_trim.py {output} {sample} {ref_alt}'.format(python=self.strPython, output=self.strOutSampleDir,
                                                                                     sample=self.strSample, ref_alt=self.strTargetRefAlt)
        sp.call(strCmd, shell=True)


def Main():
    print 'BaseEdit program start: %s' % datetime.now()

    sCmd = "BaseEdit frequency analyzer\n./Run_BaseEdit_freq.py -t 15 -w 16-48 --indel_check_pos 39-40 --target_ref_alt A,T --PAM_seq NGG --PAM_pos 43-45 --Guide_pos 23-42"
    sCmd += " --gap_open 20 --gap_extend 1 --end_open 20 --end_extend 1\n\n"
    sCmd += "1: Barcode\n"
    sCmd += "2: Base target window\n"
    sCmd += "3: Indel check pos\n"
    sCmd += "4: PAM pos\n"
    sCmd += "5: Guide pos (without PAM)\n"
    sCmd += "TATCTCTATCAGCACACAAGCATGCAATCACCTTGGGTCCAAAGGTCC\n"
    sCmd += "<------1------><----------------2--------------->\n"
    sCmd += "                                     <3>  <4>   \n"
    sCmd += "                      <---------5-------->      \n\n"

    parser = OptionParser(sCmd)

    parser.add_option("-t", "--thread", default="1", type="int", dest="multicore", help="multiprocessing number")
    parser.add_option('--gap_open', default='-10', type='float', dest='gap_open', help='gap open: -100~0')
    parser.add_option('--gap_extend', default='1', type='float', dest='gap_extend', help='gap extend: 1~100')
    parser.add_option("-w", "--target_window", type="str", dest="target_window", help="a window size for target sequence : 20-48")
    parser.add_option("--indel_check_pos", type="str", dest="indel_check_pos", help="indel check position to filter : 39-40; insertion 39, deletion 39 & 40")
    parser.add_option("--target_ref_alt", type="str", dest="target_ref_alt", help="Ref 'A' is changed to Alt 'T': A,T")
    parser.add_option("--PAM_seq", type="str", dest="PAM_seq", help="PAM sequence: NGG, NGC ...")
    parser.add_option("--PAM_pos", type="str", dest="PAM_pos", help="PAM position range in the reference seqeunce : 43-45")
    parser.add_option("--Guide_pos", type="str", dest="Guide_pos", help="Guide position range in the reference seqeunce : 23-42")
    parser.add_option('--python', dest='python', help='The python path including the CRISPResso2')
    parser.add_option('--user', dest='user_name', help='The user name with no space')
    parser.add_option('--project', dest='project_name', help='The project name with no space')
    parser.add_option('--ednafull', dest='ednafull', help='The nucleotide alignment matrix')

    options, args = parser.parse_args()

    InstInitFolder = InitialFolder(options.user_name, options.project_name, os.path.basename(__file__))
    InstInitFolder.MakeDefaultFolder()
    InstInitFolder.MakeInputFolder()
    InstInitFolder.MakeOutputFolder()

    logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                        level=logging.DEBUG,
                        filename=InstInitFolder.strLogPath,
                        filemode='a')
    logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

    logging.info('Program start')
    if options.multicore > 15:
        logging.warning('Optimal treads <= 15')
    logging.info(str(options))

    with open(InstInitFolder.strProjectFile) as Sample_list:

        listSamples        = Helper.RemoveNullAndBadKeyword(Sample_list)

        strInputProject = './Input/{user}/Query/{project}'.format(user=options.user_name, project=options.project_name)

        @CheckProcessedFiles
        def RunPipeline(**kwargs):
            for strSample in listSamples:
                if strSample[0] == '#': continue

                tupSampleInfo = Helper.SplitSampleInfo(strSample)
                if not tupSampleInfo: continue
                strSample, strRef, strExpCtrl = tupSampleInfo

                InstBaseEdit = clsBaseEditRunner(strSample, strRef, options, InstInitFolder)
                InstBaseEdit.MakeReference()

                listCmd = InstBaseEdit.MakeIndelSearcherCmd()
                ###print(lCmd[:5])
                RunMulticore(listCmd, options.multicore)  ## from CoreSystem.py

                InstBaseEdit.MakeMergeTarget()

        RunPipeline(InstInitFolder=InstInitFolder,
                    strInputProject=strInputProject,
                    listSamples=listSamples,
                    logging=logging)

    print 'BaseEdit program end: %s' % datetime.now()


if __name__ == '__main__':
    Main()
