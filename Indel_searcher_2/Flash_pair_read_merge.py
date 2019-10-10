import os, sys
import subprocess as sp

strUser    = sys.argv[1]
strProject = sys.argv[2]
strFlash   = sys.argv[3]
strThread  = sys.argv[4]


def RunFlash():

    strFlashDir   = '../{flash}'.format(flash=strFlash)
    strProjectDir = './Input/{user}/FASTQ/{project}'.format(user=strUser, project=strProject)

    for strSampleDir in os.listdir(strProjectDir):
        strSamplePath = os.path.join(strProjectDir, strSampleDir)

        if os.path.isdir(strSamplePath):

            listPairFiles = []

            for strFile in os.listdir(os.path.join(strProjectDir, strSampleDir)):
                if '_1.fastq.gz' in strFile or '_2.fastq.gz' in strFile:
                    listPairFiles.append(strFile)

            strForward = os.path.join(strSamplePath, listPairFiles[0])
            strReverse = os.path.join(strSamplePath, listPairFiles[1])
            strOutput  = os.path.join(strSamplePath, listPairFiles[0].replace('_1.fastq.gz', ''))

            strLog = './Output/{user}/{project}/Log'.format(user=strUser,
                                                            project=strProject)

            if not os.path.isdir(strLog): os.makedirs(strLog)

            strCmd = '{flash_dir}/flash -m 10 -M 400 -O -o {output} -t {thread} {r1} {r2} >{log}/flash.log 2>&1 '.format(
                flash_dir=strFlashDir,
                output=strOutput,
                thread=strThread,
                r1=strForward,
                r2=strReverse,
                log=strLog)

            print(strCmd)
            sp.call(strCmd, shell=True)
            print('complete, {fow} {rev} are moved to project folder'.format(fow=listPairFiles[0], rev=listPairFiles[1]))
            sp.call('mv {sample_path}/*.fastq.gz {project_dir} &&'
                    ' rm {sample_path}/*hist* {project_dir} &&'
                    ' rm {sample_path}/*notCombined* {project_dir}'.format(sample_path=strSamplePath,
                                                                           project_dir=strProjectDir), shell=True)


def Main():
    RunFlash()


Main()
