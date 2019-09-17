import os, sys, logging
import pandas as pd
import subprocess as sp
from pdb import set_trace

sOutput_dir = sys.argv[1]
strSample  = sys.argv[2]
strLogPath  = sys.argv[3]

logging.basicConfig(format='%(process)d %(levelname)s %(asctime)s : %(message)s',
                    level=logging.DEBUG,
                    filename=strLogPath,
                    filemode='a')
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def Parsing_summary():

    dfSummary                = pd.read_table('{outdir}/Tmp/{sample}_Summary.txt'.format(sample=strSample, outdir=sOutput_dir), header=None)
    dfSummary.columns        = ['Barcode', 'Total', 'Insertion', 'Deletion', 'Complex']
    dfSummary                = dfSummary.groupby(['Barcode']).sum()
    dfSummary['Total_indel'] = dfSummary['Insertion'] + dfSummary['Deletion'] + dfSummary['Complex']
    dfSummary['IND/TOT']     = dfSummary['Total_indel'] / dfSummary['Total']
    dfSummary['IND/TOT'].fillna(0, inplace=True)
    dfSummary.to_csv('{outdir}/Result/{sample}_Summary_result.tsv'.format(sample=strSample, outdir=sOutput_dir), sep='\t')

def Annotate_final_result():

    dfCount_INDEL = pd.read_table('{outdir}/Tmp/{sample}_Indel_summary.txt'.format(sample=strSample, outdir=sOutput_dir), header=None)
    dfSummary     = pd.read_table('{outdir}/Result/{sample}_Summary_result.tsv'.format(sample=strSample, outdir=sOutput_dir), index_col='Barcode')

    dfCount_INDEL.set_index(0, inplace=True)
    dfConcat_result  = pd.concat([dfCount_INDEL, dfSummary.loc[:,['Total_indel', 'Total', 'IND/TOT']]],axis=1)
    dfConcat_result.dropna(inplace=True)
    dfConcat_result  = dfConcat_result.reset_index()
    dfConcat_result  = dfConcat_result.loc[:,['index','Total_indel', 'Total', 'IND/TOT', 1,2]]
    dfConcat_result.columns = ['Barcode', 'Total_indel', 'Total', 'IND/TOT', 'Match','Info']
    dfConcat_result  = dfConcat_result.round(2)
    dfConcat_result.to_csv('{outdir}/Result/{sample}_Final_indel_result.tsv'.format(sample=strSample, outdir=sOutput_dir), sep='\t', index=False)

if __name__ == '__main__':
    logging.info('Make a summary result.')
    Parsing_summary()
    Annotate_final_result()
    logging.info('The summary result has been completed.\n\n')
