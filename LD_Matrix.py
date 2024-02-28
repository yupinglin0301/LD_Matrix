import sys,time,re,os
import numpy as np
import logger
import multiprocessing
from sklearn.covariance import ledoit_wolf
import plink2SnpReader
import bgen_reader
import pandas as pd



def main_with_args(args):
    startTime=time.time()
    
    print('=============plinkLD====================')
    helpText = '''
=================HELP=================
--bfile: Binary data file (Default: test)
--bed: Binary data file (Genotypes; Default: test.bed)
--bim: Binary data file (SNP info; Default: test.bim)
--fam: Binary data file (Individual info; Default: test.fam)
--block: Block file (By default, all SNPs are in one block)
--snplist: SNP list file (By default, all SNP pairs are calculated)
--output: output filename (Default: LD.h5)
--method: Correlation estimation method, including: Pearson, LW (Default:Pearson)
--thread: Thread number for calculation (Default: Total CPU number)
--compress: compression level for output (Default: 9) 
--log: log file (Default: plinkLD.log)
--help: Help
======================================'''

    arg={'--bfile':None,
         '--bed':'test.bed',
         '--bim':'test.bim',
         '--fam':'test.fam',
         '--bgen':None,
        '--block':None,
        '--snplist':None, 
        '--output':'LD.h5',
        '--log': 'plinkLD.log',
        '--thread':multiprocessing.cpu_count(), 
        '--method':'Pearson', 
        '--compress': 9}
    

    
    