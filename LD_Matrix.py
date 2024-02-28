import sys,time,re,os
import numpy as np
import logger
import multiprocessing
from sklearn.covariance import ledoit_wolf
import pandas as pd


# python plinkLD.py --bfile 1000G.EUR.QC --block fourier_ls-all.bed --out LD.h5



def main_with_args(args):
    startTime=time.time()
    
    print('=============plinkLD====================')
    helpText = '''
=================HELP=================
--bfile: Binary data file (Default: test)
--bed: Binary data file (Genotypes; Default: test.bed)
--bim: Binary data file (SNP info; Default: test.bim)
--fam: Binary data file (Individual info; Default: test.fam)
--output: output filename (Default: LD.h5)
--method: Correlation estimation method, including: Pearson, LW (Default:Pearson)
--thread: Thread number for calculation (Default: Total CPU number)
--compress: compression level for output (Default: 9) 
--log: log file (Default: plinkLD.log)
--help: Help
======================================'''

    arg={'--bfile':None,
         '--bed':'data/geno.bed',
         '--bim':'data/geno.bim',
         '--fam':'data/geno.fam',
        '--output':'LD.h5',
        '--log': 'plinkLD.log',
        '--thread':multiprocessing.cpu_count(), 
        '--method':'Pearson', 
        '--compress': 9}
    
    sys.stdout = logger.logger(arg['--log'])
    
    cpuNum = multiprocessing.cpu_count()
    try:
        threadNum = round(float(arg['--thread']))
    except ValueError:
        print('Warning: --thread must be a numeric value')
        threadNum = cpuNum
    
    if threadNum<=0: 
        threadNum=cpuNum
    print(cpuNum, 'CPUs detected, using', threadNum, 'thread...' )
  
    print(arg['--method'],'method for LD calculation')
    
    
    
    if arg['--bfile']!=None:
        arg['--bed']=arg['--bfile']+'.bed'
        arg['--bim']=arg['--bfile']+'.bim'
        arg['--fam']=arg['--bfile']+'.fam'
    
    print('Start loading SNP information...')
    try:
        #with open(arg['--bim'],'r') as f:
        #snpInfo = [i.strip() for i in f.readlines() if len(i.strip())!=0]
        snpInfo = pd.read_table(arg['--bim'], sep='\s+', names=['CHR','SNP','GD','BP','A1','A2'], dtype={'SNP':str,'CHR':str,'A1':str,'A2':str})
    except:
        print("Could not read SNP Info file:", arg['--bim'])
        exit()
    

    snpNum = len(snpInfo)
    print(snpNum,'SNPs.')

    #Change A1 and A2 to lower case for easily matching
    snpInfo['SNP'] = snpInfo['SNP'].str.lower()
    snpInfo['A1'] = snpInfo['A1'].str.lower()
    snpInfo['A2'] = snpInfo['A2'].str.lower()
    
    ch = snpInfo['CHR'].tolist()#[0]*snpNum  #Chromosome
    snpID = snpInfo['SNP'].tolist()#['']*snpNum #SNP ID
    bp = snpInfo['BP'].tolist()#[0]*snpNum  #Base position
    chSet = list(set(ch))
    
    
    print("Start loading individual information... ")
    try:
        with open(arg['--fam'],'r') as f:
            indInfo = [i.strip() for i in f.readlines()]  
    except IOError:
        print("Could not read Individual Info file", arg['--fam'])
        exit()
    
    indNum = len(indInfo)
    print(indNum,'individuals.')
    
    
  
if __name__ == '__main__':
    main_with_args(sys.argv) 
    