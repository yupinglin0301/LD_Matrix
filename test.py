import sys,time,os
import numpy as np
import logger
import multiprocessing
import pandas as pd
from pandas_plink import read_plink


def calBlockCorr(blockGenotype):
    # Standardize
    indNum, blockSNPnum = blockGenotype.shape
    af = np.nanmean(blockGenotype, axis=0)/2.
    expectation = np.outer(np.ones(indNum), 2.*af)
    scale = np.nanstd(blockGenotype, axis=0) #np.sqrt(2*af*(1-af))
    scaleMat = np.outer(np.ones(indNum), scale)
    blockGenotype_norm = (blockGenotype-expectation)/scaleMat
    mask = (~np.isnan(blockGenotype_norm)).astype(int)
    blockGenotype_norm[mask==0] = 0
        
    blockLD = np.around(np.dot(blockGenotype_norm.T, blockGenotype_norm)/indNum, decimals=3)
   
    return blockLD, af



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
         '--output':'results/LD.h5',
         '--log': 'plinkLD.log',
         '--thread':multiprocessing.cpu_count(), 
         '--method':'Pearson', 
         '--compress': 9}
    
    sys.stdout = logger.logger(arg['--log'])
    

    
    # Output blockSNPinfo and blockLD into HDFStore  
    filename = arg['--output']
    dirname = os.path.dirname(filename)
    if not  os.path.exists(dirname): print("123")

if __name__ == '__main__':
    main_with_args(sys.argv) 