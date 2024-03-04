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
    print('=============plinkLD====================')
    print(helpText)
    
    cpuNum = multiprocessing.cpu_count()
    try:
        threadNum = round(float(arg['--thread']))
    except ValueError:
        print('Warning: --thread must be a numeric value')
        threadNum = cpuNum
    
    if threadNum<=0: 
        threadNum = cpuNum
    print(cpuNum, 'CPUs detected, using', threadNum, 'thread...' )
  

    try:
        complevel = round(float(arg['--compress']))
    except ValueError:
        print('Warning: --compress must be a numeric value')
        complevel = 0
    
   
    print('Start loading SNP information...')
    if arg['--bfile']!=None:
        arg['--bed']=arg['--bfile']+'.bed'
        arg['--bim']=arg['--bfile']+'.bim'
        arg['--fam']=arg['--bfile']+'.fam'
    
    try:
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
    
    
    print("Start building block file... ")
    blockNum = len(chSet)
    blockCH = ['']*blockNum
    blockStart = [0]*blockNum
    blockStop = [0]*blockNum
    print(blockNum, 'Blocks (Each chromosome is an unique block).')
    
    for i in range(0,blockNum):
        blockCH[i] = chSet[i]
        tmpBP = [bp[j] for j in range(0, snpNum) if ch[j]==blockCH[i]]
        blockStart[i] = min(tmpBP)
        blockStop[i] = max(tmpBP)

    
    
    print("Reading BED file...")
    bim, fam, bed = read_plink("data/geno*.bed")
        
    print("LD calculation... ")
    pool = multiprocessing.Pool(processes = threadNum) 
    effectiveI = []
    tmpResults = []
    snpInfoList = []
    totalSNPnum =0
    for i in range(blockNum):
        #SNP index in the block
        if i==0 or not blockCH[i] in blockCH[0:i]:
            idx = [j for j in range(0, snpNum) if ch[j]==blockCH[i] \
                        and blockStart[i]<= bp[j] and bp[j]<=blockStop[i]]
          
        else:
            idx = [j for j in range(0, snpNum) if ch[j]==blockCH[i] \
                        and blockStart[i]< bp[j] and bp[j]<=blockStop[i]]
          
    
        if len(idx)==0: continue
        if len(idx)==1 : 
            print('Block '+ str(i)+ ' : Only '+str(len(idx))+ ' SNP in the block [Ignored]') 
            continue
        print('Block '+ str(i)+ ' : '+str(len(idx))+ ' SNPs [Calculating LD ...]')
        totalSNPnum += len(idx)
        
        blockGenotype = np.transpose(bed.compute()[idx, :])
        blockSNPinfo = snpInfo.iloc[idx]
        
        effectiveI.append(i)
        snpInfoList.append(blockSNPinfo)
        tmpResults.append(pool.apply_async(calBlockCorr, args=(blockGenotype,)))    
    
    # Output blockSNPinfo and blockLD into HDFStore  
    filename = arg['--output']
    dirname = os.path.dirname(filename)
    if not (dirname=='' or os.path.exists(dirname)): os.makedirs(dirname)
    try:
        store = pd.HDFStore(filename, 'w', complevel=complevel)
    except:
        print('Unable to create file', filename)
        exit()
        
    for k in range(len(effectiveI)):
        i = effectiveI[k]
        blockSNPinfo = snpInfoList[k]
        blockLD, af = tmpResults[k].get()
        print('Block '+ str(i)+ ' : '+str(blockLD.shape[0])+ ' SNPs [Finished]')
        blockSNPinfo.insert(6,'F',af)
        blockLD = pd.DataFrame(data=blockLD)
        # Output Large file
        store.put('SNPINFO'+str(i),value=blockSNPinfo, complevel=complevel, format='fixed')
        store.put('LD'+str(i), value=blockLD, complevel=complevel, format='fixed')
    
    pool.close()
    pool.join()
    store.close()
    totalTime = time.time() - startTime
    print("==========================")
    print('Finish Calculation, Total_time =', float(totalTime)/60, 'min.') 

if __name__ == '__main__':
    main_with_args(sys.argv) 
    