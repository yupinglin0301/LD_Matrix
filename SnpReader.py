import sys, os
import numpy as np


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

    
    filename = arg['--output']
    dirname = os.path.dirname(filename)
    
    print("Reading BED file...")
    bim, fam, bed = read_plink("/Users/yu-pinglin/Desktop/LD_Matrix/data/geno*.bed")
        
    
    print("LD calculation... ")
    pool = multiprocessing.Pool(processes = threadNum) 
    effectiveI = []
    tmpResults = []
    snpInfoList = []
    totalSNPnum =0
    for i in [0]:
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
        tmpResults.append(pool.apply_async(calBlockCorr, args=(blockGenotype)))    
        


   
    

    
    
    
    
   
    
        
    

    
 
    
  
if __name__ == '__main__':
    main_with_args(sys.argv) 
    