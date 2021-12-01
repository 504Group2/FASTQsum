import gzip
import numpy as np
import pandas as pd
from Bio import SeqIO

# Goal : Convert input .fastq into dataframe
# Steps : 
#1.Read fastq 
#2.seq_record object to lists 
#3.lists to series 
#4.series to dataframe  

#1.Read fastq 
# Extract,Read,Parse ont.exp2.fastq.gz into list seqList
def getFq(filePath) :
    
    seqList=[]
    with gzip.open(filePath,'rt') as f:
    #with gzip.open('../fastq/ont.exp2.fastq.gz','rt') as f: # Decompress ont.exp2.fastq.gz 
        for idx, seq_record in enumerate(SeqIO.parse(f, "fastq")): # Read and parse 
            seqList.append(seq_record)
            if idx == 1000: # set how many reads we want
                break
    return seqList
#2.seq_record object to lists             
# Breakdown list seqList into lists of each column
# List for column Read_id   
def colToList(seqList): 
    index=range(len(seqList))
    idcodel=[]
    for i in seqList:
        idcodel.append(i.id)

# List for column Sequence_length_template
    seql = []
    for i in seqList:
        seql.append(len(i.seq))

# description from fastq is a string
# split string into 2 columns: Start_time, Barcode_arrangement
    descriptionl=[]
    
    barcodel=[]
    for i in seqList:
        descriptionl.append(i.description)
   
    for des in descriptionl:
        dl=des.split()
       
        barcodel.append(dl[6][8:])

# List for column Mean_qscore_template
    qual = []
    for i in seqList:
        qual.append(np.mean(i.letter_annotations['phred_quality']))
        mycolList=[index,idcodel,seql,qual,barcodel]
    return mycolList


def listToDf(mycolList): 
#3.lists to series 
#Convert lists to pandas Series  
#   
    indexse = pd.Series(mycolList[0], name='index')
    readid = pd.Series(mycolList[1], name='Read_ID')
    
    seqs = pd.Series(mycolList[3], name='Sequence_length_template')
    quality = pd.Series(mycolList[4], name='Mean_qscore_template')
    barcodes = pd.Series(mycolList[5], name='Barcode_arrangement')

#4.series to dataframe 

    fastqdf=pd.DataFrame(dict(Index=indexse,Read_ID=readid,Sequence_length_template=seqs,Mean_qscore_template=quality,Barcode_arrangement=barcodes)).set_index(['Index'])
    return fastqdf

#5.convert dataframe to csv file
def dfToCsv(fastqdf):
    csvLocation='../test.csv'
    fastqdf.to_csv(csvLocation, index=None) #Read_id
    #fastqdf.to_csv('/Users/naphat/Desktop/504/test-1.csv', index=None)
    pd.read_csv(csvLocation)
    #success="Create .csv successfully"
    return csvLocation

def fqToCsv(filePath):
    seqList=getFq(filePath)
    newcolList=colToList(seqList)
    mydf=listToDf(newcolList)
    print(dfToCsv(mydf))
    return dfToCsv(mydf) #.csv location




