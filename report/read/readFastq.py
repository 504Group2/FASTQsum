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
# Extract,Read,Parse ont.exp2.fastq.gz into list bc01
def getFq() :
    
    bc01=[]
    with gzip.open('../fastq/ont.exp2.fastq.gz','rt') as f: # Decompress ont.exp2.fastq.gz 
        for idx, seq_record in enumerate(SeqIO.parse(f, "fastq")): # Read and parse 
            bc01.append(seq_record)
            if idx == 100: # set how many reads we want
                break
    return bc01
#2.seq_record object to lists             
# Breakdown list bc01 into lists of each column
# List for column Read_id   
def colToList(bc01): 
    idcodel=[]
    for i in bc01:
        idcodel.append(i.id)

# List for column Sequence_length_template
    seql = []
    for i in bc01:
        seql.append(len(i.seq))

# description from fastq is a string
# split string into 2 columns: Start_time, Barcode_arrangement
    descriptionl=[]
    start_timel=[]
    barcodel=[]
    for i in bc01:
        descriptionl.append(i.description)
   
    for des in descriptionl:
        dl=des.split()
        start_timel.append(dl[5][11:])
        barcodel.append(dl[6][8:])

# List for column Mean_qscore_template
    qual = []
    for i in bc01:
        qual.append(np.mean(i.letter_annotations['phred_quality']))
    mycolList=[idcodel,start_timel,seql,qual,barcodel]
    return mycolList


def listToDf(mycolList): 
#3.lists to series 
#Convert lists to pandas Series     
    readid = pd.Series(mycolList[0], name='Read_ID')
    start_times = pd.Series(mycolList[1], name='Start_time')
    seqs = pd.Series(mycolList[2], name='Seq_length_template')
    quality = pd.Series(mycolList[3], name='Mean_qscore_template')
    barcodes = pd.Series(mycolList[4], name='Barcode_arrangement')

#4.series to dataframe 

    fastqdf=pd.DataFrame(dict(Read_ID=readid,Start_time=start_times,Sequence_length_template=seqs,Mean_qscore_template=quality,Barcode_arrangement=barcodes)).set_index(['Read_ID'])
    return fastqdf

#5.convert dataframe to csv file
def dfToCsv(fastqdf):
    csvLocation='../test-1.csv'
    fastqdf.to_csv(csvLocation, index=None)
    #fastqdf.to_csv('/Users/naphat/Desktop/504/test-1.csv', index=None)
    pd.read_csv(csvLocation)
    #success="Create .csv successfully"
    return csvLocation

def fqToCsv():
    bc01=getFq()
    newcolList=colToList(bc01)
    mydf=listToDf(newcolList)
    print(dfToCsv(mydf))
    return dfToCsv(mydf)




