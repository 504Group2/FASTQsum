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
bc01=[]
with gzip.open('fastq/ont.exp2.fastq.gz','rt') as f: # Decompress ont.exp2.fastq.gz 
    for idx, seq_record in enumerate(SeqIO.parse(f, "fastq")): # Read and parse 
        bc01.append(seq_record)
        if idx == 100: # set how many reads we want
            break

#2.seq_record object to lists             
# Breakdown list bc01 into lists of each column
# List for column Read_id   
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
    
#3.lists to series 
# Convert lists to pandas Series  
readid = pd.Series(idcodel, name='Read_ID')
start_times = pd.Series(start_timel, name='Start_time')
seqs = pd.Series(seql, name='Seq_length_template')
quality = pd.Series(qual, name='Mean_qscore_template')
barcodes = pd.Series(barcodel, name='Barcode_arrangement')

#4.series to dataframe 
fastqdf=pd.DataFrame(dict(Read_ID=readid,Start_time=start_times,Sequence_length_template=seqs,Mean_qscore_template=quality,Barcode_arrangement=barcodes)).set_index(['Read_ID'])
print(fastqdf)

#5.convert dataframe to csv file
fastqdf.to_csv('test-1.csv', index=None)
pd.read_csv('test-1.csv') 
