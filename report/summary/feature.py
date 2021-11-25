from report.read.readFastq import fqToCsv
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
#16nov
import gzip
from os import name
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
import csv
from scipy import stats
from math import fabs, log
from Bio import SeqIO

def numberReads(csvOut):
    numberreads = len(csvOut)
    return numberreads

## Calculate Total bases
def totalBase(templatelength):
    totalbase = sum(templatelength)
    return totalbase

## Calculate Mean read length
def meanLength(templatelength):
    meanlength = statistics.mean(templatelength)
    return meanlength

## Calculate Median read length
def medianLength(templatelength):
    medianlength = statistics.median(templatelength)    
    return medianlength

def longestread(templatelength):
    lengthsort = templatelength
    lengthsort.sort()
    maxread = max(lengthsort)
    return maxread
## Calculate Read length (N50)



def lenSum(csv): #birth
    print("This is read length summary")
    csvOut = pd.read_csv(csv)
    templatelength = []
    templatelength = csvOut['Sequence_length_template']
    templatelength = list(templatelength)
    totalbase = sum(templatelength)
    halflength = totalbase/2
    lengthsort = templatelength
    lengthsort.sort(reverse=True)
    total = 0
    for i in lengthsort:
        total += i
        if total >= halflength:
            break

    N50size = i
    print("Number of reads    : ", numberReads(csvOut))
    print("Totoal base        : ", totalBase(templatelength))
    print("Mean read length   : ", meanLength(templatelength))
    print("Median read length : ", medianLength(templatelength))
    print("Read length (N50)  : ", N50size)
    print("Longest pass read  : ", longestread(templatelength))
    ## Write basecalled read length summary to CSV
    
    totalbase = sum(templatelength)
    meanlength = statistics.mean(templatelength)
    medianlength = statistics.median(templatelength)
    halflength = totalbase/2
    lengthsort = templatelength
    lengthsort.sort(reverse=True)
    total = 0
    
    for i in lengthsort:
        total += i
        if total >= halflength:
            break
    N50size = i
    lengthsort = templatelength
    lengthsort.sort()
    maxread = max(lengthsort)

    #header = ['Number of reads','Total bases','Mean read length','Median read length','Read length (N50)','Longest pass read']
    #data = [numberReads(), totalbase, meanlength, medianlength, N50size, maxread]
    #with open('readlength.csv','w') as lengthOut:
    #    writer = csv.writer(lengthOut)
    #    writer.writerow(header)
    #    writer.writerow(data)
"""""    
def scSum(csv):
    print("This is score summary")
    sum = pd.read_csv(csv)
    
    #Summary table
    grouptable = sum.groupby('Barcode_arrangement')
    qsummarytable = grouptable['Mean_qscore_template'].describe()
    print(qsummarytable)
    
    #Basecalled reads PHRED quality
    qfig = ff.create_distplot([sum[sum['Barcode_arrangement']=='barcode02']['Mean_qscore_template'],
                         sum[sum['Barcode_arrangement']=='barcode01']['Mean_qscore_template']], 
                         ['barcode02', 'barcode01'],
                         colors = ['#F66095', '#2BCDC1'],
                         show_hist=False
                        )
    qfig.update_layout(title="Basecalled reads PHRED quality", xaxis_title="Reads quality scores", yaxis_title="Read density", legend_title="Barcode")
    qfig.add_vline(x=8.0, line_width=1.5, line_dash="dash", line_color="red", annotation_text="Cut-off line", annotation_font_color="red")

    #Number of reads per quality score
    PassReads =  sum.query('Mean_qscore_template >= 8')
    FailReads =  sum.query('Mean_qscore_template < 8')

    qlabels = ['Pass Reads','Fail Reads']
    qvalues = [len(PassReads), len(FailReads)]

    qPiefig = go.Figure(data=[go.Pie(labels=qlabels, values=qvalues, textinfo='label+percent',
                             insidetextorientation='radial', pull=[0, 0.2]
                            )])
    qPiefig.update_layout(title="Number of reads per quality score", legend_title="Type of reads")
"""    
def scVsLen(csv):
    #Bam
    print("This is score vs summary summary")
    c=pd.read_csv(csv)
    #c1=barcode
    #c2=barcode
    fig1 =px.density_heatmap(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig2 =px.scatter(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig1.write_html("../density.html")
    fig2.write_html("../scatter.html")
    #annotated
    # write-html.py
    # how to combine fig?
    #f = open('../lenvsquality.html','w')

    #message = """
    #<html>
    #    <head> FASTQsum : Read Length vs. Quality</head>
    #    <body>
    #        <p>This is the third sections</p>
    #    </body>
    #</html>"""

    #f.write(message)
    #f.close()
    
    # Report summary should say 
    # the read len and quality of box with most count
    # percentage in the each box

def csvToHtml(csv):

    lenSum(csv) #birth
    #scSum(csv) #pe
    scVsLen(csv) #bam

def fqToHtml(filePath) :
    csv=fqToCsv(filePath)  #csvlocation
    csvToHtml(csv)
