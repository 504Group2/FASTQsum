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

def lenSum(csv): #birth
    
    csvOut = pd.read_csv('../test.csv')
    templatelength = []
    templatelength = csvOut['Sequence_length_template']
    templatelength = list(templatelength)
    numberreads = len(csvOut)
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
     
    #print("Number of reads    : ", numberreads)
    #print("Total base        : ", totalbase)
    #print("Mean read length   : ", meanlength)
    #print("Median read length : ", medianlength)
    #print("Read length (N50)  : ", N50size)
    #print("Longest pass read  : ", maxread)
     
    ## --------------------------------------------------------------------------------------
    figTable = go.Figure(data=[go.Table(header=dict(values=['Number of reads','Total bases','Mean read length','Median read length','Read length (N50)','Longest pass read'], fill_color= 'lavender'),
                cells=dict(values=[numberreads, totalbase, meanlength, medianlength, N50size, maxread], fill_color= 'lightcyan'))], 
                layout=go.Layout(title=go.layout.Title(text="Basecall summary")
                , height=240
                ))
    
     
    #figTable.show()
    # figTable.write_image("figTable.png")
    # figTable.write_html('figTable.html')
     
    ## Graph: All barcodes
    ## --------------------------------------------------------------------------------------
    group_label = ['readlength'] 
    
    lenfig = ff.create_distplot([csvOut['Sequence_length_template']], group_label, colors = ['#17202A'], show_hist = False)
    lenfig.update_layout(
        title="Basecalled reads length for all barcodes",
        xaxis_title="Basedcall length",
        yaxis_title="Read density",
    )
    lenfig.add_vline(x=meanlength, line_width=2, line_dash="dash", line_color="Red", annotation_text="Mean", annotation_position="top right")
    lenfig.add_vline(x=medianlength, line_width=2, line_dash="dash", line_color="Blue", annotation_text="Median", annotation_position="top left")
    lenfig.update_xaxes(type='log')
     
    #lenfig.show()
    lenfig.write_image("lenfig.png")
    ## Graph: Split barcode by facet
    ## --------------------------------------------------------------------------------------
    csvOut['log_Sequence_length_template']=np.log10(csvOut['Sequence_length_template'])
    lenfigSplit = px.histogram(data_frame=csvOut, x=csvOut['log_Sequence_length_template'], facet_col='Barcode_arrangement', color='Barcode_arrangement', title="Basecalled reads length for each barcode", log_x=False)
    '''
    lenfigSplit.update_layout(
    xaxis_title="Basedcall length",
    xaxis2=dict(title="Basedcall length"),
    yaxis_title="Read count",
    legend_title="Barcode"
    )
    '''
    lenfigSplit.update_xaxes(title='')
    lenfigSplit.update_layout(xaxis2=dict(title="log10 Basedcall length"), yaxis_title="Read count", legend_title="Barcode")    
    lenfigSplit.for_each_annotation(lambda a: a.update(text=a.text.replace("Barcode_arrangement=", "")))    
  
    #lenfigSplit.show()
    lenfigSplit.write_image("lenfigSplit.png") 
    
     
    ## Write HTML
    html_template = '''<!doctype html>
    <html lang="en">
    <head>
    </head>
    <body>
    {results}
    </body>
    </html>
    '''
    
    figTable_html = '<div><h2>Basecall summary</h2>'+figTable.to_html(full_html=False, include_plotlyjs='cdn')+'<p style="color:SlateGray;"><b>Note: </b>In computational biology, N50 is statistics of a set of contig or scaffold lengths. The N50 is similar to a mean or median of lengths, but has greater weight given to the longer contigs. It is used widely in genome assembly, especially in reference to contig lengths within a draft assembly.</p></div>'
    lenfig_html = '<div><h2>Basecalled reads length</h2><h3 style="color:Navy;"><ins>Basecalled reads length for all barcodes</ins></h3>'+'<p style="text-align:center;"><img src="lenfig.png" width="1000" height="600"></p>'+'<p style="color:DodgerBlue;">Blue line: Median</p><p style="color:Tomato;">Red line: Mean</p><p><strong>Explanation: </strong>Basecalled reads length represents distribution plot (distplot) to show the relationship between read density on y-axis and basecall length as a logarithmic scale on x-axis for all barcodes in FASTQ file.</p><p style="color:SlateGray;"><b>Note: </b><br>- <ins>A distplot or distribution plot</ins> depicts the variation in the data distribution.<br>- <ins>A logarithmic scale (or log scale)</ins> is a way of displaying numerical data over a very wide range of values in a compact way.<br></p></div>'
    lenfigSplit_html = '<div><h3 style="color:Navy;"><ins>Basecalled reads length for each barcode</ins></h3>'+'<p style="text-align:center;"><img src="lenfigSplit.png" width="1000" height="600" class="center"></p>'+'<p><strong>Explanation: </strong>Basecalled reads length represents histogram plot (histplot) to show the relationship between count as number of reads on y-axis and basecall length as a logarithmic scale on x-axis for each barcode in FASTQ file.</p><p style="color:SlateGray;"><b>Note: </b><br>- <ins>A histplot or histogram plot</ins> is an excellent tool for visualizing and understanding the probabilistic distribution of numerical data or image data.<br>- <ins>A logarithmic scale (or log scale)</ins> is a way of displaying numerical data over a very wide range of values in a compact way.<br></p></div>'
     
    results = figTable_html+lenfig_html+lenfigSplit_html
    print("Read Length Summary Complete")
    return results
    
def scSum(csv):
    
    sum = pd.read_csv(csv)
    
    #Summary table
    grouptable = sum.groupby('Barcode_arrangement')
    qsummarytable = grouptable['Mean_qscore_template'].describe()
    qfigTable = ff.create_table(qsummarytable, index='Barcode_arrangement')

    #Basecalled reads PHRED quality (Distribution plot - All reads)
    group_label = ['Mean_qscore_template']
    qDfig = ff.create_distplot([sum['Mean_qscore_template']], group_label, show_hist = False)
    qDfig.update_layout(xaxis_title="Reads quality scores", yaxis_title="Read density")
    qDfig.add_vline(x=8.0, line_width=1.5, line_dash="dash", line_color="red")
    qDfig.write_image("qDfig.png")     

    #Basecalled reads PHRED quality (Histrogram plot - Each of barcode)
    qHfig = px.histogram(sum, x=sum['Mean_qscore_template'], color=sum['Barcode_arrangement'], facet_col=sum['Barcode_arrangement'], marginal="rug")
    qHfig.update_xaxes(title='')
    qHfig.update_layout(xaxis2=dict(title="Reads quality scores"), yaxis_title="Read count", legend_title="Barcode")
    qHfig.for_each_annotation(lambda a: a.update(text=a.text.replace("Barcode_arrangement=", "")))
    qHfig.add_vline(x=8.0, line_width=1.5, line_dash="dash", line_color="red")
    qHfig.write_image("qHfig.png")

    #Number of reads per quality score
    PassReads =  sum.query('Mean_qscore_template >= 8')
    FailReads =  sum.query('Mean_qscore_template < 8')

    qlabels = ['Passed Reads','Failed Reads']
    qvalues = [len(PassReads), len(FailReads)]

    qPiefig = go.Figure(data=[go.Pie(labels=qlabels, values=qvalues, textinfo='label+percent', insidetextorientation='radial', pull=[0, 0.2])])
    qPiefig.update_layout(legend_title="Type of reads")

    #HTML
    qfigTable_html = '<div><h2>Quality score summary</h2>'+qfigTable.to_html(full_html=False, include_plotlyjs='cdn')+'<p><strong>Explanation: </strong>The quality score summary table shows the descriptive statistics information in each barcode arrangement.</p></div>'
    qfig_html = '<div><h2>Basecalled reads PHRED quality</h2><p style="text-align:center;"><img src="qDfig.png" width="1000" height="600"></p><p style="color:Tomato;">Red line: Cut-off line suggestion (Mean quality score at 8.0)</p><p><strong>Explanation: </strong>Basecalled reads PHRED quality distplot represents the distribution of mean quality score for all reads.</p><p style="text-align:center;"><img src="qHfig.png" width="1000" height="600"></p><p style="color:Tomato;">Red line: Cut-off line suggestion (Mean quality score at 8.0)</p><p><strong>Explanation: </strong>Basecalled reads PHRED quality histplot represents the frequency distribution of mean quality score in each barcode arrangement.</p></div>'
    qPiefig_html = '<div><h2>Number of reads per quality score</h2>'+qPiefig.to_html(full_html=False, include_plotlyjs='cdn')+'<p><strong>Explanation: </strong>Number of reads per quality score plot represents the proportion between the number of the passed and failed reads.</p></div>'

    results = qfigTable_html+qfig_html+qPiefig_html
    print("Quality Score Summary Complete")
    return results
   
def scVsLen(csv):
    #Bam
    
    c=pd.read_csv(csv)
    c['log_Sequence_length_template']=np.log10(c['Sequence_length_template'])
    heatlog=px.density_heatmap(data_frame=c,x=c.loc[:,'log_Sequence_length_template'],y=c.loc[:,'Mean_qscore_template'],facet_col='Barcode_arrangement',nbinsx=5,nbinsy=10,histnorm='percent',range_color=[0,100] ,title='Density Heatmap Read Length vs PHRED Score')
    heatlog.update_xaxes(title='') 
    heatlog.update_layout(xaxis2=dict(title="Reads Length"), yaxis_title="Quality Score", legend_title="Barcode")    
    heatlog.for_each_annotation(lambda a: a.update(text=a.text.replace("Barcode_arrangement=", "")))  
    '''
    heatlog.update_layout(
    xaxis_title="Basedcall length",
    xaxis2=dict(title="Basedcall length"),
    yaxis_title="Quality Score",
    legend_title="Barcode"
    )
    '''
    
      
    scatterlog=px.scatter(data_frame=c,x=c.loc[:,'log_Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'],facet_col='Barcode_arrangement',opacity=0.2,title='Scatter plot Read Length vs Quality Score')
    scatterlog.update_xaxes(title='')  
    scatterlog.update_layout(xaxis2=dict(title="Reads Length"), yaxis_title="Quality Score", legend_title="Barcode")    
    scatterlog.for_each_annotation(lambda a: a.update(text=a.text.replace("Barcode_arrangement=", "")))    
  
    '''
    scatterlog.update_layout(
    xaxis_title="Basedcall length",
    xaxis2=dict(title="Basedcall length"),
    yaxis_title="Quality Score",
    legend_title="Barcode"
    )
    '''
    scatterlog.write_image("scatterlog.png")
    heatdiv='''
                <div>
                    <p><strong>Explanation: </strong> The heatmaps share the same x-y axis with scatter plot. These colorful interactive plot illustrate how many reads are in each range of read length and quality score. Heatmap draws a grid that count reads in each partition. To inteprete, yellow color represent high frequency (many reads). In contrast, purple color means the lower percentage. 
                    The percentage of each shade is as displayed in the color bar.</p></div>
                </div>'''
    scatterdiv='''
                <div>
                <p> This last section illustrates the relationship between read length and quality score. In other word, you could conveniently see how many good read there are at each length.</p>
                <h2> Scatter plot </h2>
                    <p style="text-align:center;"><img src="scatterlog.png" /></p>
                    <p><strong>Explanation: </strong>Each dot represents each read in the fastq. The X-axis is the log10-transformed read length, while the y-axis is the quality score. Reads with different barcode are plot separately. Scatter plot provide a simple explicit visualization of how length and quality distributed across dataset. Each dot is slightly opage. Thus you can observe where the dot overlap and condense.</p></div>
                </div>'''
    results = '<div><h2>Read Length vs Quality Score Summary</h2>'+scatterdiv+'<h2> Heat Map </h2>'+heatlog.to_html(full_html=False, include_plotlyjs='cdn')+heatdiv+'</div>'
    print("Length vs Quality Complete")

    return results
    
    


def csvToHtml(csv):
    

    html_template = '''<!doctype html>
    <html>
    <head>
    <p><img src="flogo.png" width="200"/></p>    
    </head>

    {birth}
    
    {pe}
    
    {bam}
    

    </html>
    '''

    with open('report.html','w') as outf:
        outf.write(html_template.format(birth=lenSum(csv),pe=scSum(csv),bam=scVsLen(csv)))

def fqToHtml(filePath) :
    csv=fqToCsv(filePath)  #csvlocation
    csvToHtml(csv)
