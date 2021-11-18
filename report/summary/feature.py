from report.read.readFastq import fqToCsv
import pandas as pd
import plotly.express as px
#16nov

def lenSum(csv):
    print("This is read length summary")
def scSum(csv):
    print("This is score summary")
def scVsLen(csv):
    #Bam
    print("This is score vs summary summary")
    c=pd.read_csv(csv)
    #print(c['Read_ID'])
    #print (pd.DataFrame.head(c))
    #print(c.loc[:,'Sequence_length_template'])
    fig1 =px.density_heatmap(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    #fig2 =px.scatter(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig1.write_html("../file.html")
    #fig2.write_html("../file.html")

def csvToHtml(csv):

    lenSum(csv)
    scSum(csv)
    scVsLen(csv)

def fqToHtml() :
    csv=fqToCsv()  
    csvToHtml(csv)
