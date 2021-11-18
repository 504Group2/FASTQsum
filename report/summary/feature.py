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
   
    fig1 =px.density_heatmap(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig2 =px.scatter(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig1.write_html("../density.html")
    fig2.write_html("../scatter.html")
    # write-html.py

    f = open('../lenvsquality.html','w')

    message = """
    <html>
        <head> FASTQsum : Read Length vs. Quality</head>
        <body>
            <p>This is the third sections</p>
        </body>
    </html>"""

    f.write(message)
    f.close()
    # Report summary should say 
    # the read len and quality of box with most count
    # percentile in the selected range

def csvToHtml(csv):

    lenSum(csv)
    scSum(csv)
    scVsLen(csv)

def fqToHtml() :
    csv=fqToCsv()  
    csvToHtml(csv)
