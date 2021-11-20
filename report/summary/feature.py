from report.read.readFastq import fqToCsv
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
#16nov

def lenSum(csv):
    print("This is read length summary")
def scSum(csv):
    print("This is score summary")
    sum = pd.read_csv(csv)
    #Summary table
    qsummarytable = sum.groupby('Barcode_arrangement').describe()
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
def scVsLen(csv):
    #Bam
    print("This is score vs summary summary")
    c=pd.read_csv(csv)
   
    fig1 =px.density_heatmap(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig2 =px.scatter(data_frame=c,x=c.loc[:,'Sequence_length_template'], y=c.loc[:,'Mean_qscore_template'])
    fig1.write_html("../density.html")
    fig2.write_html("../scatter.html")
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

    lenSum(csv)
    scSum(csv)
    scVsLen(csv)

def fqToHtml() :
    csv=fqToCsv()  
    csvToHtml(csv)
