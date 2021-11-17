from report.read.readFastq import fqToCsv

#16nov

def lenSum(csv):
    print("This is read length summary")
def scSum(csv):
    print("This is score summary")
def scVsLen(csv):
    print("This is score vs summary summary")

def csvToHtml(csv):

    lenSum(csv)
    scSum(csv)
    scVsLen(csv)

def fqToHtml() :
    csv=fqToCsv()  
    csvToHtml(csv)
