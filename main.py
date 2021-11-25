from report.summary.feature import fqToHtml
from report.read.readFastq import fqToCsv
def argparserLocal():
    from argparse import ArgumentParser
    parser=ArgumentParser(prog='fastqsum',description='Summarize fastq to html report')
    
    subparsers = parser.add_subparsers(
        title='commands',description='Please choose command below:',
        dest='command'    
    )
    subparsers.required = True
    sum_command = subparsers.add_parser('sum',help='Generate html report')
    sum_command.add_argument("-r","--report",type=str,default=None,
                            help="generate complete html report")
    sum_command.add_argument("-c","--save2csv",action='store_true',default=None,
                            help="generate only .csv summary")
  
    #print(parser.print_help())
    return parser
def main():
    parser = argparserLocal()
    args=parser.parse_args()
    filePath=args.report
    if args.command == 'sum':
        if args.save2csv :
            fqToCsv(filePath)
        else :
            fqToHtml(filePath)
            
            

if __name__ == "__main__":
    main()
    
    


