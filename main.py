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
    sum_command.add_argument("-c","--save2csv",type=str,default=None,
                            help="generate only .csv summary")
  
    #print(parser.print_help())
    return parser
def main():
    parser = argparserLocal()
    args=parser.parse_args()
    
    if args.report:
        
        fqToHtml(args.report)
    
    elif args.save2csv:
        
        fqToCsv(args.save2csv)
    
            
            
            

if __name__ == "__main__":
    main()
    
    


