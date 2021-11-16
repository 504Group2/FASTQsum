from report.summary.feature import fqToHtml
from report.summary.readFastq import fqToCsv
def argparserLocal():
    from argparse import ArgumentParser
    parser=ArgumentParser(prog='fastqsum',description='Summarize fastq to html report')
    
    subparsers = parser.add_subparsers(
        title='commands',description='Please chose command below:',
        dest='command'    
    )
    subparsers.required = True
    sum_command = subparsers.add_parser('sum',help='Generate html report')
    sum_command.add_argument("-c","--save2csv",type=str,default=None,
                            help="generate only .csv summary")
   
    #print(parser.print_help())
    return parser
def main():
    parser = argparserLocal()
    args=parser.parse_args()
    if args.command == 'sum':
        if args.csv :
            fqToCsv()
        else :
            fqToHtml()

if __name__ == "__main__":
    main()
    
    


