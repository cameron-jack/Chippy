from cogent.parse.fasta import MinimalFastaParser, Info
from parse_psl import MinimalPslParser
from parse_fastq import MinimalFastqParser

def get_corrupt_seq_names(psl_name, test_run):
    psl_parser = MinimalPslParser(psl_name)
    psl_parser.next()
    psl_parser.next()
    
    num = 0
    to_trim = set()
    for record in psl_parser:
        name = record[9]
        num += 1
        to_trim.update([name])
        if test_run and num >= 1000:
            break
    print to_trim
psl_name = open("/Users/hussain/Python_prog/Trehmethick_Data/result_6_100.psl")
get_corrupt_seq_names(psl_name, test_run)
 
