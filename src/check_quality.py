from cogent.parse.psl import MinimalPslParser
from cogent.parse.fasta import MinimalFastaParser

def eval_psl(psl_name):
    to_trim = set()
    psl_parser = MinimalPslParser(psl_name)
    psl_parser.next()
    psl_parser.next()
    num = 0
    for record in psl_parser:
        name = record[9]
        num += 1
        to_trim.update([name])
    
    print 'Number seqs with hits to primer', num
    print 'Number unique seqs with hits to primer', len(to_trim)

def eval_fasta(fasta_name):
    num = 0
    for line in open(fasta_name):
        if line[0] == '>':
            num += 1
    print 'Total number seqs', num

psl_name = '../data/result_s_6.psl'
# eval_psl(psl_name)


fasta_name = '../data/s_6.fasta'
eval_fasta(fasta_name)
