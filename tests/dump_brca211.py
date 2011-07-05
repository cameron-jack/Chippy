"""dumps exon 11 sequence as a fast file"""
import sys
sys.path.append('../src')
from light_seq import LightSeq

from cogent.db.ensembl import HostAccount, Genome
account = HostAccount('cg.anu.edu.au', 'gavin', 'gavin')
mouse = Genome('mouse', Release=58, account=account)
brca2 = mouse.getGeneByStableId('ENSMUSG00000041147')
for exon in brca2.CanonicalTranscript.Exons:
    if exon.StableId == 'ENSMUSE00000258226':
        break

print exon
seq = exon.Seq
seq_qual_default = 'h'*seq_length

seq_length = 75
name_template = '%s_%s' % (brca2.StableId, exon.Location.Start) + '_%s_%s'
fastq = []
num_plus = 0
num_minus = 0
for slice_start in range(0, len(exon), seq_length):
    fwd = seq[slice_start: slice_start+seq_length]
    rev = seq[slice_start: slice_start+seq_length].rc()
    if len(fwd) < 50:
        continue
    
    name = name_template % ('plus', slice_start)
    seq_object = LightSeq(str(fwd), name, seq_qual_default[:len(fwd)])
    fastq.append(seq_object.toFastq())
    num_plus += 1
    
    name = name_template % ('minus', slice_start)
    seq_object = LightSeq(str(rev), name, seq_qual_default[:len(fwd)])
    fastq.append(seq_object.toFastq())
    num_minus += 1


outfile = open('data/brca2-11.fastq', 'w')
outfile.write('\n'.join(fastq) + '\n')
outfile.close()

outfile = open('data/brca2-11.fasta', 'w')
outfile.write(seq.toFasta())
outfile.close()

print num_plus, num_minus

print 'Done'
