"""given an rdump, find the within chromosome gene indexes for the listed
genes"""
from __future__ import division
import os
from cogent import LoadTable

from util import data_dir
from parse_r_dump import RDumpToTable

# we first get the rdumped genes and the ensembl transcript to gene id table
def get_gene_ids(rdump_filename):
    """returns gene IDs ordered by their rank in the rdump file"""
    rdump = RDumpToTable(rdump_filename)
    rdump = rdump.filtered(lambda x: '---' not in x, columns='ENSEMBL')
    transcript_ids = rdump.getRawData('ENSEMBL')
    transcript_to_gene = LoadTable(os.path.join(data_dir, 'ensembl_transcript_gene_ids.txt'),sep='\t')
    transcript_to_gene = transcript_to_gene.getRawData(['TranscriptId', 'GeneId'])
    transcript_to_gene = dict(transcript_to_gene)
    rows = []
    problem = 0
    for index, transcript_group in enumerate(transcript_ids):
        gene_ids = set()
        for transcript_id in transcript_group:
            try:
                gene_id = transcript_to_gene[transcript_id]
            except KeyError:
                continue
            gene_ids.update([gene_id])
        
        if len(gene_ids) != 1:
            problem += 1
            print '\nThe following transcript ids: %s\nMatch multiple genes: %s' % (str(transcript_group), gene_ids)
            continue
        rows += [list(gene_ids)[0]]
    print 'There were %d problem transcripts:' % (problem)
    print 'Total valid IDS: %s' % len(rows)
    return rows

def dump_gene_ids(rdump_filename, outfile_name):
    """dumps Ensembl StableId's to outfile_name"""
    gene_ids = get_gene_ids(rdump_filename)
    outfile = open(outfile_name, 'w')
    outfile.write('\n'.join(gene_ids)+'\n')
    outfile.close()

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

    script_info = {}
    descr = "Given an R-output file generate a file of Ensembl IDs."

    script_info['brief_description']= descr
    script_info['script_description'] = descr
    script_info['version'] = '0.1.alpha'
    script_info['script_usage']=[]
    script_info['script_usage'].append(
        ("Example 1","""Control lane 7; Treatment Lane 8""",
        """python segment_count.py -r r_gene_file -o outfile_name"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-r','--r_gene_file',
                    help='The R generated output file containing gene info.'),
        make_option('-o','--outfile_name',
                    help='Name of output file'),
        ]

    parser, opts, args = parse_command_line_parameters(**script_info)

    dump_gene_ids(opts.r_gene_file, opts.outfile_name)
    
