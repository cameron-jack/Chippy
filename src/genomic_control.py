import re
from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Species, Genome

from light_seq import LightSeq

def get_gene_coords(infile_name, window_size, seq_length, chrom_name=None):
    """Given a file with stableIDs, start and end coordinates for a set of genes
    and a window size of interest, and the length of each read, get appropriate
    coordinates"""

    table = LoadTable(infile_name, sep='\t')
    if chrom_name is not None:
        table = table.filtered(lambda x: str(x) == str(chrom_name),
                               columns='CoordName')
    table = table.sorted(columns=['CoordName', 'Start'])
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
        else:
            TSS = row['Start']

        start_coord = TSS - (window_size + (seq_length-1))
        end_coord = TSS + (window_size + (seq_length-1))
        rows += [[row['StableId'], row['CoordName'], start_coord, end_coord]]

    return rows

def find_all_instances(findThis, inThis):
    """Provide indices of where findThis occurs in inThis"""

    indices = []
    for m in re.finditer(findThis, inThis):
        indices += [m.start()]

    if len(indices) is 0:
        indices = None

    return indices

def generate_control_seqs(infile_name, window_size=10000, seq_length=75,
                          outfile_name=None, specie='Mouse',chrom_name=None,
                          account=None):
    """After determining the appropriate coordinates around transcription start
    sites, download the sequences from Ensembl, and generate a set of seqs
    based on a sliding window of 1 for the entire window period."""

    coordinates = get_gene_coords(infile_name, window_size, seq_length,
                                  chrom_name)

    if account is None:
        account = HostAccount('cg.anu.edu.au','compgen','compgen')

    specie = Genome(Species=specie, Release = 58, account=account)

    if outfile_name is None:
        outfile_name = 'control_seqs'
    if chrom_name is None:
        outfile_name = outfile_name+'.fastq'
    else:
        outfile_name = outfile_name + '_chr%s.fastq'%(str(chrom_name))

    outfile = open(outfile_name, 'w')
    seq_qual_default = 'h'*seq_length

    for [stable_id, chrom, start, end] in coordinates:
        print 'Processing %s on Chromosome %s'%(stable_id, str(chrom))
        gene = specie.getRegion(CoordName=str(chrom), Start=start, End=end)
        isDegenerate = gene.Seq.isDegenerate()

        num = 0
        for seq in gene.Seq.slidingWindows(seq_length,1,0):
            num += 1
            seq_str = str(seq)
            seq_qual = seq_qual_default

            if isDegenerate:
                indices = find_all_instances('N', seq_str)
                if indices is not None:
                    qualList = list(seq_qual)
                    for index in indices:
                        qualList[index] = 'C'
                    seq_qual = ''.join(qualList)

            seq_object = LightSeq(seq_str, stable_id+'_%d'%num, seq_qual)
            outfile.write(seq_object.toFastq() + '\n')
    outfile.close()

#generate_control_seqs_2('../data/mouse_gene_coords.txt', 10000, 75)

if __name__ == "__main__":
    from cogent.util.misc import parse_command_line_parameters
    from optparse import make_option

    script_info = {}
    descr = "Given a window size, download portions of a gene around the TSS "\
            "for each gene on every chromosome. Then for each portion, create "\
            "a set of sequences of a particular length, using a sliding "\
            "window approach (sliding by 1 bp)"

    script_info['brief_description']= descr
    script_info['script_description'] = descr
    script_info['version'] = '0.1.alpha'
    script_info['script_usage']=[]
    script_info['script_usage'].append(
        ("Example 1","""Mouse with default options:""",
        """python genomic_control.py -i mouse_coords.txt"""))
    script_info['script_usage'].append(
        ("Example 2","""Mouse, only chromosome 1:""",
        """python genomic_control.py -i mouse_coords.txt -c 1"""))

    script_info['help_on_no_arguments'] = True
    script_info['required_options'] = [
        make_option('-i','--input_file',
                    help='A gene coordinates file')]

    script_info['optional_options'] = [
        make_option('-w','--window_size', type='int',
                    dest='window_size', default = 10000,
                    help='Number of bases on either site of TSS that you are '\
                    'interested in ' +'[default: %default]'),
        make_option('-o','--output_fileroot', default = 'control_seqs',
                    help='Path and fileroot for the output fastq file(s) '\
                    'Note that if a chromosome is specified, _chrN will be '\
                    'added to the fileroot. ' + '[default: %default]'),
        make_option('-c', '--chrom_name', help="the chromosome name"),
        make_option('-l', '--length', type='int', dest='seq_length',
                    default=75, help="Length of each generated "\
                    "sequence. " + '[default: %default]'),
        make_option('-s', '--specie', default='Mouse',
                    help="Ensembl recognisable specie (common name)" +
                    '[default: %default]')]

    parser, opts, args = parse_command_line_parameters(**script_info)

    generate_control_seqs(opts.input_file, opts.window_size, opts.seq_length,
                          opts.output_fileroot, opts.specie, opts.chrom_name)


