import re
from cogent import LoadTable
from cogent.db.ensembl import HostAccount, Species, Genome

from light_seq import LightSeq

def get_gene_coords(infile_name, window_size, seq_length):
    """Given a file with stableIDs, start and end coordinates for a set of genes
    and a window size of interest, and the length of each read, get appropriate
    coordinates"""

    table = LoadTable(infile_name, sep='\t')
    table = table.sorted(columns=['CoordName', 'Start'])
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
        else:
            TSS = row['Start']

        start_coord = TSS - (window_size + seq_length)
        end_coord = TSS + (window_size + seq_length)
        rows += [[row['StableId'], row['CoordName'], start_coord, end_coord]]

    return rows

def generate_control_seqs(infile_name, window_size, seq_length, outfile_name=None,
                          account=None, specie='Mouse'):
    """After determining the appropriate coordinates around transcription start
    sites, download the sequences from Ensembl, and save them to outfile as a
    fastq formatted file"""

    coordinates = get_gene_coords(infile_name, window_size, seq_length)

    if account is None:
        account = HostAccount('cg.anu.edu.au','compgen','compgen')

    mouse = Genome(Species=specie, Release = 58, account=account)
    if outfile_name is None:
        outfile_name = 'mouse_control_seqs.fastq'
    outfile = open(outfile_name, 'w')

    for [stable_id, chrom, start, end] in coordinates:
        print 'Processing %s on Chromosome %s'%(stable_id, str(chrom))
        gene = mouse.getRegion(CoordName=str(chrom), Start=start, End=end)
        isDegenerate = gene.Seq.isDegenerate()

        num = 0
        for seq in gene.Seq.slidingWindows(seq_length,1,0):
            num += 1
            seq_str = str(seq)
            seq_qual = 'h'*seq_length

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


def find_all_instances(findThis, inThis):
    """Provide indices of where findThis occurs in inThis"""

    indices = []
    for m in re.finditer(findThis, inThis):
        indices += [m.start()]

    if len(indices) is 0:
        indices = None

    return indices


#test = 'AAACCACAAGAGTGACTTGG'
#print find_all_instances('N', test)
generate_control_seqs('../data/mouse_gene_coords.txt', 10000, 75)

