from cogent import LoadTable

from segment_count import get_gene_coords
from make_counts import get_read_counts_bowtie, get_file_length

# read in the chromosome lengths and create a dictionary
chrom_lengths = LoadTable('../data/mouse_chrom_lengths_release_58.txt',
                          sep='\t')
chrom_lengths = dict(chrom_lengths.getRawData(['chrom', 'length']))

# input files
treatment_file = '../../../tremethick/s_7.map'
control_file = '../../../tremethick/s_8.map'
gene_coords_file = '../data/mouse_gene_coords.txt'

# The area from the TSS that we are interested in
window_size = 2000

print 'Control Analysis...'
get_file_length(control_file)
for chrom in chrom_lengths:

    chrom_str = 'chr' + str(chrom)
    print 'Getting Counts for ' + chrom_str
    counter = get_read_counts_bowtie(control_file, chrom_str,
                                     chrom_lengths[chrom])
    mouse_gene_coords = get_gene_coords(gene_coords_file, chrom)
    print 'Saving binary for ' + chrom_str
    counter.save('s_8_counts_%s'%chrom_str, mouse_gene_coords, window_size)

print 'Treatment Analysis...'
get_file_length(treatment_file)
for chrom in chrom_lengths:
    chrom_str = 'chr' + str(chrom)
    print 'Getting Counts for ' + chrom_str
    counter = get_read_counts_bowtie(treatment_file, chrom_str,
                                     chrom_lengths[chrom])
    mouse_gene_coords = get_gene_coords(gene_coords_file, chrom)
    print 'Saving binary for ' + chrom_str
    counter.save('s_7_counts_%s'%chrom_str, mouse_gene_coords, window_size)

