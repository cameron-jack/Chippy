from cogent import LoadTable

from segment_count import get_gene_coords
from make_counts import get_read_counts_bowtie

# read in the chromosome lengths and create a dictionary
chromLengths = LoadTable('../data/mouse_chrom_lengths_release_58.txt', sep='\t')
chromLengthsDict = dict(chromLengths.getRawData(['chrom', 'length']))
chroms = chromLengthsDict.keys()

# input files
treatmentFile = '../../../tremethick/s_7.map'
controlFile = '../../../tremethick/s_8.map'
geneCoordsFile = '../data/mouse_gene_coords.txt'

# The area from the TSS that we are interested in
window_size = 2000

print 'Control Analysis...'
for chrom in chroms:

    chrom_str = 'chr' + str(chrom)
    print 'Getting Counts for ' + chrom_str
    counter = get_read_counts_bowtie(controlFile, chrom_str, chromLengthsDict[chrom])
    mouseGeneCoords = get_gene_coords(geneCoordsFile, chrom)
    print 'Saving binary for ' + chrom_str
    counter.save('s_8_counts_%s'%chrom_str, mouseGeneCoords, window_size)

print 'Treatment Analysis...'
for chrom in chroms:
    chrom_str = 'chr' + str(chrom)
    print 'Getting Counts for ' + chrom_str
    counter = get_read_counts_bowtie(treatmentFile, chrom_str, chromLengthsDict[chrom])
    mouseGeneCoords = get_gene_coords(geneCoordsFile, chrom)
    print 'Saving binary for ' + chrom_str
    counter.save('s_7_counts_%s'%chrom_str, mouseGeneCoords, window_size)




