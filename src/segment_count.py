import numpy
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator
from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from region_count import RegionCounts

def get_gene_coords(infile_name, chrom_name, window, just_strand=1):
    """returns 5', TSS, 3', strand"""
    table = LoadTable(infile_name, sep='\t')
    table = table.filtered(lambda x: x == chrom_name, columns='CoordName')
    table = table.sorted()
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
            five = TSS + window
            three = TSS - window
        else:
            TSS = row['Start']
            five = TSS - window
            three = TSS + window
        
        if strand == just_strand:
            rows += [(chrom_name, TSS, five, three, strand)]
    return rows

def get_binned_counts(a, bin_size=150):
    """returns counts in a bin"""
    result = [a[i: i+bin_size].sum() for i in range(a.shape[0]-bin_size)]
    return numpy.array(result)

@display_wrap
def run(ui):
    """docstring for run"""
    win_size = 2000
    gene_TSS = get_gene_coords('../data/mouse_gene_coords.txt', 1, win_size)
    treatment = RegionCounts(filename='../data/s_7_pristine-chr1.npy')
    control = RegionCounts(filename='../data/s_8_pristine-chr1.npy')
    trt_binned = None
    ctl_binned = None
    for i in ui.series(range(len(gene_TSS))):
        chrom, tss, five, three, strand = gene_TSS[i]
        trt_counts = treatment.getCounts(five, three, None)
        ctl_counts = control.getCounts(five, three, None)
        if ctl_binned is None:
            trt_binned = get_binned_counts(trt_counts)
            ctl_binned = get_binned_counts(ctl_counts)
        else:
            trt_binned += get_binned_counts(trt_counts)
            ctl_binned += get_binned_counts(ctl_counts)
    
    # averaged by number of genes
    trt_binned /= len(gene_TSS)
    ctl_binned /= len(gene_TSS)
    
    diff = trt_binned - ctl_binned
    
    sd = numpy.sqrt(trt_binned + ctl_binned)
    norm_diff = diff / sd
    norm_diff /= len(gene_TSS)
    x = numpy.arange(-win_size+75, win_size-75)
    fig = pyplot.figure(figsize=(10,5))
    minor_locator = MultipleLocator(100)
    pyplot.plot(x, norm_diff)
    pyplot.xlabel('Position relative to TSS')
    pyplot.ylabel('Normalised difference (averaged by gene)')
    pyplot.title('Mouse Chrom 1 Genes + strand')
    pyplot.gca().xaxis.set_minor_locator(minor_locator)
    pyplot.grid(True)
    pyplot.savefig('mouse_chrom1_plus_strand.pdf')
    print '\n\nDone!'

if __name__ == "__main__":
    run()
