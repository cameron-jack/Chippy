import numpy
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator
from cogent import LoadTable
from cogent.util.progress_display import display_wrap

from region_count import RegionCounts
from get_gene_data import get_gene_indices

def get_gene_coords(infile_name, chrom_name):
    """returns TSS, strand for the nominated chromosome"""
    table = LoadTable(infile_name, sep='\t')
    table = table.filtered(lambda x: x == chrom_name, columns='CoordName')
    table = table.sorted(columns=['CoordName', 'Start'])
    rows = []
    for row in table:
        strand = row['Strand']
        if strand == -1:
            TSS = row['End']
        else:
            TSS = row['Start']
        
        rows += [(TSS, strand)]
    
    return rows

def get_binned_counts(a, bin_size=150):
    """returns counts in a bin"""
    result = [a[i: i+bin_size].sum() for i in range(a.shape[0]-bin_size)]
    return numpy.array(result)

@display_wrap
def run(rdumpfilename, title, ui):
    win_size = 2000
    genes = get_gene_indices(rdumpfilename)
    num_genes = len(genes)
    treatment_template = '../../RegionAnalysis/s_7-window_5000-chr%s.npy'
    control_template = '../../RegionAnalysis/s_8-window_5000-chr%s.npy'
    treatments = {}
    controls = {}
    trt_binned = None
    ctl_binned = None
    for gene in genes:
        if gene.chrom in treatments:
            treatment = treatments[gene.chrom]
            control = controls[gene.chrom]
        else:
            treatment = numpy.load(treatment_template % gene.chrom)
            control = numpy.load(control_template % gene.chrom)
            treatments[gene.chrom] = treatment
            controls[gene.chrom] = control
            if trt_binned is None:
                trt_binned = numpy.zeros(control.shape[1], float)
                ctl_binned = numpy.zeros(control.shape[1], float)
        ctl_binned += control[gene.index]
        trt_binned += treatment[gene.index]
    
    trt_binned /= num_genes
    ctl_binned /= num_genes
    
    diff = trt_binned - ctl_binned
    
    sd = numpy.sqrt(trt_binned + ctl_binned)
    norm_diff = diff / sd
    norm_diff /= num_genes
    norm_diff = get_binned_counts(norm_diff)
    x = numpy.arange(-(control.shape[1]/2)+75, (control.shape[1]/2)-75)
    fig = pyplot.figure(figsize=(10,5))
    minor_locator = MultipleLocator(100)
    pyplot.plot(x, norm_diff)
    pyplot.xlabel('Position relative to TSS')
    pyplot.ylabel('Normalised difference (averaged by gene)')
    pyplot.title("%s: N=%s" % (title, len(genes)))
    pyplot.gca().xaxis.set_minor_locator(minor_locator)
    pyplot.setp(pyplot.gca(), ylim=(-0.5,2.5))
    pyplot.grid(True)
    pyplot.savefig('%s.pdf' % rdumpfilename)
    print '\n\nDone!'

if __name__ == "__main__":
    rdump_filename1 = '../g2-vs-g1/paa.G2.gt.G1.sig.sort.2FC.txt'
    title1 = 'G2 > G1'
    rdump_filename2 = '../g2-vs-g1/paa.G2.lt.G1.sig.sort.2FC.txt'
    title2 = 'G2 < G1'
    
    run(rdump_filename1, title1)
