import util

from cogent import LoadTable

from region_count import CacheLaneCounts
from gene_data import GeneIndexes

def get_counts(stable_id_file, ctl_path, ctl_lane, trt_path, trt_lane):
    """returns a Table instance with column headers the stable IDs and rows
    the tag counts"""
    # get the genes, sorted by chromosome
    stable_ids = [i.strip() for i in open(stable_id_file).readlines()]
    genes = sorted(map(GeneIndexes, stable_ids))
    trt_cache = CacheLaneCounts(trt_lane, trt_path)
    ctl_cache = CacheLaneCounts(ctl_lane, ctl_path)
    rows = []
    chrom = None
    header = None
    for gene in genes:
        if gene.chrom != chrom:
            del(trt_cache[chrom])
            del(ctl_cache[chrom])
            chrom = gene.chrom
        
        # get the saved counts for this gene and convert to list
        gene_trt_counts = trt_cache[chrom][gene.index].tolist()
        gene_ctl_counts = ctl_cache[chrom][gene.index].tolist()
        
        if header is None:
            # start with the range of values
            length = len(gene_trt_counts)
            header = map(str, range(-(length/2), (length/2)))
            header = ['stable_id', 'c/t'] + header
        
        # insert critical identifying info into each list
        gene_trt_counts = [gene.stable_id, 't'] + gene_trt_counts
        gene_ctl_counts = [gene.stable_id, 'c'] + gene_ctl_counts
        rows.append(gene_trt_counts)
        rows.append(gene_ctl_counts)
    table = LoadTable(header=header, rows=rows)
    return table

if __name__ == "__main__":
    import os
    path_join = os.path.join
    
    parent_dir = '/Users/gavin/DevRepos/RegionAnalysis/'
    stable_id_file = path_join(parent_dir, 'results/g2-vs-g1/ensembl_ids.G2.lt.G1.txt')
    ctl_dir=path_join(parent_dir, "results/run2/counts/control")
    trt_dir=path_join(parent_dir, "results/run2/counts/treatment")
    ctl_lane=8
    trt_lane=7
    
    result = get_counts(stable_id_file, ctl_dir, ctl_lane,
                    trt_dir, trt_lane)
    print result[:10, :10]
