from __future__ import division
import os
import util

from cogent import LoadTable

from region_count import CacheLaneCounts
from subregion_map import MapScores
from gene_data import GeneIndexes

def region_counts(lane, data_path, coordinates):
    """returns counts for each region in coordinates"""
    data_path = os.path.abspath(data_path)
    cache = CacheLaneCounts(lane, data_path)
    
    # sort coords, improves data loading performance
    coordinates = sorted(coordinates)
    chrom = coordinates[0].chrom
    rows = []
    for coord in coordinates:
        if coord.chrom != chrom:
            del(cache[chrom])
            chrom = coord.chrom
        counts = cache[chrom][coord.index]
        rows.append(counts)
    
    return rows

def simple_get_tss_counts(stable_ids, lane, data_path, type_, window_size=None):
    stable_ids = util.unique_records(stable_ids, 'StableId')
    data_path = os.path.abspath(data_path)
    genes = []
    for index, stable_id in enumerate(stable_ids):
        try:
            gene = GeneIndexes(stable_id)
            gene.order = index
            genes += [gene]
        except KeyError:
            pass
    
    # we sort the genes
    genes = sorted(genes)
    cache = CacheLaneCounts(lane, data_path, window_size=window_size)
    rows = []
    chrom = genes[0].chrom
    header = None
    for gene in genes:
        if gene.chrom != chrom:
            del(cache[chrom])
            chrom = gene.chrom
        
        # get the saved counts and mappability for this gene and convert to list
        gene_counts = cache[chrom][gene.index].tolist()
        
        if header is None:
            # start with the range of values
            length = len(gene_counts)
            if window_size is not None:
                start, end = util.get_centred_coords(length, window_size)
            else:
                start = 0
                end = length
            
            header = map(str, range(-(length//2), (length//2)))
            header = ['Order', 'StableId', 'Type'] + header[start: end]
            
        
        # insert critical identifying info into each list
        gene_counts = [gene.order, gene.stable_id, type_] + gene_counts[start: end]
        rows.append(gene_counts)
    table = LoadTable(header=header, rows=rows)
    table = table.sorted('Order').getColumns(table.Header[1:])
    return table

def get_counts(stable_id_file, ctl_path, ctl_lane, trt_path, trt_lane,
               score_path, window_size=None):
    """returns a Table instance with column headers the stable IDs and rows
    the tag counts
    
    If window_size provided, will return data +/- that size from middle of counts"""

    # get the genes, sorted by chromosome
    try:
        stable_ids = LoadTable(stable_id_file, sep='\t').getRawData('StableId')
    except:
        stable_ids = [i.strip() for i in open(stable_id_file).readlines()]
    
    stable_ids = util.unique_records(stable_ids, 'StableId')
    
    genes = []
    for index, stable_id in enumerate(stable_ids):
        try:
            gene = GeneIndexes(stable_id)
            gene.order = index
            genes += [gene]
        except KeyError:
            pass
    
    # we sort the genes
    genes = sorted(genes)
    trt_cache = CacheLaneCounts(trt_lane, trt_path, window_size=window_size)
    ctl_cache = CacheLaneCounts(ctl_lane, ctl_path, window_size=window_size)
    mapscores = MapScores(score_path, window_size=window_size)
    rows = []
    chrom = genes[0].chrom
    header = None
    for gene in genes:
        if gene.chrom != chrom:
            del(trt_cache[chrom])
            del(ctl_cache[chrom])
            del(mapscores[chrom])
            chrom = gene.chrom
        
        # get the saved counts and mappability for this gene and convert to list
        gene_trt_counts = trt_cache[chrom][gene.index].tolist()
        gene_ctl_counts = ctl_cache[chrom][gene.index].tolist()
        gene_mapscores = mapscores[chrom][gene.index].tolist()
        
        if header is None:
            # start with the range of values
            length = len(gene_trt_counts)
            if window_size is not None:
                start, end = util.get_centred_coords(length, window_size)
            else:
                start = 0
                end = length
            
            header = map(str, range(-(length//2), (length//2)))
            header = ['Order', 'StableId', 'Type'] + header[start: end]
            
        
        # insert critical identifying info into each list
        gene_trt_counts = [gene.order, gene.stable_id, 't'] + gene_trt_counts[start: end]
        gene_ctl_counts = [gene.order, gene.stable_id, 'c'] + gene_ctl_counts[start: end]
        gene_mapscores =  [gene.order, gene.stable_id,  'm'] + gene_mapscores[start: end]
        rows.append(gene_trt_counts)
        rows.append(gene_ctl_counts)
        rows.append(gene_mapscores)
    table = LoadTable(header=header, rows=rows)
    table = table.sorted('Order').getColumns(table.Header[1:])
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
