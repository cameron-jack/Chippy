from cogent import LoadTable

genes = LoadTable('../data/mouse_gene_coords.txt', sep='\t')
genes = genes.sorted(columns=['CoordName', 'Start'])
print genes[:10]
genes.writeToFile('../data/mouse_gene_coords.txt', sep='\t')
