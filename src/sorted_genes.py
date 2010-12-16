from cogent import LoadTable

genes = LoadTable('../data/mouse_gene_coords.txt', sep='\t')
genes = genes.sorted(columns=['CoordName', 'Start'])
columns = [h for h in genes.Header if h != 'Index']
genes = genes.getColumns(columns)
print genes[:10]


def Counter():
    global count
    count = -1
    def call(val):
        global count
        count += 1
        return count
    
    return call

tables = []
for chrom_name in range(1,20)+['X', 'Y']:
    chrom = genes.filtered(lambda x: x == chrom_name,
                                    columns='CoordName')
    chrom = chrom.withNewColumn('Index', Counter(), columns)
    tables.append(chrom)

table = tables[0].appended(None, tables[1:])
table.writeToFile('../data/mouse_gene_coords.txt', sep='\t')
