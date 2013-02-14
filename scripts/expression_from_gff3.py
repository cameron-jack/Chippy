__author__ = "Cameron Jack"
__copyright__ = "Copyright 2011-2012, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"
__version__ = '638'

import sys

usage_note = '\nUsage: python expression_from_gff3.py <data.gff> [data_out.exp]\n'

if len(sys.argv) == 2:
    if sys.argv[1].lower() == '-h' or sys.argv[1].lower() == '--help':
        print ''
        print 'Expression_from_gff3.py'
        print '-----------------------'
        print 'Expression_from_gff3.py takes a standard gff3 formatted file '
        print 'containing Ensembl standard gene entries and log transformed '
        print 'expression values in the SCORE field, and converts these to '
        print "ChipPy's own expression data standard. This way anyone can "
        print 'write simple scripts to convert their own data to work with '
        print 'ChipPy. Please read the code in this file for information on '
        print 'how to format expression data for ChipPy.'
        print usage_note
        print 'If you do not specify an optional output, the new file will '
        print "have the same name as the original, but with the '.exp' extension"
        print ''
        sys.exit()
if len(sys.argv) != 2 and len(sys.argv) != 3:
    print usage_note
    sys.exit()

# GFF3 example:
# chrV    Hobson_2012     WT_avg_RNAPII_gene_and_100bp_upstream   85675   86149
# 1166.518261     +       0       ID=YEL034W
try:
    infile = open(sys.argv[1], 'r')
except IOError:
    print 'File not found:', sys.argv[1]
    print usage_note
    sys.exit()

outfile = None

# create output file
if len(sys.argv) > 2:
    outfile = open(sys.argv[2], 'w')
else:
    fn = sys.argv[1]
    fn_parts = fn.split('.')
    fn_parts[-1] = 'exp'
    fn = '.'.join(fn_parts)
    outfile = open(fn, 'w')

# write header
outfile.write('ENSEMBL\tprobeset\texp\n')

# grab and write data
for i, line in enumerate(infile):
    line_parts = line.split('\t')
    attributes = line_parts[8].split(';')
    # we need the gene_id
    id_field = [val for val in attributes if 'ID=' or 'id=' in val][0]
    # ENSEMBL field
    ensembl_id = id_field.split('=')[1].strip()
    # exp field
    exp = line_parts[5]
    # Need to create:
    # ENSEML    probeset    exp
    outfile.write(ensembl_id + '\tprobe_'+str(i)+'\t' + exp + '\n')

outfile.close()
print 'Wrote',i+1,'lines to:',sys.argv[2]

### End of expression_from_gff3.py


