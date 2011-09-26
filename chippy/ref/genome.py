import os
from optparse import make_option
from cogent.util.misc import parse_command_line_parameters
from cogent.db.ensembl import Genome, HostAccount

from chippy.ref.util import chroms

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'


def get_chrom_seqs(species, release, account=None, debug=False):
    """yields sequence objects for the indicated chromosomes from Ensembl"""
    genome = Genome(species, Release=release, account=account)
    for chrom in chroms[species]:
        region = genome.getRegion(CoordName=chrom)
        seq = region.Seq
        name = 'chr_%s' % chrom
        seq.Name = name
        if debug:
            print name
            print repr(seq)
        
        yield seq

script_info = {}
script_info['script_description'] = "Dump species chromosome sequences in fasta format."
script_info['output_description']= "fasta formatted sequences"
script_info['version'] = __version__
script_info['authors'] = 'Gavin Huttley'

script_info['required_options'] = [
 make_option('-s',
            '--species',
            type='choice', 
            help='Species to query [default: %default]',
            default='mouse',
            choices=['mouse', 'human']),
 make_option('-r','--release',
             help='The Ensembl release number [default: %default]',
             default=None),
]

script_info['optional_options'] = [
 make_option('-o','--outdir',
             help='Directory to write files to [default: %default]',
             default='.'),
 make_option('-c','--by_chrom', action='store_true', default=False,
             help='Write separate files for each chromosome [default: %default]'),
 make_option('-H','--hostname',
             help='MySQL server hostname. Will look for ENSEMBL_ACCOUNT by default.',
             default=None),
 make_option('-U',
             '--user',
             help='MySQL username. Will look for ENSEMBL_ACCOUNT by default.',
             default=None),
 make_option('-P',
             '--passwd',
             help='MySQL password. Will look for ENSEMBL_ACCOUNT by default.',
             default=None),
 make_option('-t',
             '--test_run',
             action='store_true',
             help='Do not write any output [default: %default]',
             default=False),
]

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    if None in (opts.hostname, opts.user, opts.passwd):
        assert len(set((opts.hostname, opts.user, opts.passwd))) == 1,\
            'You must provide all MySQL options, or none at all.'
    
    if opts.hostname is not None:
        account = HostAccount(opts.hostname,opts.user,opts.passwd)
    elif 'ENSEMBL_ACCOUNT' in os.environ:
        h, u, p = os.environ['ENSEMBL_ACCOUNT'].split()
        account = HostAccount(h,u,p)
    else:
        account = None
    
    if opts.test_run:
        print account
    
    outdir = os.path.abspath(opts.outdir)
    if not os.path.exists(outdir):
        print 'FAIL: %s directory does not exist' % outdir
        exit(-1)
    
    if not opts.by_chrom:
        outfile_name = os.path.join(outdir, '%s-%s.fasta' % (opts.species, opts.release))
        if not opts.test_run:
            outfile = open(outfile_name, 'w')
    
    if opts.test_run:
        print 'Will write to: %s' % outdir
        if not opts.by_chrom:
            print outfile_name
    
    for chrom in get_chrom_seqs(opts.species, opts.release, account,
                                debug=opts.test_run):
        fasta = chrom.toFasta()
        
        if opts.by_chrom:
            outfile_name = os.path.join(outdir, '%s.fasta' % chrom.Name)
        
        if opts.test_run:
            print 'Will write to: %s' % outfile_name
            break
        
        if opts.by_chrom:
            outfile = open(outfile_name, 'w')
        
        outfile.write(fasta+'\n')
        
        if opts.by_chrom:
            outfile.close()


if __name__ == "__main__":
    main()
