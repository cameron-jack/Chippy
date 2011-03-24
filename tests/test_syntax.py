"""simple test to catch syntax errors from incorrect editing of conflicts"""
import sys
sys.path.extend(['../src', '..'])

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod


_modules_to_import = [
    'ensembl_query',
    'ensembl_transcript_gene_idmap',
    'fastq_to_fasta',
    'get_pristine_seqs',
    'inline_stats',
    'chippy.parse.light_seq',
    'plot_quality',
    'make_counts',
    'chippy.parse.fastq',
    'chippy.parse.sam',
    'region_analysis',
    'region_count',
    'segment_count',
    ]

for mod in _modules_to_import:
    my_import(mod)
