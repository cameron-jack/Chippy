from __future__ import division

import os

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

# bad design, dependency on specific path relationships needs to be completely
# removed
project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
project_dir = project_dir.split('chippy')[0]
data_dir = os.path.join(project_dir, 'data')
src_dir = os.path.join(project_dir, 'src')

DEFAULT_USER_DIR = os.path.expanduser('~/Desktop')

import sys
sys.path.append(src_dir)

def make_even_groups(data, num_per_group, limit=None):
    """returns data split into groups of size num_per_group"""
    grouped = []
    if limit is None:
        limit = len(data)
    
    for start in range(0, limit, num_per_group):
        end = start + num_per_group
        if end > limit:
            break
        group = data[start: end]
        grouped += [group]
    
    return grouped

def get_centred_coords(length, window_size):
    """returns start, end for a window centred on a fragment with length"""
    assert 2 * window_size <= length, 'Invalid window size'
    start = (length//2)- window_size
    end = start + 2 * window_size
    return start, end

def unique_records(data, column=None):
    """returns unique records in original order"""
    try: # it's a Table
        stable_ids = data.getDistinctValues(column)
    except AttributeError:
        stable_ids = set(data)
    
    def is_unique(x, stable_ids=stable_ids):
        if x in stable_ids:
            stable_ids.remove(x)
            return True
        return False
    
    try: # it's a Table
        result = data.filtered(is_unique, columns=column)
    except AttributeError:
        result = filter(is_unique, data)
    
    return result

# following needs to be moved to client code
def get_num_lines(filename):
    # from http://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    f = open(filename, 'U')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read # loop optimization

    buf = read_f(buf_size)
    while buf:
        lines += buf.count('\n')
        buf = read_f(buf_size)
    
    f.close()
    return lines

def create_path(path):
    """creates dir path"""
    try:
        os.makedirs(path)
    except OSError, e:
        pass

def make_cl_command(args):
    """returns string that would be used to invoke script on the command
    line."""
    new = []
    for arg in args:
        if ' ' in arg:
            arg = '"%s"' % arg
        new.append(arg)
    return ' '.join(new)


def dirname_or_default(path):
    """returns dirname of path or, if not specified, a default save
    location"""
    if '/' not in path:
        dirname = DEFAULT_USER_DIR
    else:
        dirname = os.path.dirname(path)
    return dirname

def just_filename(path):
    """whether path is just a filename"""
    return os.path.basename(path) == path

def safe_save_path(path):
    """makes a saf save path if required"""
    if not just_filename(path):
        return path
    
    return os.path.join(DEFAULT_USER_DIR, path)

def grouped_by_chrom(genes):
    """returns dict with genes grouped into chromosome
    
    Assumes gene instances have a coord_name attribute corresponding to chrom"""
    chrom_ordered = {}
    for gene in genes:
        try:
            chrom_ordered[gene.coord_name].append(gene)
        except KeyError:
            chrom_ordered[gene.coord_name] = [gene]
    
    return chrom_ordered

