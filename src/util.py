from __future__ import division

import os

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(project_dir, 'data')
src_dir = os.path.join(project_dir, 'src')

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
