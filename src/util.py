from __future__ import division

import os

project_dir = os.path.dirname(os.path.dirname(__file__))
data_dir = os.path.join(project_dir, 'data')
src_dir = os.path.join(project_dir, 'src')

import sys
sys.path.append(src_dir)

def make_even_groups(data, num_per_group):
    """returns data split into groups of size num_per_group"""
    grouped = []
    N = len(data)
    for start in range(0, N, num_per_group):
        end = start + num_per_group
        if end > N:
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

