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

