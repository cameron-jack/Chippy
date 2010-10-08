import os

project_dir = os.path.dirname(os.path.dirname(__file__))
data_dir = os.path.join(project_dir, 'data')
src_dir = os.path.join(project_dir, 'src')

import sys
sys.path.append(src_dir)
