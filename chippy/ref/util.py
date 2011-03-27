__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

# we ignore mitochondria
mouse_chroms = map(str, range(1,20)+['X', 'Y'])
human_chroms = map(str, range(1,23)+['X', 'Y'])
chroms = dict(mouse=mouse_chroms, human=human_chroms)
