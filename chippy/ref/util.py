__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

mouse_chroms = tuple(map(str, range(1,20)+['X', 'Y', 'MT']))
human_chroms = tuple(map(str, range(1,23)+['X', 'Y', 'MT']))
chroms = dict(mouse=mouse_chroms, human=human_chroms)
