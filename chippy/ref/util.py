import sys
sys.path.extend(['..'])
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "pre-release"
__version__ = '0.1'

### This file is deprecated. A Chroms table is now kept in the ChippyDB

mouse_chroms = tuple(map(str, range(1,20)+['X', 'Y', 'MT']))
human_chroms = tuple(map(str, range(1,23)+['X', 'Y', 'MT']))
cerevisiae_chroms = tuple(['I','II','III','IV','V','VI','VII','VIII',
        'IX','X','XI','XII','XIII','XIV','XV','XVI','MT'])
chroms = dict(mouse=mouse_chroms, human=human_chroms,
        yeast=cerevisiae_chroms)

# implement as a singleton, so we can use this anywhere
# this implementation may need to be changed if converted to Python3
class Chromosomes:
    __single = None
    def __init__(self, species, chroms):
        # chroms should be a list
        if Chromosomes.__single:
            raise Chromosomes.__single
        Chromosomes.__single = self
        self.species = species
        self.chroms = chroms

def chromHandle(species='mouse', chroms=mouse_chroms):
    """chrom_handle provides the safe way of creating an instance of
    Chrom and returning THE existing instance of Chromosomes. """
    rr = RunRecord('chromHandle')
    try:
        chromsInstance = Chromosomes(species, chroms)
        rr.addInfo('Chromosomes instance created', True)
        rr.addInfo('Chromosomes instance species', species)
        rr.addInfo('Chromosomes instance chromosomes', chroms)
    except Chromosomes, c:
        chromsInstance = c
        if species or chroms:
            rr.addWarning('provided species', species)
            rr.addWarning('provided chromosomes', chroms)
    return chromsInstance
