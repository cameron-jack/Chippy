import sys
sys.path.append('../src')

import numpy

from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from export_counts import get_counts

# we should be able to load counts from saved numpy files
# export a subset save to file, reload and the results match the original
# numpy data
class ExportCountsTests(TestCase):
    def test_one(self):
        pass

if __name__ == "__main__":
    main()