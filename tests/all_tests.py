#!/usr/bin/env python
#
# suite of unit tests.
# run suite by executing this file
#

import sys, os
os.chdir('tests')
sys.path.append(os.path.expanduser('~cggroup/repos/PyCogent-trunk'))

import doctest, cogent.util.unit_test as unittest, sys, os

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

def suite():
    modules_to_test = [
        'test_bowtie_output',
        'test_inline_stats',
        'test_light_seq',
        'test_parse_fastq',
        'test_parse_psl',
        'test_region_count',
        'test_syntax',
        'test_cache_lane_counts',
        'test_jackknife',
        'test_subregion_map'
        ]

    alltests = unittest.TestSuite()

    for module in modules_to_test:
        test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

class BoobyTrappedStream(object):
    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError, "Output not allowed in tests"


if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        orig = sys.stdout
        if '--output-ok' in sys.argv:
            sys.argv.remove('--output-ok')
        else:
            sys.stdout = BoobyTrappedStream(orig)
        try:
            unittest.main(defaultTest='suite', argv=sys.argv)
        finally:
            sys.stdout = orig
