#!/usr/bin/env python
#
# suite of unit tests.
# run suite by executing this file
#

import sys, os

p = os.path.abspath(__file__)
p = p.split('tests')[0]
sys.path.append(p)
sys.path.append(os.path.expanduser('~cggroup/repos/PyCogent-trunk'))
os.chdir(os.path.join(p, 'tests'))

import doctest, cogent.util.unit_test as unittest, sys, os
from cogent.util.misc import app_path

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
        'test_read_count',
        'test_collection',
        'test_db',
        'test_parse_bed',
        'test_parse_fastq',
        'test_parse_rdump',
        'test_parse_sam',
        'test_syntax',
        'test_util'
        ]
    
    # CAUTION: dodgy hack ahead .. !
    # dodgy because more than one app required
    apps = [('bwa', 'test_fastq_pipe')]
    for app, test_name in apps:
        if app_path(app):
            modules_to_test.append(test_name)
        else:
            print >> sys.stderr, "Can't find %s executable: skipping test" % app
    

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
        raise RuntimeError, "Output not allowed in tests, use --output-ok"


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
