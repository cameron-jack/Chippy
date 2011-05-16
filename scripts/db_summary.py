#!/usr/bin/env python
"""reports a summary of the database"""

from __future__ import division

import os, sys, glob
sys.path.extend(['../../src'])

from optparse import make_option

from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import parse_command_line_parameters

from chippy.express import db_query, db_schema
from chippy.core.read_count import get_combined_counts

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

ensembl_release='58'
session = db_query.make_session('sqlite:///%s' % db_path)
samples = db_query.get_samples(session)
if not samples:
    samples = [None]

script_info = {}

script_info['title'] = 'Report what files have been related to a sample'
script_info['script_description'] = "Prints a table showing what files have"\
                                    " been related to a sample."
script_info['version'] = __version__

script_info['required_options'] = [
make_option('-c', '--sample', type='choice',
           help='Choose the expression study [default: %default]',
           choices=[str(s) for s in samples]),
]

script_info['authors'] = __author__

def _samples_name(sample):
    for s in samples:
        if str(s) == sample:
            return s

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)
    if opts.sample is None:
        raise RuntimError('No samples available')
    
    sample = _samples_name(opts.sample)
    reffiles = session.query(db_schema.ReferenceFile).join(db_schema.Sample).\
            filter(db_schema.Sample.sample_id==sample.sample_id).all()
    
    if not reffiles:
        print 'No reference files have been linked to this sample'
    else:
        rows = []
        for reffile in reffiles:
            rows.append([reffile.date, reffile.name])
        table = LoadTable(header=['Date added', 'Reference file path'],
                rows=rows, title='Files related to "%s"' % sample)
        print table

if __name__ == "__main__":
    main()

