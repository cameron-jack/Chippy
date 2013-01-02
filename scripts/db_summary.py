#!/usr/bin/env python
"""reports a summary of the database - needs more work"""

from __future__ import division

import os, sys, glob
sys.path.extend(['../../src'])

from cogent import LoadTable
from cogent.util.progress_display import display_wrap
from cogent.util.misc import parse_command_line_parameters

from chippy.express import db_query, db_schema
from chippy.core.read_count import get_combined_counts
from chippy.util.command_args import Args

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012, Gavin Huttley, Cameron Jack, Anuj Pahwa"
__credits__ = ["Gavin Huttley, Cameron Jack"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Development"

def set_environment():
    pos_args = ['db_path']
    req_args = ['sample']
    args = Args(positional_args=pos_args, required_args=req_args)

    script_info = {}

    script_info['title'] = 'Report what files have been related to a sample'
    script_info['script_description'] = "Prints a table showing what files have"\
                                    " been related to a sample."
    script_info['version'] = __version__
    script_info['authors'] = __author__
    return args, script_info

def main():
    args, scripts_info = set_environment()
    session = db_query.make_session('sqlite:///' + args.db_path)

    sample = args.sample
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

