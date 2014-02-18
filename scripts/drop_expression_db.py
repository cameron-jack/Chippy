#!/usr/bin/env python
import sys
sys.path.extend(['..'])

import chippy.express.db_query as db_query
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.2'

script_info = {}
script_info['title'] = 'Remove an expression study'
script_info['script_description'] = 'Remove an expression study and all '+\
            'associated linked objects.'
script_info['brief_description'] = 'Remove an expression study and all '+\
            'associated linked objects.'
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= 'None.'

# Process command-line arguments
req_args = ['sample']
opt_args = []
pos_args = ['db_path']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args,
    positional_args=pos_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def main():
    rr = RunRecord('drop_expression_db')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse(window_title='Drop Expression Data')
    session = db_query.make_session(args.db_path)

    if db_query.drop_sample_records(session, args.sample):
        rr.addInfo('Removing ' + args.sample, 'Success')
    else:
        rr.addWarning('Removing ' + args.sample, 'Failure')

    rr.display()

if __name__ == "__main__":
    main()

