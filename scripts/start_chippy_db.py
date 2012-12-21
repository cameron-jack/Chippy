#!/usr/bin/env python
"""creates the gene feature database - must be run first"""
import sys
sys.path.extend(['..', '../src'])

import os
from cogent.db.ensembl import HostAccount
from chippy.util.util import create_path

from chippy.express.db_schema import make_session
from chippy.express.db_populate import add_ensembl_gene_data, \
        create_dummy_expr
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012"
__credits__ = ["Gavin Huttley, Cameron Jack, Anuj Pahwa"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"

def set_environment():
    script_info = {}
    script_info['title'] = 'Creates a chippy project'
    script_info['script_description'] = "Makes a chippy SQLite database."
    script_info['version'] = __version__

    # Process command-line arguments
    req_args = ['save_db_path', 'ensembl_release', 'species',
            'hostname', 'username', 'password', 'port']
    inputs = Args(required_args=req_args)

    return inputs.parsed_args, script_info

def main():
    args, script_info = set_environment()
    rr = RunRecord()
    
    if args.species != 'mouse':
        raise RuntimeError('Currently only support mouse, sorry!')

    create_path(args.save_db_path)

    if not os.path.isdir(args.save_db_path):
        sys.stderr.write('The save_db_path must be a directory.\n')
        return

    chippy_db_name = 'chippy_%d_%s.db' % (args.ensembl_release, args.species)
    db_path = os.path.join(args.save_db_path, chippy_db_name)
    session = make_session("sqlite:///%s" % db_path)
    hostname = args.hostname
    username = args.username
    password = args.password
    
    account = HostAccount(hostname, username, password, port=args.port)
    #add_ensembl_gene_data(session, args.species,
    #        ensembl_release=args.ensembl_release, account=account)

    success, rr = create_dummy_expr(session, rr=rr)
    if success:
        print 'Dummy data added successfully'
    else:
        print 'Dummy data failed to upload to DB. Expect bigger problems'

    rr.addInfo('start_chippy_db' ,'Chippy DB written to:', db_path)
    rr.display()

if __name__ == "__main__":
    main()
