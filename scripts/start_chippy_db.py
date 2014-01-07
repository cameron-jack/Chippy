#!/usr/bin/env python
'''creates the gene feature database - must be run first'''
import sys
sys.path.extend(['..', '../src'])

import os
from cogent.db.ensembl import HostAccount
from chippy.util.util import create_path

from chippy.express.db_query import make_session
from chippy.express.db_populate import add_ensembl_gene_data, \
        create_dummy_expr
from chippy.util.run_record import RunRecord
from chippy.util.command_args import Args

__author__ = 'Gavin Huttley, Cameron Jack'
__copyright__ = 'Copyright 2011-2014, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack', 'Anuj Pahwa']
__license__ = 'GPL'
__version__ = '0.9.dev'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'Pre-release'

script_info = {}
script_info['title'] = 'Creates a chippy project'
script_info['script_description'] = 'Makes a chippy SQLite database.'
script_info['brief_description'] =\
    'Makes a ChipPy SQLite DB from an Ensembl source'
script_info['output_description'] = 'chippy sqlite database .db file'
script_info['version'] = __version__

# Process command-line arguments
req_args = ['save_db_dir', 'save_db_prefix', 'ensembl_release', 'species',
            'hostname', 'username', 'password', 'port']
opt_args = ['dummy_data']

script_info['args'] = Args(required_args=req_args, optional_args=opt_args)
script_info['required_options'] = script_info['args'].getReqCogentOpts()
script_info['optional_options'] = script_info['args'].getOptCogentOpts()

def main():
    rr = RunRecord('start_chippy_db')
    rr.addCommands(sys.argv)

    args = script_info['args'].parse()
    create_path(args.save_db_dir)

    if not os.path.isdir(args.save_db_dir):
        sys.stderr.write('The save_db_dir must be an existing directory.\n')
        return

    release = args.ensembl_release
    species = args.species
    chippy_db_name = args.save_db_prefix + '_chippy_' + str(release) +\
            '_' + species + '.db'
    db_path = os.path.join(args.save_db_dir, chippy_db_name)
    if not os.path.exists(db_path):
        session = make_session(db_path)

        hostname = args.hostname
        username = args.username
        password = args.password
    
        account = HostAccount(hostname, username, password, port=args.port)
        add_ensembl_gene_data(session, args.species,
                ensembl_release=args.ensembl_release, account=account)

        if args.dummy_data:
            success = create_dummy_expr(session)
            if success:
                print 'Dummy data added successfully. Expr=1.'
            else:
                print 'Dummy data failed to upload to DB. Expect bigger problems'

        rr.addInfo('Chippy DB written', db_path)
    else:
        rr.addError('Chippy DB with this name already exists', db_path)
    rr.display()

if __name__ == '__main__':
    main()
