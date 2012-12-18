#!/usr/bin/env python
"""creates the gene feature database - must be run first"""
import sys
sys.path.extend(['..', '../src'])

import os
from optparse import make_option

from cogent import LoadTable
from cogent.db.ensembl import HostAccount
from cogent.util.misc import parse_command_line_parameters
from chippy.util.util import create_path

from chippy.express.db_schema import make_session
from chippy.express.db_populate import add_ensembl_gene_data, upload_data
from chippy.express.db_query import get_stable_id_genes_mapping
from chippy.express.util import sample_types
from chippy.util.run_record import RunRecord

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2012"
__credits__ = ["Gavin Huttley, Cameron Jack, Anuj Pahwa"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Cameron Jack"
__email__ = "cameron.jack@anu.edu.au"
__status__ = "Pre-release"

script_info = {}

script_info['title'] = 'Creates a chippy project'
script_info['script_description'] = "Makes a chippy SQLite database."
script_info['version'] = __version__

script_info['required_options'] = [
    make_option('-S','--save_db_path',
        help='path to directory where chippy.db will be saved.'),
    make_option('-R','--ensembl_release', type='int',
        help='Ensembl release to use.'),
    make_option('-s','--species', type='choice', default='mouse',
        choices=['mouse', 'human'],
        help='Create for species'),
]

_mysql_msg = '(uses ENSEMBL_ACCOUNT if not provided)'
script_info['optional_options'] = [
    make_option('--hostname', default=None,
            help='hostname for MySQL Ensembl server %s' % _mysql_msg),
    make_option('--username', default=None,
            help='username MySQL Ensembl server %s' % _mysql_msg),
    make_option('--password', default=None,
            help='password for MySQL Ensembl server %s' % _mysql_msg),
    make_option('--port', default=None,
        help='Port for MySQL Ensembl server')
]

def create_dummy_expr(session, rr=RunRecord()):
    """ create flat and spread dummy data """
    genes_dict = get_stable_id_genes_mapping(session)

    # flat expression dummy
    expr_table_rows = []
    expr_table_rows.append(['gene', 'probeset', 'exp']) # header
    for i, gene_id in enumerate(genes_dict):
        expr_table_rows.append(['gene_id', 'P'+str(i), '1'])

    success, rr = upload_data(session, 'dummy_flat',
            'each gene has expression score of 1',
            'None', expr_table_rows, gene_id_heading='gene',
            probeset_heading='probeset', expr_heading='exp',
            sample_type=sample_types['exp_absolute'],
            reffile1=None, reffile2=None, rr=rr)
    if not success:
        return success, rr

    # spread expression dummy
    expr_table_rows = []
    expr_table_rows.append(['gene', 'probeset', 'exp']) # header
    for i, gene_id in enumerate(genes_dict):
        expr_table_rows.append(['gene_id', 'P'+str(i), i])

    success, rr = upload_data(session, 'dummy_spread',
            'each gene has unique expression score',
            'None', expr_table_rows, gene_id_heading='gene',
            probeset_heading='probeset', expr_heading='exp',
            sample_type=sample_types['exp_absolute'],
            reffile1=None, reffile2=None, rr=rr)
    if not success:
        return success, rr

    return success, rr

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)
    rr= RunRecord()
    
    if opts.species != 'mouse':
        raise RuntimeError('Currently only support mouse, sorry!')

    create_path(opts.save_db_path)

    if not os.path.isdir(opts.save_db_path):
        sys.stderr.write('The save_db_path must be a directory.\n')
        return

    chippy_db_name = 'chippy_%d_%s.db' % (opts.ensembl_release, opts.species)
    db_path = os.path.join(opts.save_db_path, chippy_db_name)
    session = make_session("sqlite:///%s" % db_path)
    hostname = opts.hostname
    username = opts.username 
    password = opts.password
    if None in (hostname, username, password):
        try:
            hostname, username, password = os.environ['ENSEMBL_ACCOUNT'].split()
        except KeyError:
            sys.stderr.write('Must provide account settings or have valid '\
                    'ENSEMBL_ACCOUNT environment variable\n')
            return
    
    account = HostAccount(hostname, username, password, port=opts.port)
    add_ensembl_gene_data(session, opts.species,
            ensembl_release=opts.ensembl_release, account=account)

    success, rr = create_dummy_expr(session, rr=rr)
    if success:
        print 'Dummy data added successfully'
    else:
        print 'Dummy data failed to upload to DB. Expect bigger problems'

    rr.addInfo('start_chippy_db' ,'Chippy DB written to:', db_path)
    rr.display()

if __name__ == "__main__":
    main()
