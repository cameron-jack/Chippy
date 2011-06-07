#!/usr/bin/env python
"""creates the gene feature database"""
import sys, os
from optparse import make_option

from cogent import LoadTable
from cogent.db.ensembl import HostAccount
from cogent.util.misc import parse_command_line_parameters

from chippy.express.db_schema import make_session
from chippy.express.db_populate import add_ensembl_gene_data

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.9.dev"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

script_info = {}

script_info['title'] = 'Creates a chippy project'
script_info['script_description'] = "Makes a chippy SQLite database."
script_info['version'] = __version__

script_info['required_options'] = [
    make_option('-S','--save_db_path',
        help='path to directory where chippy.db will be saved.'),
    make_option('-R','--ensembl_release',
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

script_info['authors'] = __author__

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)
    
    if opts.species != 'mouse':
        raise RuntimeError('Currently only support mouse, sorry!')
    
    if not os.path.isdir(opts.save_db_path):
        sys.stderr.write('The save_db_path must be a directory.\n')
        return
    
    db_path = os.path.join(opts.save_db_path, 'chippy.db')
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
    

if __name__ == "__main__":
    main()
