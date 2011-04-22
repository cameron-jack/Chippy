#!/usr/bin/env python
import sys, os

from optparse import make_option
from sqlalchemy import and_
from cogent.util.misc import parse_command_line_parameters

from chippy.express.db_schema import Expression, ReferenceFile, make_session
from chippy.express.db_query import get_gene_expression_query, get_samples

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

if 'CHIPPY_DB' in os.environ:
    db_path = os.environ['CHIPPY_DB']
else:
    raise RuntimeError('You need to set an environment variable CHIPPY_DB '\
                       'that indicates where to find the database')

session = make_session( "sqlite:///%s" % db_path)

script_info = {}


script_info['title'] = 'Remove an expression study'
script_info['script_description'] = "Remove an expression study and all "\
                                    "associated linked objects."

script_info['version'] = __version__

_samples = get_samples(session)
choice_to_reffile = {}
if _samples:
    choices = []
    for s in _samples:
        name = s.name
        for reffile in s.reference_files:
            basename = os.path.basename(reffile.name)
            display = '%s - %s' % (name, basename)
            choices.append(display)
            choice_to_reffile[display] = reffile

if not choice_to_reffile:
    choices = [None]
    choice_to_reffile[None] = None

existing_sample = make_option('-s','--sample_reffile', type='choice',
    choices=choices, help='Select an sample+reffile combo to drop')

script_info['required_options'] = [existing_sample]

def _get_name_description(value):
    """returns the name and description of a : separated sample"""
    new_sample = value.split(':')
    assert len(new_sample) == 2, "Only one ':' allowed!"
    name, description = new_sample
    name = name.strip()
    description = description.strip()
    return name, description

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    reffile = choice_to_reffile[opts.sample_reffile]
    records = session.query(Expression).join(ReferenceFile).filter(
        and_(Expression.sample_id==reffile.sample_id,
             Expression.reffile_id==reffile.reffile_id)).all()
    for record in records:
        session.delete(record)
    session.delete(reffile)
    session.commit()
    

if __name__ == "__main__":
    main()
