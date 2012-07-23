#!/usr/bin/env python
import sys, os

from optparse import make_option
from sqlalchemy import and_
from cogent.util.misc import parse_command_line_parameters

from chippy.express.db_schema import Expression, ExpressionDiff, ReferenceFile, \
    make_session
from chippy.express.db_query import get_gene_expression_query, get_samples
from chippy.util.run_record import RunRecord
from chippy.util.definition import LOG_DEBUG, LOG_INFO, LOG_WARNING, \
    LOG_ERROR, LOG_CRITICAL

__author__ = "Gavin Huttley, Cameron Jack"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley, Cameron Jack"
__credits__ = ["Gavin Huttley, Cameron Jack"]
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
NO_FILES = 'No files'

if _samples:
    choices = []
    for s in _samples:
        name = s.name
        this_sample = []
        for reffile in s.reference_files:
            basename = os.path.basename(reffile.name)
            display = '%s - %s' % (name, basename)
            this_sample.append(display)
            choice_to_reffile[display] = reffile
        if not this_sample:
            display = '%s - %s' % (name, NO_FILES)
            this_sample = [display]
            choice_to_reffile[display] = s
        
        choices.extend(this_sample)

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
    
    rr = RunRecord()
    object = choice_to_reffile[opts.sample_reffile]
    if isinstance(object, ReferenceFile):
        reffile = object
        records = session.query(Expression).join(ReferenceFile).filter(
            and_(Expression.sample_id==reffile.sample_id,
                 Expression.reffile_id==reffile.reffile_id)).all()
        if len(records) == 0:
            records = session.query(ExpressionDiff).join(ReferenceFile).filter(
            and_(ExpressionDiff.sample_id==reffile.sample_id,
                 ExpressionDiff.reffile_id==reffile.reffile_id)).all()
        rr.addMessage('drop_expression_db', LOG_INFO,
            'No. dropped expression records', len(records))
        for record in records:
            session.delete(record)
        # and if this is the last file with the sample, delete it
        sample = reffile.sample
        if len(sample.reference_files) < 2:
            rr.addMessage('drop_expression_db', LOG_INFO,
                'Dropped sample', str(sample))
            session.delete(sample)
        
        session.delete(reffile)
        rr.addMessage('drop_expression_db', LOG_INFO,
            'Dropped reffile basename', os.path.basename(reffile.name))
    else:
        rr.addMessage('drop_expression_db', LOG_INFO,
            'Dropped sample', str(object))
        session.delete(object) # a sample
    session.commit()
    rr.display()

if __name__ == "__main__":
    main()

