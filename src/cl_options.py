#!/usr/bin/env python

from __future__ import division

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2009, Modelling the impact of DNA methylation"
__credits__ = ["Gavin Huttley", "Hua Ying"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin Huttley"
__status__ = "Prototype"


from optparse import OptionParser

usage_str = """usage: %prog [options]

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
"""

def parse_command_line_parameters(add_options_func, required_options=[]):
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    add_options_func(parser)

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False)

    opts,args = parser.parse_args()

    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option)


    return opts,args


def config_options(parser):
    """sets up the command line options"""

    parser.add_option('-i', '--input_file',
       help="name of the input file")

    parser.add_option('-o', '--output_file',
       help="name of the output results file")


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose