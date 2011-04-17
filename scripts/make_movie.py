# because of a bug in PyCogent.app, we run movie making as a subprocess
import subprocess, os, warnings, sys
from optparse import make_option
import cogent
from cogent.util.misc import parse_command_line_parameters

from chippy.draw import movie
from chippy.util.util import make_cl_command, safe_save_path
from chippy.prep.command_line import run_command

__author__ = 'Gavin Huttley'
__copyright__ = 'Copyright 2011, Gavin Huttley'
__credits__ = ['Gavin Huttley']
__license__ = 'GPL'
__maintainer__ = 'Gavin Huttley'
__email__ = 'Gavin.Huttley@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'


script_info = movie.script_info
movie_path = movie.__file__
if movie_path.endswith('.pyc'):
    movie_path = movie_path[:-1]

cogent_path = os.path.dirname(os.path.dirname(cogent.__file__))


def run_movie_command(args):
    """runs the command to make a movie as an external process due to a bug
    in PyCogent.app"""
    command = ' '.join(['PYTHONPATH=%s python' % cogent_path, movie_path,
                        args])
    retcode, stdout, stderr = run_command(command, False)
    if retcode != 0:
        sys.stdout.writelines(stdout)
        sys.stderr.writelines(stderr)
    
    return

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    cl_args = script_info['command_line_args']
    
    # make sure a user has entered a safe path for movie_name
    movie_path_value_index = cl_args.index('--movie_name') + 1
    movie_path = safe_save_path(cl_args[movie_path_value_index])
    if movie_path != cl_args[movie_path_value_index]:
        warnings.warn('no absolute path specified, saving to %s' % movie_path,
        UserWarning, 2)
    
    cl_args[movie_path_value_index] = movie_path
    # make_cl_command
    
    run_movie_command(' '.join(cl_args))
    

if __name__ == "__main__":
    main()

