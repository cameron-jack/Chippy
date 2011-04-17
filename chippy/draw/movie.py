import sys, glob, os
from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

try:
    from AppKit import NSImage
    from QTKit import QTMovie, QTMakeTime, QTAddImageCodecType
except ImportError:
    raise NotImplementedError("Currently restricted to OSX systems")


__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2011, Anuj Pahwa, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "alpha"
__version__ = '0.1'

def make_movie(image_dir, mov_name, suffix, time_val=1, time_scale=1):
    """I am creating a new image sequence (QuickTime movie).
    time_scale is spf or fps"""
    
    if not os.path.exists(image_dir):
        raise RuntimeError('Could not find image_dir: %s' % image_dir)
    
    mov, err = QTMovie.alloc().initToWritableFile_error_(mov_name, None)
    
    imgpaths = glob.glob1(image_dir, suffix)
    
    if mov == None:
        raise IOError("Could not create movie file: %s" % (mov_name))
    
    duration = QTMakeTime(time_val, time_scale)
    # you can also use "tiff"
    attrs = {QTAddImageCodecType: "tiff"}
    
    for imgpath in imgpaths:
        img = NSImage.alloc().initWithContentsOfFile_(os.path.join(image_dir,
                                                        imgpath))
        mov.addImage_forDuration_withAttributes_(img, duration, attrs)
        mov.updateMovieFile()
    


script_info = {}
descr = "Convert directory of images into a Quicktime movie"
script_info['script_description'] = descr
script_info['version'] = __version__

script_info['required_options'] = [
    make_option('-i','--image_dir',
            help='Directory containing image files'),
    make_option('-m','--movie_name',
                help='Name of the movie file'),
    make_option('-s','--image_suffix', type='choice',
                choices=['pdf', 'png'],
                help='Filename suffix of image files'),
    ]

script_info['optional_options'] = [\
    make_option('-T','--time_value', type='int', default=1,
                help='Time value'),
    make_option('-S','--time_scale', type='int', default=1,
                help='Time scale, in frames per second'),
                ]

def main():
    parser, opts, args = parse_command_line_parameters(**script_info)
    make_movie(opts.image_dir, opts.movie_name, '*.%s' % opts.image_suffix,
        opts.time_value, opts.time_scale)

if __name__ == "__main__":
    main()

