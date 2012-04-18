from __future__ import division
from math import log10, floor, ceil

import os, sys, glob
sys.path.extend(['..', '../src'])

import numpy

from optparse import make_option
from cogent.util.misc import parse_command_line_parameters

from matplotlib import pyplot
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011, Gavin Huttley, Anuj Pahwa, Cameron Jack'
__credits__ = ['Gavin Huttley, Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Gavin Huttley'
__email__ = 'Gavin.Huttley@anu.edu.au'
__status__ = 'alpha'
__version__ = '0.1'


script_info = {}
script_info['title'] = 'Boxplot'
script_info['script_description'] = "Make a boxplot of expression ranks of a subset of genes from a larger grouping"
script_info['version'] = __version__
script_info['authors'] = __author__
script_info['output_description']= "boxplot."
script_info['help_on_no_arguments'] = True

# alternate option organisation

# essential source files
opt_subset = make_option('-s', '--subset',
    help='Path to the sub set data')
opt_collection = make_option('-c', '--collection',
    help='Path to the super set of data')

script_info['required_options'] = [opt_collection, opt_subset]

def main():
    option_parser, opts, args =\
    parse_command_line_parameters(**script_info)

    collection_file = opts.collection
    subset_file = opts.subset

    collection_name_bits = collection_file.split('/')
    collection_name = collection_name_bits[-1]
    subset_name_bits = subset_file.split('/')
    subset_name = subset_name_bits[-1]

    data_collection = RegionCollection(filename=collection_file)
    data_subset = RegionCollection(filename=subset_file)

    collection_labels_ranks = dict(zip(data_collection.labels, data_collection.ranks))
    subset_labels_ranks = dict(zip(data_subset.labels, data_subset.ranks))

    intersecting_genes = filter(collection_labels_ranks.has_key, subset_labels_ranks.keys())

    print "\nelements in collection: %d\n" % len(collection_labels_ranks)
    print "\nelements in sub set: %d\n" % len(subset_labels_ranks)
    print "\nelements in intersection: %d\n" % len(intersecting_genes)

    intersecting_ranks = []
    for gene in intersecting_genes:
        if collection_labels_ranks.has_key(gene):
            intersecting_ranks.append(collection_labels_ranks.pop(gene))

    rank_array = numpy.array(intersecting_ranks)
    pyplot.boxplot(rank_array)
    title = 'Gene ranks of %s in %s' % (subset_name, collection_name)
    pyplot.suptitle(title, fontsize=8)
    pyplot.ylabel('Gene rank. Lower=higher expr')
    pyplot.savefig('/home/cameron/boxplot.png')

if __name__ == '__main__':
    main()