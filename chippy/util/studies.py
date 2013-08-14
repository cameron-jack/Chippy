from __future__ import division

###
# studies.py
# Defines combinations of sequence collection and expression study.
#
# CentredStudy has been created by Export_Centred_Counts.py
#
# MatchedStudy holds ExprGene,ChrmGene lists.
#
# _Gene, _ExprGene and _ChrmGene are defined here as abstract genes,
# and simplified representations of processed expression and mapped
# sequence data respectively.
##

import os, sys, math
sys.path.extend(['..', '../src'])

import numpy
import gzip

from chippy.express import db_query
from chippy.util.run_record import RunRecord
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev
from chippy.express.db_query import make_session, get_gene_ids
from chippy.draw.plot_data import PlotLine

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2013, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.1'

class _Gene(object):
    """ defined by a stableId in a given study """

    def __init__(self, stableId, study, *args, **kwargs):
        super(_Gene, self).__init__()
        self.stableId = stableId
        self.study = study

    def __repr__(self):
        return repr((self.stableId, self.study))

class _ChrmGene(_Gene):
    """ gene entry from a ChipPy study """
    def __init__(self, counts, *args, **kwargs):
        self.counts = counts # a numpy array
        self.feature_pos = len(self.counts)/2
        self.feature_count = self.counts[self.feature_pos]
        self.total_count = numpy.sum(self.counts)
        super(_ChrmGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr((self.counts, self.rank, self.feature_pos,
                     self.feature_count))

class _ExprGene(_Gene):
    """ gene entry from a microarray expression study """
    def __init__(self, expr, *args, **kwargs):
        self.expr = expr
        super(_ExprGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr(self.expr, self.rank)

class MatchedStudy(object):
    """ stores a expr_study_gene_list, collection_gene_list pair
    """

    def __init__(self, expr_study_name, collection_name, db_path):
        self.chrm_gene_list = self._load_chrm(collection_name)
        self.expr_gene_list = self._load_expr(expr_study_name, db_path)

    def _load_expr(self, expr_study, db_path):
        """ loads expression records from a ChippyDB """
        rr = RunRecord('load_expr')

        sample_name = expr_study.split(' : ')[0]
        session = db_query.make_session(db_path)

        expr_gene_list = []
        #sample_type == 'Expression data: absolute ranked'
        print 'Querying sample from ChippyDB'
        sample_genes = db_query.get_ranked_expression(session, sample_name,
            biotype='protein_coding', data_path=None, rank_by='mean',
            test_run=False)
        for gene in sample_genes:
            gene_record = _ExprGene(gene.MeanScore, gene.ensembl_id, sample_name)
            expr_gene_list.append(gene_record)
        rr.addInfo('genes found in ' + sample_name, len(sample_genes))

        return expr_gene_list

    def _load_chrm(self, collection):
        """ loads gene entries from a ChipPy collection """
        rr = RunRecord('load_chrm')

        chrm_gene_list = []
        if os.path.isfile(collection):
            try:
                # to load counts data from file
                file1 = gzip.GzipFile(collection, 'rb')
                data = numpy.load(file1)
                d = data.tolist()
                counts = d['counts']
                labels = d['labels']

                for count, label in zip(counts, labels):
                    gene_record = _ChrmGene(count, str(label), collection)
                    chrm_gene_list.append(gene_record)
                rr.addInfo('genes found in ' + collection, len(labels))

            except IOError: # some exception type
                rr.dieOnCritical('file found but could not be read', collection)
        else:
            rr.dieOnCritical('unrecognised collection file', collection)

        return chrm_gene_list

    def keepCommonGenes(self):
        """ keep only those genes that are common to each study pair """
        rr = RunRecord('keep_common_genes')

        # get the intersection of all available stableIds
        kept = 0
        removed = 0

        for study_pair in self.matched_studies:
            expr_id_set = set(gene.stableId for gene in study_pair[0])
            chrm_id_set = set(gene.stableId for gene in study_pair[1])
            common_id_set = expr_id_set.intersection(chrm_id_set)

            for gene in study_pair[0]:
                if gene.stableId not in common_id_set:
                    study_pair[0].remove(gene)
                    removed += 1
                else:
                    kept += 1
            for gene in study_pair[1]:
                if gene.stableId not in common_id_set:
                    study_pair[1].remove(gene)
                    removed += 1
                else:
                    kept += 1

        rr.addInfo('number of genes kept', kept)
        rr.addInfo('number of genes discarded', removed)

class CentredStudy(object):
    """ Specifies the RegionCollection associated with an expression
            data set.
    Members: collection (a RegionCollection), window_radius,
            collection_label
    Methods: filterByGenes, filterByCutoff, normaliseByBases,
            asPlotLines
    """
    def __init__(self, collection_fn, counts_func,
                 *args, **kwargs):
        super(CentredStudy, self).__init__(*args, **kwargs)
        rr = RunRecord('Study')
        # Keep the source file name for labelling purposes
        fn = collection_fn.split('/')[-1].rstrip('.gz')
        self.collection_label = fn.replace('_', ' ')
        try:
            self.data_collection = RegionCollection(filename=collection_fn)
        except IOError:
            rr.dieOnCritical('Collection will not load', collection_fn)

        # Frequency normalized counts need to be converted
        if counts_func is column_sum:
            self.data_collection = self.data_collection.asfreqs()
        self.counts_func = counts_func

        # Get feature window radius
        try:
            self.window_radius =\
                    self.data_collection.info['args']['window_radius']
        except KeyError:
            self.window_radius = len(self.data_collection.counts[0])/2

    def filterByGenes(self, db_path, chrom=None, include_sample=None,
                      exclude_sample = None):
        """ keep only results that match selected genes """
        rr = RunRecord('filterByGenes')
        if not include_sample and not exclude_sample:
            return

        rr.addInfo('Starting no. of genes', self.data_collection.N)

        session = make_session(db_path)
        if include_sample:
            include_sample = include_sample.split(':')[0].strip()
        if include_sample:
            exclude_sample = exclude_sample.split(':')[0].strip()

        filter_gene_ids = get_gene_ids(session, chrom=chrom,
            include_target=include_sample, exclude_target=exclude_sample)

        self.data_collection =\
                self.data_collection.filteredByLabel(filter_gene_ids)
        rr.addInfo('Remaining genes', self.data_collection.N)

        if self.data_collection is None or\
           self.data_collection.ranks.max() == 0:
            rr.dieOnCritical('No genes remaining after filtering', 'Failure' )

        # total_features used to normalise coloring
        total_features = self.data_collection.ranks.max()
        self.data_collection.ranks /= total_features

    def filterByCutoff(self, cutoff=None):
        """ keep only results that pass Chebyshev cutoff """
        rr = RunRecord('filterByCutoff')

        rr.addInfo('Starting no. of genes', self.data_collection.N)

        # exclude outlier genes using one-sided Chebyshev
        if cutoff is not None and cutoff != 0.0:
            try:
                cutoff = float(cutoff)
                if cutoff < 0.0 or cutoff >= 1.0:
                    rr.addError('Cutoff out of range', cutoff)
                    rr.addInfo('Cutoff set to default', 0.05)
                    cutoff = 0.05
            except ValueError:
                rr.addError('Cutoff not given as float', cutoff)
                rr.addInfo('Cutoff set to default', 0.05)
                cutoff = 0.05
                # Do Chebyshev filtering

            self.data_collection =\
                    self.data_collection.filteredChebyshevUpper(p=cutoff)
            rr.addInfo('Used Chebyshev filter cutoff', cutoff)
            rr.addInfo('No. genes after normalisation filter',
                    self.data_collection.N)
        else:
            rr.addInfo('Outlier cutoff filtering', 'Off')

        if self.data_collection is None or\
                self.data_collection.ranks.max() == 0:
            rr.dieOnCritical('No data after filtering', 'Failure')

        # total_features used to normalise coloring
        total_features = self.data_collection.ranks.max()
        self.data_collection.ranks /= total_features

    def normaliseByRPM(self):
        """ This requires 'mapped tags', 'tag count' or 'base count' to be present
            in the collection and gives counts per mapped million tags/bases.
            Mapped tags is the total experimental mapped tags.
            Tag count and base count are region specific.
        """
        rr = RunRecord('normaliseByRPM')
        try:
            norm_RPM = self.data_collection.info['args']['mapped tags']
            rr.addInfo("'mapped tags' value", norm_RPM)
        except KeyError:
            rr.addError('Info field not found', 'mapped tags')
            return
        norm_factor = 1000000.0/norm_RPM
        rr.addInfo('normalising by RPMs', norm_factor)
        normalised_counts = []
        for c in self.data_collection.counts:
            c = c * norm_factor
            normalised_counts.append(c)
        self.data_collection.counts = numpy.array(normalised_counts)

    def _groupAllGeneCounts(self):
        """ Group counts for all genes and return as a single PlotLine.
            Called by asPlotLines or _groupNGeneCounts().
            Returns a list.
        """
        rr = RunRecord('_groupAllGeneCounts')
        counts, ranks, se = self.data_collection.transformed(\
            counts_func=self.counts_func)
        if not len(counts):
            rr.dieOnCritical('No counts data in', 'Study.groupAllGeneCounts')

        # Always name single lines by their collection name
        label = self.collection_label
        plot_lines = [PlotLine(counts, ranks, label, study=label, stderr=se)]
        return plot_lines

    def _groupNoGeneCounts(self):
        """ Don't group counts. Simply return a PlotLine for each set of
            counts.
            Called by asPlotLines()
        """
        rr = RunRecord('_groupNoGeneCounts')
        counts = self.data_collection.counts
        ranks = self.data_collection.ranks
        labels = self.data_collection.labels
        plot_lines = []
        for c,r,l in zip(counts, ranks, labels):
            if self.counts_func == stdev:
                stdev_ = c.std()
                if stdev_ > 0:
                    c = (c - c.mean()) / stdev_
                    plot_lines.append(PlotLine(c, r , l,
                        study=self.collection_label))
            else:
                plot_lines.append(PlotLine(c, r , l,
                    study=self.collection_label))

        # If no data was returned default to groupAllCollectionCounts
        if not len(plot_lines):
            rr.dieOnCritical('No data in collection', 'Failure')

        # If a single line is created label it with the collection name
        if len(plot_lines) == 1:
            plot_lines[0].label = [self.collection_label]

        return plot_lines

    def _groupNGeneCounts(self, group_size):
        """ Group counts for N genes and return as PlotLines. Defaults to
            _groupAllGeneCounts() if group size is too large.
            Called by asPlotLines()
        """
        rr = RunRecord('_groupNGeneCounts')
        plot_lines = []
        for index, (c,r,l,se) in enumerate(self.data_collection.\
        iterTransformedGroups(group_size=group_size,
            counts_func=self.counts_func)):
            plot_lines.append(PlotLine(c, rank=index, label=l,
                study=self.collection_label, stderr=se))

        # If no data was returned default to groupAllCollectionCounts
        if not len(plot_lines):
            rr.addWarning('Defaulting to ALL features. Not enough '+\
                          'features for group of size', group_size)
            plotLines = self.groupAllGeneCounts()
            return plotLines

        # If a single line is created label it with the collection name
        if len(plot_lines) == 1:
            plot_lines[0].label = [self.collection_label]

        return plot_lines

    def asPlotLines(self, studies, group_size, group_location):
        """ returns a list of PlotLine objects from this study """
        rr = RunRecord('asPlotLines')

        if type(group_size) is str and group_size.lower() == 'all':
            plot_lines= self._groupAllGeneCounts()
        elif type(group_size) is int:
            if group_size == 1:
                plot_lines = self._groupNoGeneCounts()
            else:
                plot_lines = self._groupNGeneCounts(group_size)
        else:
            rr.dieOnCritical('group_size, wrong type or value',
                [type(group_size), group_size])

        if group_location:
            rr.addInfo('grouping genes from location', group_location)
            if group_location.lower() == 'top':
                plot_lines = [plot_lines[0]]
            elif group_location.lower() == 'middle':
                plot_lines = [plot_lines[len(plot_lines)/2]]
            elif group_location.lower() == 'bottom':
                plot_lines = [plot_lines[-1]]

        rr.addInfo('Plottable lines from study', len(plot_lines))
        return plot_lines

