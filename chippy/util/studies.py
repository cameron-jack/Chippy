from __future__ import division

"""
 studies.py
 Defines combinations of sequence collection and expression study.

 CentredStudy has been created by Export_Centred_Counts.py

 MatchedStudy holds ExprGene,ChrmGene lists.

 _Gene, _ExprGene and _ChrmGene are defined here as abstract genes,
 and simplified representations of processed expression and mapped
 sequence data respectively.
"""

import os, sys
sys.path.extend(['..'])

import numpy
import gzip

from chippy.express import db_query
from chippy.util.run_record import RunRecord
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev
from chippy.express.db_query import make_session, get_gene_ids
from chippy.draw.plot_data import PlotLine, PlotPoint

__author__ = 'Cameron Jack'
__copyright__ = 'Copyright 2011-2014, Gavin Huttley, Cameron Jack, Anuj Pahwa'
__credits__ = ['Gavin Huttley', 'Cameron Jack']
__license__ = 'GPL'
__maintainer__ = 'Cameron Jack'
__email__ = 'cameron.jack@anu.edu.au'
__status__ = 'pre-release'
__version__ = '0.1'



class Gene(object):
    """ defined by a stableId in a given study """

    def __init__(self, stableId, study, *args, **kwargs):
        super(Gene, self).__init__()
        self.stableId = stableId
        self.study = study

    def __repr__(self):
        return repr((self.stableId, self.study))

class CountsGene(Gene):
    """ gene entry from a ChipPy study """
    Rank = 0 # low rank means higher score
    Score = 0 # whichever counts feature is chosen
    def __init__(self, counts, *args, **kwargs):
        self.counts = counts # a numpy array

        if 'window_upstream' in kwargs.keys() and\
                'window_downstream' in kwargs.keys():
            window_upstream = kwargs['window_upstream']
            window_downstream = kwargs['window_downstream']

            self.promoter_counts = numpy.sum(self.counts[:window_upstream])
            self.coding_counts = numpy.sum(self.counts[window_downstream:])
            self.feature_pos = window_upstream + 1
            self.feature_counts = self.counts[self.feature_pos]
        else:
            self.promoter_counts = numpy.sum(self.counts[:len(counts)/2])
            self.coding_counts = numpy.sum(self.counts[len(counts)/2:])
            self.feature_pos = len(self.counts)/2
            self.feature_counts = self.counts[self.feature_pos]

        self.total_counts = numpy.sum(self.counts)
        super(CountsGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr((self.counts, self.Rank, self.feature_pos,
                     self.feature_count))

class CountsStudy(object):
    """ A collection of CountsGene objects """
    def __init__(self, collection):
        """ loads gene entries from a ChipPy collection """
        rr = RunRecord('load_counts')

        print 'Loading counts collection file', collection
        self.counts_genes = []
        if os.path.isfile(collection):
            try:
                # to load counts data from file
                file1 = gzip.GzipFile(collection, 'rb')
                data = numpy.load(file1)
                d = data.tolist()
                info = d['info']['args']
                # get collection file metadata
                self.window_upstream = info.get('window_upstream', None)
                self.window_downstream = info.get('window_downstream', None)
                self.window_radius = info.get('window_radius', None)
                self.feature_type = info.get('feature_type', None)
                self.sample_name = info['sample_name']
                self.species = info['species']
                self.tag_count = info['tag count']
                self.base_count = info.get('base count', None)
                self.mapped_tags = info['mapped tags']

                counts = d['counts']
                labels = d['labels']

                for count, label in zip(counts, labels):
                    gene_record = CountsGene(count, str(label), collection,
                            self.window_upstream, self.window_downstream)
                    self.counts_genes.append(gene_record)
                rr.addInfo('genes found in ' + collection, len(labels))

            except IOError: # some exception type
                rr.dieOnCritical('file found but could not be read', collection)
        else:
            rr.dieOnCritical('unrecognised collection file', collection)

    def normaliseByRPM(self):
        """ Normalise counts by per million mapped tags """
        rr = RunRecord('normaliseByRPM')
        norm_factor = 1000000.0/self.mapped_tags
        rr.addInfo('normalising by RPMs', norm_factor)

        for gene in self.counts_genes:
            gene.promoter_counts *= norm_factor
            gene.coding_counts *= norm_factor
            gene.feature_counts *= norm_factor
            gene.total_counts *= norm_factor

    def scoresAsRankedArray(self, metric, log2=False):
        """
            sort the gene scores based on some metric and return ordered
            Numpy array. Return as the log base 2 scores if log2==True
        """
        if metric == 'feature':
            scores = sorted([gene.feature_counts \
                    for gene in self.counts_genes])
        elif metric == 'total':
            scores = sorted([gene.total_counts \
                    for gene in self.counts_genes])
        elif metric == 'promoter':
            scores = sorted([gene.promoter_counts \
                    for gene in self.counts_genes])
        elif metric == 'coding':
            scores = sorted([gene.coding_counts \
                    for gene in self.counts_genes])

        scores = numpy.array(scores)
        if log2:
            scores[scores > 0] = numpy.log2(scores[scores > 0])

        return scores

class ExprGene(Gene):
    """ gene entry from a microarray expression study """
    def __init__(self, score, rank, *args, **kwargs):
        self.Score = score
        self.Rank = rank
        super(ExprGene, self).__init__(*args, **kwargs)

    def __repr__(self):
        return repr((self.Score, self.Rank))

class ExprStudy(object):
    """ A collection of ExprGene objects """
    def __init__(self, sample_name, db_path, include_targets=None,
                  exclude_targets=None):
        """
            loads expression records from a ChippyDB and also
            ranks by expr
        """
        rr = RunRecord('load_expr')
        session = db_query.make_session(db_path)

        self.expr_genes = []
        #sample_type == 'Expression data: absolute ranked'
        print 'Querying sample from ChippyDB'

        sample_genes = db_query.get_genes_by_ranked_expr(session, sample_name,
            biotype='protein_coding', data_path=None, rank_by='mean',
            include_targets=include_targets, exclude_targets=exclude_targets)

        for gene in sample_genes:
            gene_record = ExprGene(gene.MeanScore, gene.Rank,
                gene.ensembl_id, sample_name)
            self.expr_genes.append(gene_record)
        rr.addInfo('genes found in ' + sample_name, len(sample_genes))

    def scoresAsRankedArray(self, log2=False):
        """
            sort the gene scores based on some metric and return ordered
            Numpy array. Return as the log base 2 scores if log2==True
        """
        scores = sorted([gene.Score for gene in self.expr_genes])

        scores = numpy.array(scores)
        if log2:
            scores[scores > 0] = numpy.log2(scores[scores > 0])

        return scores

class _MatchedGene(object):
    """ represents a single gene with both expression and counts properties """
    def __init__(self, counts_gene, expr_gene):
        self.counts_gene = counts_gene
        self.expr_gene = expr_gene
        self.id = counts_gene.stableId

class MatchedStudy(object):
    """ stores a expr_study_gene_list, collection_gene_list pair
    """

    def __init__(self, expr_study_name, collection_name, db_path,
            region_feature='total_counts', include_targets=None,
            exclude_targets=None):
        """
            Loads counts and expr. Calculates counts rank based on 
            region_feature. Creates a matched_gene list.
        """
        self.load_expr(expr_study_name, db_path,
                include_targets=include_targets,
                exclude_targets=exclude_targets)

        self.load_counts(collection_name)

        # Calculate collection count data ranks
        if region_feature.lower() == 'feature':
            ranked_counts = sorted(self.counts_genes,
                    key=lambda counts: counts.total_counts, reverse=True)
        elif region_feature.lower() == 'promoter':
            ranked_counts = sorted(self.counts_genes,
                    key=lambda counts: counts.total_counts, reverse=True)
        elif region_feature.lower() == 'coding':
            ranked_counts = sorted(self.counts_genes,
                    key=lambda counts: counts.total_counts, reverse=True)
        else: # total_counts
            ranked_counts = sorted(self.counts_genes,
                    key=lambda counts: counts.total_counts, reverse=True)

        self.counts_genes = []
        for i, counts in enumerate(ranked_counts):
            counts.Rank = i
            if region_feature.lower() == 'feature':
                counts.Score = counts.feature_counts
            elif region_feature.lower() == 'promoter':
                counts.Score = counts.promoter_counts
            elif region_feature.lower() == 'coding':
                counts.Score = counts.coding_counts
            else: # total_counts
                counts.Score = counts.total_counts
            self.counts_genes.append(counts)

        self.common_genes = self.keepCommonGenes()

    def load_expr(self, expr_study, db_path, include_targets=None,
            exclude_targets=None):
        """
            loads expression records from a ChippyDB and also
            ranks by expr
        """
        rr = RunRecord('load_expr')

        sample_name = expr_study.split(' : ')[0]
        session = db_query.make_session(db_path)

        self.expr_genes = []
        #sample_type == 'Expression data: absolute ranked'
        print 'Querying sample from ChippyDB', sample_name

        sample_genes = db_query.get_genes_by_ranked_expr(session, sample_name,
                biotype='protein_coding', data_path=None, rank_by='mean',
                include_targets=include_targets,
                exclude_targets=exclude_targets)

        for gene in sample_genes:
            gene_record = ExprGene(gene.MeanScore, gene.Rank,
                    gene.ensembl_id, sample_name)
            self.expr_genes.append(gene_record)
        rr.addInfo('genes found in ' + sample_name, len(sample_genes))

    def load_counts(self, collection):
        """ loads gene entries from a ChipPy collection """
        rr = RunRecord('load_counts')

        print 'Loading counts collection file', collection
        self.counts_genes = []
        if os.path.isfile(collection):
            try:
                # to load counts data from file
                file1 = gzip.GzipFile(collection, 'rb')
                data = numpy.load(file1)
                d = data.tolist()
                counts = d['counts']
                labels = d['labels']

                for count, label in zip(counts, labels):
                    gene_record = CountsGene(count, str(label), collection)
                    self.counts_genes.append(gene_record)
                rr.addInfo('genes found in ' + collection, len(labels))

            except IOError: # some exception type
                rr.dieOnCritical('file found but could not be read', collection)
        else:
            rr.dieOnCritical('unrecognised collection file', collection)

    def keepCommonGenes(self):
        """ keep only those genes that are common to each study pair """

        print 'Keeping common genes'

        expr_id_obj = {}
        counts_id_obj = {}

        # build common id list
        expr_id_set = set([gene.stableId for gene in self.expr_genes])
        counts_id_set = set([gene.stableId for gene in self.counts_genes])
        common_id_set = expr_id_set.intersection(counts_id_set)

        for counts in self.counts_genes:
            if counts.stableId in common_id_set:
                counts_id_obj[counts.stableId] = counts

        for expr in self.expr_genes:
            if expr.stableId in common_id_set:
                expr_id_obj[expr.stableId] = expr

        self.matched_genes = []

        for id in common_id_set:
            mg = _MatchedGene(counts_id_obj[id], expr_id_obj[id])
            self.matched_genes.append(mg)

    def get_matched_genes_as_xy_plotpoints(self, x_axis_type, expr_ranks=False,
            counts_ranks=False):
        """ return a set of points ready for plotting """

        if x_axis_type.lower() == 'expression':
            if expr_ranks and counts_ranks:
                return [PlotPoint(gene.expr_gene.Rank, gene.counts_gene.Rank) \
                        for gene in self.matched_genes]
            elif expr_ranks:
                return [PlotPoint(gene.expr_gene.Rank, gene.counts_gene.Score)\
                        for gene in self.matched_genes]
            elif counts_ranks:
                return [PlotPoint(gene.expr_gene.Score, gene.counts_gene.Rank)\
                        for gene in self.matched_genes]
            else:
                return [PlotPoint(gene.expr_gene.Score, gene.counts_gene.Score)\
                        for gene in self.matched_genes]

        else: # counts is X-axis
            if expr_ranks and counts_ranks:
                return [PlotPoint(gene.counts_gene.Rank, gene.expr_gene.Rank)\
                        for gene in self.matched_genes]
            elif expr_ranks:
                return [PlotPoint(gene.counts_gene.Score, gene.expr_gene.Rank)\
                        for gene in self.matched_genes]
            elif counts_ranks:
                return [PlotPoint(gene.counts_gene.Rank, gene.expr_gene.Score)\
                        for gene in self.matched_genes]
            else:
                return [PlotPoint(gene.counts_gene.Score, gene.expr_gene.Score)\
                        for gene in self.matched_genes]


class RegionStudy(object):
    """ Specifies the RegionCollection associated with an expression
            data set. Used to collate data for plot_counts.py.

    Members: collection (a RegionCollection), window_start, window_end,
            collection_label
    Methods: filterByGenes, filterByCutoff, normaliseByBases,
            asPlotLines
    """
    def __init__(self, collection_fn, counts_func,
                 *args, **kwargs):
        super(RegionStudy, self).__init__(*args, **kwargs)
        rr = RunRecord('Study')
        # Keep the source file name for labelling purposes
        self.collection_path = collection_fn
        fn = collection_fn.split('/')[-1].rstrip('.gz')
        self.collection_label = ' '.join(fn.replace('_', ' ').split('.')[:-1])
        try:
            self.data_collection = RegionCollection(filename=collection_fn)
        except IOError:
            rr.dieOnCritical('Collection will not load', collection_fn)

        # Frequency normalized counts need to be converted
        if counts_func is column_sum:
            self.data_collection = self.data_collection.asfreqs()
        self.counts_func = counts_func

        # Get feature window start and end
        try:
            self.window_upstream =\
                    self.data_collection.info['args']['window_upstream']
        except KeyError:
            rr.dieOnCritical('Collection value not defined', 'window_upstream')

        try:
            self.window_downstream =\
                    self.data_collection.info['args']['window_downstream']
        except KeyError:
            rr.dieOnCritical('Collection value not defined', 'window_downstream')

        try:
            self.feature_type =\
                    self.data_collection.info['args']['feature_type']
        except KeyError:
            self.feature_type = 'Unknown'

    def filterByGenes(self, db_path, chrom=None, include_samples=None,
                      exclude_samples=None):
        """ keep only results that match selected genes """

        rr = RunRecord('filterByGenes')
        if not include_samples and not exclude_samples and not chrom:
            return

        rr.addInfo('Starting no. of genes', self.data_collection.N)

        session = make_session(db_path)
        if include_samples:
            for sample in include_samples:
                rr.addInfo('Restricting plot by include sample', sample)

        if exclude_samples:
            for sample in exclude_samples:
                rr.addInfo('Restricting plot by exclude sample', sample)

        if not chrom is None:
            rr.addInfo('Restricting plot to chromosome', chrom)

        filter_gene_ids = get_gene_ids(session, chrom=chrom,
                include_targets=include_samples,
                exclude_targets=exclude_samples)

        self.data_collection =\
                self.data_collection.filteredByLabel(filter_gene_ids)
        rr.addInfo('Remaining genes', self.data_collection.N)

        if self.data_collection is None or\
           len(self.data_collection.ranks) == 0:
            rr.dieOnCritical('Genes remaining after filtering', '0' )

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
            c2 =  c * norm_factor
            normalised_counts.append(c2)
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
            rr.dieOnCritical('No counts data in', 'Study._groupAllGeneCounts')

        ranks = 0 # rank is irrelevant for 'all' genes

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

    def _groupNGeneCounts(self, group_size, p=None):
        """ Group counts for N genes and return as PlotLines. Defaults to
            _groupAllGeneCounts() if group size is too large.
            Called by asPlotLines()
        """
        rr = RunRecord('_groupNGeneCounts')
        plot_lines = []
        for index, (c,r,l,se) in enumerate(self.data_collection.\
                iterTransformedGroups(group_size=group_size,
                counts_func=self.counts_func, p=p)):
            plot_lines.append(PlotLine(c, rank=r, label=l,
                study=self.collection_label, stderr=se))

        # If no data was returned default to groupAllCollectionCounts
        if not len(plot_lines):
            rr.addWarning('Defaulting to ALL features. Not enough '+\
                          'features for group of size', group_size)
            plotLines = self._groupAllGeneCounts()
            return plotLines

        return plot_lines

    def asPlotLines(self, group_size, group_location, p=None):
        """
            Returns a list of PlotLine objects from this study.
            'p' is the Chebyshev cut-off if not None
        """
        rr = RunRecord('asPlotLines')
        if p is not None:
            rr.addInfo('Applying per-line Chebyshev filtering', p)

        if type(group_size) is str and group_size.lower() == 'all':
            plot_lines= self._groupAllGeneCounts()
        elif type(group_size) is int:
            if group_size == 1:
                plot_lines = self._groupNoGeneCounts()
            else:
                plot_lines = self._groupNGeneCounts(group_size,p=p)
        else:
            rr.dieOnCritical('group_size, wrong type or value',
                [type(group_size), group_size])

        if group_location.lower() != 'all':
            rr.addInfo('grouping genes from location', group_location)
            plot_lines.sort(key=lambda x: x.rank)
            if group_location.lower() == 'top':
                plot_lines = [plot_lines[0]]
            elif group_location.lower() == 'middle':
                plot_lines = [plot_lines[int(len(plot_lines)/2)]]
            elif group_location.lower() == 'bottom':
                plot_lines = [plot_lines[-1]]

        rr.addInfo('Plottable lines from study', len(plot_lines))
        return plot_lines

