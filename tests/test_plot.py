import sys
sys.path.extend(['../src', '..'])

import numpy

from scripts.plot_counts import load_studies, div_plots, set_plot_colors
from chippy.draw.plot_data import PlotLine
from chippy.draw.plottable import PlottableGroups
from chippy.core.collection import RegionCollection, column_sum, column_mean, stdev

from cogent.util.unit_test import TestCase, main

class PlotCountTests(TestCase):
    def generate_expected_counts(self):
        """
            These are the functions used to generate the original data files.
        """
        pos_c = numpy.zeros([200,100], dtype=numpy.float32)
        for i in xrange(100):
            for j in xrange(200):
                if j < 100:
                    pos_c[i][j] = 100 + i + j
                else:
                    pos_c[i][j] = 200 + i - (100 - j)

        neg_c = numpy.zeros([200,100], dtype=numpy.float32)
        for i in xrange(100):
            for j in xrange(200):
                neg_c[i][j] = 400 - i - j
        return pos_c, neg_c

    def test_load_studies(self):
        # Test loading one sample with mean counts

        studies, window_upstream, window_downstream =\
                load_studies(['data/plot_data_pos_corr.chp'], column_mean)

        self.assertEqual(len(studies), 1)
        self.assertEqual(window_upstream, 100)
        self.assertEqual(window_downstream, 100)
        self.assertEqual(studies[0].collection_label, 'plot data pos corr')

        # Test loading two samples with mean counts
        studies, window_upstream, window_downstream =\
                load_studies(['data/plot_data_pos_corr.chp',
                              'data/plot_data_neg_corr.chp'], column_mean)

        self.assertEqual(len(studies), 2)
        self.assertEqual(window_upstream, 100)
        self.assertEqual(window_downstream, 100)
        self.assertEqual(studies[0].collection_label, 'plot data pos corr')
        self.assertEqual(studies[1].collection_label, 'plot data neg corr')

    def test_div_plots(self):
        """
            Divide a dataset by itself to get unity
        """
        sample_studies, window_up, window_down =\
                load_studies(['data/plot_data_pos_corr.chp'], column_mean)
        div_studies, window_up, window_down =\
                load_studies(['data/plot_data_pos_corr.chp'], column_mean)
        div_studies[0].collection_label += '_div'

        group_size = 10
        group_location = 'all'
        # Need to combine studies so we get one set of plot lines
        # They need to have different names for div to work
        study_lines = sample_studies[0].asPlotLines(group_size, group_location, p=0.0)
        div_lines = div_studies[0].asPlotLines(group_size, group_location, p=0.0)
        lines = []
        for s,d in zip(study_lines, div_lines):
            lines.append(s)
            lines.append(d)
        div_plots(lines, div_studies[0].collection_label, div_by='all')
        pass

    def test_plot_colors(self):
        """
            Check that colours are assigned in the correct direction,
            dependant upon line rank and not upon counts.
            red ~= rank 1, blue ~= rank -1
            The color space is not spaced linearly so only relative color
            can be checked.
        """
        # positive correlation between expression and counts
        pos_studies, window_upstream, window_downstream =\
                load_studies(['data/plot_data_pos_corr.chp'], column_mean)
        group_size = 10
        group_location = 'all'
        study_lines = pos_studies[0].asPlotLines(group_size, group_location, p=0.0)

        lines = set_plot_colors(study_lines, pos_studies, None, 'black', False,
                            restrict_colors=None)

        # bluest should have lowest expression (highest rank)
        rgb_lines = sorted(lines, key=lambda l: l.color) # blue to red
        self.assertTrue(rgb_lines[0].rank > rgb_lines[-1].rank)
        # and lowest counts
        self.assertTrue(rgb_lines[0].getMaxCount() < rgb_lines[-1].getMaxCount())

        # negative correlation between expression and counts
        neg_studies, window_upstream, window_downstream =\
                load_studies(['data/plot_data_neg_corr.chp'], column_mean)
        group_size = 10
        group_location = 'all'
        study_lines = neg_studies[0].asPlotLines(group_size, group_location, p=0.0)

        lines = set_plot_colors(study_lines, neg_studies, None, 'black', False,
                restrict_colors=None)

        # bluest should have lowest expression (highest rank)
        rgb_lines = sorted(lines, key=lambda l: l.color) # blue to red
        self.assertTrue(rgb_lines[0].rank > rgb_lines[-1].rank)
        # and highest counts
        self.assertTrue(rgb_lines[0].getMaxCount() > rgb_lines[-1].getMaxCount())

    def test_plot_colors_multistudy(self):
        """
            Check that assigned line colours are correct when multiple
            studies are plotted together.
            There should be ten lines for each study in this test
            Each study is assigned a different color
            Base color should be the same for all lines from the same study
            Only alpha should change with rank
        """
        # Test loading two samples with mean counts
        studies, window_upstream, window_downstream =\
                load_studies(['data/plot_data_pos_corr.chp',
                'data/plot_data_neg_corr.chp'], column_mean)

        group_size = 10
        group_location = 'all'
        study_lines = []
        for study in studies:
            lines = study.asPlotLines(group_size, group_location, p=0.0)
            for l in lines:
                study_lines.append(l)

        lines = set_plot_colors(study_lines, studies, None, 'black', False,
                            restrict_colors=None)
        groups = {}
        for line in lines:
            if line.study not in groups.keys():
                groups[line.study] = []
            groups[line.study].append(line)

        # There should be 10 lines in each study group
        for g in groups.keys():
            self.assertEqual(len(groups[g]), 10)

        # each study has a different color
        k1, k2 = groups.keys()
        self.assertTrue(groups[k1][0].color != groups[k2][0].color)

        # Each group should have its own color
        # and alpha should be different
        for g in groups.keys():
            group_color = groups[g][0].color
            group_rgb = (group_color[0], group_color[1], group_color[2])
            set_alpha = set()
            for l in groups[g]:
                self.assertFloatEqual(group_rgb[0], l.color[0])
                self.assertFloatEqual(group_rgb[1], l.color[1])
                self.assertFloatEqual(group_rgb[2], l.color[2])
                set_alpha.add(l.alpha)
            self.assertEqual(len(set_alpha), 10)

        # check that alpha falls with increasing rank
        for g in groups.keys():
            ranked_lines = sorted([l for l in groups[g]], key=lambda x: x.rank)
            for i in range(len(ranked_lines)):
                if i == 0:
                    continue
                self.assertTrue(ranked_lines[i].alpha < ranked_lines[i-1].alpha)

class LineTests(TestCase):
    def simple_lines(self):
        """ return a list of simple plottable lines """
        counts_flat = numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        counts_jiggle = numpy.array([0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0])
        counts_ascending = numpy.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        counts_descending = numpy.array([8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0])
        lines = []
        lines.append(PlotLine(counts_flat, rank=1))
        lines.append(PlotLine(counts_jiggle, rank=1))
        lines.append(PlotLine(counts_ascending, rank=1))
        lines.append(PlotLine(counts_descending, rank=1))
        return lines

    def test_smoothing(self):
        """ Smooth some plot lines """
        plot_lines = self.simple_lines()
        smoothing = 4
        for line in plot_lines:
            line.applySmoothing(smoothing)
        expected_counts = []
        expected_counts.append(numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
        expected_counts.append(numpy.array([0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
        expected_counts.append(numpy.array([0.0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]))
        expected_counts.append(numpy.array([9.0, 7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 1.5]))

        for pl, el in zip(plot_lines, expected_counts):
            self.assertFloatEqual(pl.counts, el)

    def test_binning(self):
        """ Bin some plot lines """
        plot_lines = self.simple_lines()
        binning = 4
        for line in plot_lines:
            line.applyBinning(binning)
        expected_counts = []
        expected_counts.append(numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
        expected_counts.append(numpy.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]))
        expected_counts.append(numpy.array([2.5, 2.5, 2.5, 2.5, 6.5, 6.5, 6.5, 6.5]))
        expected_counts.append(numpy.array([6.5, 6.5, 6.5, 6.5, 2.5, 2.5, 2.5, 2.5]))

        for pl, el in zip(plot_lines, expected_counts):
            self.assertFloatEqual(pl.counts, el)

    def test_getMaxCounts(self):
        plot_lines = self.simple_lines()
        expected = [1.0, 2.0, 8.0, 8.0]
        expected2 = [3.0, 4.0, 10.0, 10.0]
        for i, line in enumerate(plot_lines):
            max_ = line.getMaxCount()
            max_plus_stderr = line.getMaxCount(include_stderr=True)
            self.assertFloatEqual(max_, expected[i])

            # no line.stderr
            self.assertFloatEqual(max_plus_stderr, expected[i])

            # Now add an artificial stderr
            line.stderr = 2.0
            max_plus_stderr = line.getMaxCount(include_stderr=True)
            self.assertFloatEqual(max_plus_stderr, expected2[i])

    def test_getMinCounts(self):
        plot_lines = self.simple_lines()
        expected = [1.0, 0.0, 1.0, 1.0]
        expected2 = [-1.0, -2.0, -1.0, -1.0]
        for i, line in enumerate(plot_lines):
            min_ = line.getMinCount()
            min_minus_stderr = line.getMinCount(include_stderr=True)
            self.assertFloatEqual(min_, expected[i])

            # no line.stderr
            self.assertFloatEqual(min_minus_stderr, expected[i])

            # Now add an artificial stderr
            line.stderr = 2.0
            min_minus_stderr = line.getMinCount(include_stderr=True)
            self.assertFloatEqual(min_minus_stderr, expected2[i])

class PlottableTests(TestCase):
    """
        TODO: tests for auto-figuring of axes and grid lines
    """
    def test_axes(self):
        pass

    def test_check_y_axis_scale(self):
        pass

    def test_auto_y_lims(self):
        pass

    def test_auto_grid(self):
        pass

if __name__ == "__main__":
    main()
