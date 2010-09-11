import numpy as np
from cogent import LoadTable

class JackknifeStats(object):
    """Computes the jackknife statistic for a particular statistical function
    as outlined by 'Tukey's Jackknife Method' Biometry by Sokal/Rohlf."""

    def __init__(self, data, axis, statistic):

        """Initialise the class:
        data = set of data you want to perform the statistic on
        axis = axis of correlation. If axis is 0, statistic is computed across
               the row; If axis is 1 statistic is computed accross columns.
               Example 1 (Correlation between x and y; axis 0):
               [[x1, y1]
                [x2, y2]
                [x3, y3]
                [x4, y4]]
               Example 2 (Correlation between x and y; axis 1):
               [[x1, x2, x3, x4]
                [y1, y2, y3, y4]]

        statistic = the function that will compute the statistic on the data
                    provided. The statistical value returned should be a
                    value"""

        super(JackknifeStats, self).__init__()
        self.data = data
        self.axis = axis
        self.n = data.shape[axis] # sample size depends on axis.
        self.n_minus_1 = self.n-1
        self.pseudovalues = np.zeros(self.n)
        self.subset_statistic = np.zeros(self.n)
        self.statistic = statistic
        self.jackknifed_stat = 0.0
        self.standard_error = 0.0
        self.run = False


    def jackknife(self):
        """Performs the grunt of the statistical work and returns the
        jackknifed statistic"""

        # compute the statistic in question on the whole data set
        sample_statistic = self.statistic(self.data)
        n_sample_statistic = self.n*sample_statistic

        # compute the jackknife static for the data by removing an element in
        # each iteration and computing the statistic.
        for index in range(self.n):
            indices = range(self.n)
            indices.remove(index)
            subset_data = self.data.take(indices, self.axis)
            self.subset_statistic[index] = self.statistic(subset_data)
            self.pseudovalues[index] = n_sample_statistic -\
                (self.n_minus_1)*self.subset_statistic[index]

        self.jackknifed_stat = self.pseudovalues.mean()

        """Compute the approximate standard error of the jackknifed estimate
        of the statistic"""
        variance = np.sum(np.square(self.pseudovalues - self.jackknifed_stat))
        variance_norm = np.divide(variance, self.n*self.n_minus_1)
        self.standard_error = np.sqrt(variance_norm)

        self.run = True

        return self.jackknifed_stat

    def standardError(self):
        return self.standard_error

    def observationTable(self):
        """Create a summary table of the subset statistic and pseudovalues"""

        # if the statistics haven't been run yet.
        if self.run == False:
            self.jackknife()

        header = ['i', 'Statistic-i', 'Pseudovalue-i']
        title = 'Jackknife Stat: %f; Standard Error: %f' % \
                (self.jackknifed_stat, self.standard_error)
        rows = []
        for index in range(self.n):
            rows += [[index+1, self.subset_statistic[index],
                      self.pseudovalues[index]]]

        return LoadTable(rows=rows, header=header, title=title)



