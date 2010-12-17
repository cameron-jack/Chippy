import numpy as np
from cogent import LoadTable

def default_index_gen(length):
    data = tuple(range(length))
    def gen(i):
        temp = list(data)
        temp.pop(i)
        return temp
    return gen

class JackknifeStats(object):
    """Computes the jackknife statistic for a particular statistical function
    as outlined by 'Tukey's Jackknife Method' Biometry by Sokal/Rohlf."""

    def __init__(self, length, callback_stat,
                 callback_indexGen=default_index_gen):

        """Initialise the jackknife class:

        length: The length of the data set (since data is not passed to this
                class).
        callback_stat: A callback function that computes the required statistic
                       of a defined dataset.
        callback_indexGen: A callback function that generates a list of indices
                           that are used to sub-sample the dataset."""

        super(JackknifeStats, self).__init__()
        self.n = length
        self.n_minus_1 = self.n-1
        self.calc_stat = callback_stat
        self.indexGen = callback_indexGen(self.n)
        self._subset_statistics = list()
        self._pseudovalues = list()
        self._run = False


    def jackknife(self):
        """Computes the jackknife statistics and standard error"""

        # compute the statistic in question on the whole data set
        self.sample_statistic = self.calc_stat(range(self.n))
        n_sample_statistic = self.n*self.sample_statistic

        # compute the jackknife static for the data by removing an element in
        # each iteration and computing the statistic.
        for index in range(self.n):
            self._subset_statistics.append(self.calc_stat(self.indexGen(index)))
            pseudovalue = n_sample_statistic -\
                        (self.n_minus_1)*self._subset_statistics[index]
            self._pseudovalues.append(pseudovalue)

        self._pseudovalues = np.array(self._pseudovalues)
        self._subset_statistics = np.array(self._subset_statistics)
        self._jackknifed_stat = self._pseudovalues.mean(axis=0)

        # Compute the approximate standard error of the jackknifed estimate
        # of the statistic
        variance = np.square(self._pseudovalues - self._jackknifed_stat).sum(axis=0)
        variance_norm = np.divide(variance, self.n*self.n_minus_1)
        self._standard_error = np.sqrt(variance_norm)
        self._run = True

    @property
    def JackknifedStatistic(self):
        if self._run == False:
            self.jackknife()
        return self._jackknifed_stat

    @property
    def StandardStatError(self):
        if self._run == False:
            self.jackknife()
        return self._standard_error

    @property
    def SubSampleStatistics(self):
        """Return a table of the sub-sample statistics"""

        # if the statistics haven't been run yet.
        if self._run == False:
            self.jackknife()

        # generate table
        title = 'Subsample Statistics'
        rows = []
        for index in range(self.n):
            row = []
            row.append(index)
            subset_statistics = self._subset_statistics[index]
            try:
                for value in subset_statistics:
                    row.append(value)
            except TypeError:
                row.append(subset_statistics)
            rows.append(row)

        header = ['i']
        subset_stats = self._subset_statistics[0]

        try:
            num_datasets = len(subset_stats)
            for i in range(num_datasets):
                header.append('Statistic_%s-i'%i)
        except TypeError:
            header.append('Statistic-i')

        return LoadTable(rows=rows, header=header,title=title)

    @property
    def Pseudovalues(self):
        """Return a table of the Pseudovalues"""

        # if the statistics haven't been run yet.
        if self._run == False:
            self.jackknife()

        # detailed table
        title = 'Pseudovalues'
        rows = []
        for index in range(self.n):
            row = []
            row.append(index)
            pseudovalues = self._pseudovalues[index]
            try:
                for value in pseudovalues:
                    row.append(value)
            except TypeError:
                row.append(pseudovalues)
            rows.append(row)

        header = ['i']
        pseudovalues = self._pseudovalues[0]

        try:
            num_datasets = len(pseudovalues)
            for i in range(num_datasets):
                header.append('Pseudovalue_%s-i'%i)
        except TypeError:
            header.append('Pseudovalue-i')

        return LoadTable(rows=rows, header=header,title=title)

    @property
    def SummaryStatistics(self):
        """Return a summary table with the statistic value(s) calculated for the
        the full data-set, the jackknife statistics and standard errors."""

        # if the statistics haven't been run yet.
        if self._run == False:
            self.jackknife()

        header = ['Total Statistic', 'Jackknife Statistic', 'Standard Error']
        title = 'Summary of Statistics'
        rows = np.vstack((self.sample_statistic, self._jackknifed_stat, self._standard_error))
        rows = rows.transpose()
        return LoadTable(rows=rows, header=header,title=title)
