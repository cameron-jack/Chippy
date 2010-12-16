import sys
sys.path.append('../src')

import numpy as np

from cogent.util.unit_test import TestCase, main

from jackknife import JackknifeStats

def pmcc(data):
    """Compute the Product-moment correlation coefficient.
    Expression 15.3 from Biometry by Sokal/Rohlf
    This code implementation is on the proviso that the data that is provided
    is two dimensional: [[Y1], [Y2]] (trying to determine the correlation
    coefficient between data sets Y1 and Y2"""

    mean = data.mean(axis=1)
    data_less_mean = np.array([data[0] - mean[0],
                              data[1] - mean[1]])
    sum_squares = np.sum(np.square(data_less_mean), axis=1)
    sum_products = np.sum(np.prod(data_less_mean, axis=0))
    pmcc = np.divide(sum_products, np.sqrt(np.prod(sum_squares)))
    z_trans = np.arctanh(pmcc)
    return z_trans

# test data from Box 15.2; Biometry by Sokal/Rohlf
data = np.array([[159, 179, 100, 45, 384, 230, 100, 320, 80, 220, 320, 210],
                [14.40, 15.20, 11.30, 2.50, 22.70, 14.90, 1.41, 15.81, 4.19, 15.39, 17.25, 9.52]])

class JackknifeTests(TestCase):

    def test_proper_initialise(self):
        test_knife = JackknifeStats(data=data,axis=1,statistic=pmcc)
        self.assertEqual(test_knife.data, data)
        self.assertEqual(test_knife.n, data.shape[1])
        self.assertEqual(len(test_knife.pseudovalues), data.shape[1])
        self.assertEqual(len(test_knife.subset_statistic), data.shape[1])
        self.assertEqual(test_knife.run, False)

    def test_jackknife(self):
        test_knife = JackknifeStats(data=data,axis=1,statistic=pmcc)
        self.assertAlmostEquals(test_knife.jackknife(), 1.2905845)
        self.assertAlmostEquals(test_knife.standardError(), 0.2884490)
        self.assertEqual(test_knife.run, True)

    def test_jackknife_table(self):
        test_knife = JackknifeStats(data=data,axis=1,statistic=pmcc)
        table = test_knife.observationTable()
        self.assertEqual(test_knife.run, True)

if __name__ == "__main__":
    main()