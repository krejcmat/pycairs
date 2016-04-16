__author__ = 'matt'
import sys
from math import sqrt, exp


class CovExponentional(object):
    def __init__(self, sigma=10, l_spatial=3000, l_temporal=600, gamma=1):
        """
         separable gamma-exponential

        ## Equation (4.18) in Rasmussen, C. E. and Williams, C. K. I. (2006) Gaussian processes
        ## for machine learning, Massachusett, MIT Press.
        @:param sigma -standard deviation of GP
        @:param l_spatioal -spatial correlation length
        @:param l_temporal -temporal correlation length [seconds]
        @:param gamma -exponent for smoothness in [0, 2]
        """
        self.sigma = sigma
        self.l_spatial = l_spatial
        self.l_temporal = l_temporal
        self.gamma = gamma
        if gamma < 0 or gamma > 2:
            sys.exit("Gamma must be in [0,2]")

            # TODO check if time is in milliseconds

    def __str__(self):
        str = ''
        str += "- Prior has a separable gamma-exponential covariance function with:\n"
        str += ("   standard deviation: %s\n" % self.sigma)
        str += ("   spatial correlation length: %s\n" % self.l_spatial)
        str += ("   temporal correlation length [s]: %s\n" % (self.l_temporal))
        str += ("   gamma: %s" % self.gamma)
        return str

    def f_cov(self, c1, c2):
        """
        :param c1: coordination
        :type c1: Coord
        :param c2: coordination
        :type c2: Coord
        :return: covariance between two points
        :rtype: number
        """
        # covariance
        var = self.sigma * self.sigma
        dist_spatial = sqrt((c1.x - c2.x) ** 2 + (c1.y - c2.y) ** 2)
        dist_temporal = abs(c1.time - c2.time)
        cov = var * exp(
            - (dist_spatial / self.l_spatial) ** self.gamma - (dist_temporal / self.l_spatial) ** self.gamma)
        return cov


class CovSphere(object):
    def __init__(self, sigma=10, l_spatial=3000, l_temporal=600):
        self.sigma = sigma
        self.l_spatial = l_spatial
        self.l_temporal = l_temporal
        self.var = sigma * sigma
        # check time if is in millisecond

    def __str__(self):
        str = ''
        str += ("- Prior has a spherical covariance function with:\n")
        str += ("   standard deviation:%s\n" % self.sigma)
        str += ("   spatial radius:%s\n" % self.l_spatial)
        str += ("   temporal radius [s]:%s" % (self.l_temporal))
        return str

    def f_cov(self, c1, c2):
        """
        :param c1: Coor
        :type c1:
        :param c2: Coor
        :type c2:
        :return: covariance between c1 a c2
        :rtype:
        """

        r = sqrt(((c1.x - c2.x) / self.l_spatial) ** 2 + ((c1.y - c2.y) / self.l_spatial) ** 2 + (
            (c1.time - c2.time) / self.l_temporal) ** 2)
        if r > 1:
            cov = 0
        else:
            cov = self.var * (1 - 3 / 2 * r + +1 / 2 * r ** 3)
        return cov




