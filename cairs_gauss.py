import numpy
import sys

# sys.path.remove('/usr/local/lib/python2.7/dist-packages/cubature-0.13.0_dev-py2.7-linux-x86_64.egg')
#print sys.path
#sys.path.insert(1,
#                'home/matt/apt/envVirtual/lib/python2.7/site-packages/cubature/dist/cubature-0.13.0_dev-py2.7-linux-x86_64.egg')
#sys.path.insert(1, '/home/matt/apt/envVirtual/lib/python2.7/site-packages/cubature')

from cubature import cubature
from math import sqrt, exp
from cairs_coor import Coor, Domain

# from datetime import timedelta
#import pdb
import logging

logger = logging.getLogger('cairs')


class Mean(object):
    def __init__(self, mean=0):
        self.mean = mean
        print('Prior has a constant mean of %s' % mean)
        self.extend_v = None
        self.pos_v = None
        self.d = None
        self.counter = 0
        #st=StaticContext()
        #self.REF_TIME=st.REF_TIME

    def f_mu(self, d):
        """
        #vector of domain extend
        mean function
        :param d: domain
        :type d: Domain
        """

        if type(d) is Domain:
            self.d = d
            self.extend_v = numpy.array([d.extend.x,
                                         d.extend.y,
                                         d.extend.time], dtype=numpy.float64)
            # vector of positions
            self.pos_v = numpy.array([d.position.x,
                                      d.position.y,
                                      d.position.time], dtype=numpy.float64)
            lower = []

            for p, v in zip(self.pos_v, self.extend_v):
                if v != 0:
                    lower.append(p)
            #lower=min(lower)

            upper = []
            for p, v in zip(self.pos_v, self.extend_v):
                if v != 0:
                    if v + p != 0:
                        x = p + v
                        upper.append(x)

            upper = numpy.array(upper, dtype=numpy.float64)
            lower = numpy.array(lower, dtype=numpy.float64)
            #I = integrate.fixed_quad(self.f_int,a=lower[0],b=upper[0])

            I = cubature(self.f_int, ndim=lower.shape[0],
                         fdim=upper.shape[0], maxEval=10000,
                         xmin=lower, xmax=upper, adaptive='p',
            )


            return I[0][0]
        else:
            return self.mean

    def f_int(self, x_array):
        """
        :brief : construct function to integrate over
        :param v: numbers
        :type v: list
        :return:
        :rtype:
        """
        kk = 0

        for i in range(0, 2, 1):
            if self.extend_v[i] != 0:
                self.pos_v[i] = x_array[kk]
                kk += 1

        coorTMP = Coor(self.pos_v[0], self.pos_v[1], self.pos_v[2])
        coorTMP = Coor.rotate(coorTMP, self.d.position, self.d.angle)

        return self.f_mu(coorTMP)


class Covarience(object):
    def __init__(self, sigma=10, cov_o=None, l_spatial=3000, l_temporal=60, gamma=1):

        self.sigma = sigma
        self.cov_o = cov_o

        self.extend_v = None
        self.pos_v = None
        self.d1 = None
        self.d2 = None

        self.l_spatial = l_spatial
        self.l_temporal = l_temporal
        self.gamma = gamma
        if gamma < 0 or gamma > 2:
            sys.exit("Gamma must be in [0,2]")

    def f_cov(self, d1, d2):
        """
        :param d1:
        :type d1:
        :param d2:
        :type d2:
        :return:
        :rtype:
        """
        #logger.debug("d1 = %s"%d1)
        #logger.debug("d2 = %s"%d2)
        mark = False
        if type(d1) is Coor and type(d2) is Coor:
            #logger.debug("self.c1 %s" % d1)
            #logger.debug("self.c2 %s" % d2)
            var = self.sigma ** 2
            dist_spatial = sqrt((d1.x - d2.x) ** 2 + (d1.y - d2.y) ** 2)
            dist_temporal = abs(d1.time - d2.time)
            timeDiv = dist_temporal / self.l_temporal
            expFeed = (-1) * ((dist_spatial / self.l_spatial) ** self.gamma) - (timeDiv ** self.gamma)
            cov = var * exp(expFeed)
            #logger.debug('134')
            #logger.debug("f_cov %s" % cov)
            return cov
        else:
            # -----------
            ## covariance (kernel) function for ( Domain, Domain )
            if type(d1) is Coor and type(d2) is Domain:  ## --- for ( Coor, Domain )
                d1 = Domain(d1, Coor(0, 0), 0)
                mark = True
            if type(d2) is Coor and type(d1) is Domain and not mark:
                d2 = Domain(d2, Coor(0, 0), 0)



            self.extend_v = numpy.array([d1.extend.x,
                                         d1.extend.y,
                                         d1.extend.time,
                                         d2.extend.x,
                                         d2.extend.y,
                                         d2.extend.time], dtype=numpy.float64)

            # # vector of positions
            self.pos_v = numpy.array([d1.position.x,
                                      d1.position.y,
                                      d1.position.time,
                                      d2.position.x,
                                      d2.position.y,
                                      d2.position.time], dtype=numpy.float64)
            self.d1 = d1
            self.d2 = d2

            # # construct function to integrate over
            lower = list()

            for p, v in zip(self.pos_v, self.extend_v):
                if v != 0:
                    lower.append(p)

            upper = list()
            for p, v in zip(self.pos_v, self.extend_v):
                if v != 0 and v + p != 0:
                    x = p + v
                    upper.append(x)
            ## integrate
            #sumExtend_v=sum(self.extend_v)
            #TODOif sumExtend_v != 0 and sumExtend_v >2:          # for >2-dimensional integration
            dim = 0
            for ex in self.extend_v:
                if ex != 0:
                    dim += 1

            upper = numpy.array(upper, dtype=numpy.float64)
            lower = numpy.array(lower, dtype=numpy.float64)
            if dim > 2:
                #logger.debug("dim >2")
                #    "I = cubature(self.f_int, xmin=lower,xmax=upper,ndim=lower.shape[0],fdim=upper.shape[0], maxEval=10000)")
                #I = integrate.nquad(self.f_int,[lower,upper])
                I = cubature(self.f_int, xmin=lower, xmax=upper, ndim=lower.shape[0], fdim=upper.shape[0],
                             maxEval=10000)
                #I = integrate.fixed_quad(self.f_int,a=lower[0],b=upper[0])
            else:
                #logger.debug("dim <= 2")

                #logger.debug("self.d1 %s" % self.d1)
                #logger.debug("self.d2 %s" % self.d2)
                #logger.debug("upper = %s" % upper)
                #logger.debug("lower = %s" % lower)
                #logger.debug('199')
                #if len(upper)>2:

                #I = integrate.quad(self.f_int,a=lower[0],b=upper[0])
                #logger.debug('upper.shape %s' % upper.shape)
                #.debug('upper.shape[0] %s' % upper.shape[0])

                I = cubature(self.f_int, ndim=lower.shape[0],
                             fdim=upper.shape[0], maxEval=10000,
                             xmin=lower, xmax=upper, adaptive='p')

                #logger.error(logger.findCaller())
                #logger.debug('212 cubature --------')


            #logger.debug("I[0] = %s" % I[0])
            try:
                return I[0][0]
            except:
                return I[0]
        logger.error('f_cov nothing happened')

    def f_int(self, x_array):
        #logger.debug('f_int')
        # # change the coordiates for integration
        kk = 0
        #logger.debug('x_array %s' % numpy.asarray(x_array))
        for i in range(0, 5, 1):
            if self.extend_v[i] != 0:
                if type(x_array) in [float, numpy.float64]:
                    #logger.error("x_array = %s" % x_array)
                    self.pos_v[i] = x_array
                else:
                    self.pos_v[i] = x_array[kk]
                #break

                kk += 1
        #logger.debug("len(self.pos_v) = %s"%len(self.pos_v))
        #logger.error(logger.findCaller())
        x = Coor.rotate(Coor(self.pos_v[0], self.pos_v[1], self.pos_v[2]),
                        self.d1.position,
                        self.d1.angle)
        #logger.error('x = Coor.rotate')
        y = Coor.rotate(Coor(self.pos_v[3], self.pos_v[4], self.pos_v[5]),
                        self.d2.position,
                        self.d2.angle)
        #logger.error(logger.findCaller())
        #logger.debug('242')

        cov = self.f_cov(x, y)

        r = list()
        if len(numpy.asarray(x_array)) == 1:
            return cov
        else:
            for i in range(0, len(numpy.asarray(x_array)), 1):
                r.append(cov)
        return r


    def make_cov(self, loc_1, loc_2):
        """
        :param loc_1:
        :type loc_1:
        :param loc_2:
        :type loc_2:
        :return:
        :rtype:
        """

        if loc_1 != loc_2:
            sigma = numpy.zeros(shape=(len(loc_1), len(loc_2)), dtype=numpy.float64)
            for i,l1 in enumerate(loc_1):
                for j,l2 in enumerate(loc_2):
                    sigma[i,j]=(self.f_cov(l1, l2))  # TODO check

            #logger.debug("sigma makexcov = %s"%sigma)
        else:
            n = len(loc_1)
            sigma = numpy.zeros(shape=(n, n), dtype=numpy.float64)
            for j in range(0, n, 1):
                for i in range(j, n, 1):
                    if i != j:
                        x = self.f_cov(loc_1[i], loc_2[j])
                        sigma[i, j] = sigma[j, i] = x
                    else:
                        x = self.f_cov(loc_1[i], loc_2[j])
                        sigma[i, j] = x
        return sigma

    @staticmethod
    def log_p_prior(locations, samp_dict, i_sample, mu, sigma):
        R = list()
        for loc in locations:
            a = samp_dict[loc]
            R.append(a[i_sample])
        D = list()
        for i, m in enumerate(mu):
            a = R[i] - m
            D.append(a)
        #logger.debug("R = %s" % R)
        #logger.debug("D = %s" % D)
        # compute D' inv(Sigma) * D
        neg_log_p = numpy.dot(numpy.transpose(D), numpy.linalg.inv(sigma))
        neg_log_p = numpy.dot(neg_log_p, D)  #SOLVED

        return -neg_log_p



