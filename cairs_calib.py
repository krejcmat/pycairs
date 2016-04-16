# def make_function_dict()
import numpy
from math import sqrt, exp
from sets import Set
import math
import logging
import random
from cairs_coor import Domain,Coor
from cairs_gauss import Covarience
from cairs_helpers import shaw


logger = logging.getLogger('cairs')


class Gibbs(object):
    def __init__(self, signals, prior_mean, prior_cov, n_samples, burn_in, adaption=True):
        """
        Gibbs sampling with adaptive jump distributions

        Roberts G, Rosenthal J (2009). "Examples of Adaptive MCMC."
        Computational Statistics and Data Analysis, 18, 349-367.

        signals:      vector of 'Signals'
        prior_mean:     mean function of Prior, f(c::Coor)
        prior_cov:      covariance function of Prior, f(c1::Coor, c2::Coor)
        n_sample:     number of MCMC sample
        burn_in:      number of samples to remove as burn-in
        adaption:     true/false

        FUTURE: - Ordered Overrelaxation (see MacKay)?
             - adaptive rejection sampling (http://www.stat.duke.edu/~cnk/Links/slides.pdf)?
             - slice sampling?
        """
        logger.info('Module - %s, Class - %s ' % (__name__, self.__class__.__name__))
        # numpy.set_printoptions(formatter={'float': lambda x: 'float: ' + str(x)})
        self.signals = signals  # object Signal
        self.o_mu = prior_mean  # object Mean
        self.o_cov = prior_cov  # object CovExponentional
        self.n_samples = n_samples
        self.burn_in = burn_in
        self.adaption = adaption
        #self.runGibbs()

    def runGibbs(self):
        # # -----------
        ## 1) set-up

        ## create overlaoded function of Prior
        #f_mu, f_cov = overload_GP_function(prior_mean, prior_cov)
        #showSignals(self.signals,'Signals')
        ##  create dictionary of all signals
        signal_dict = self.make_signal_dict(self.signals)
        #showSignals(signal_dict,'signal_dict')
        ## create dictionary of all point rains to sample from
        samples_dict = self.make_sample_point_dict(signal_dict, self.n_samples)
        #showSignals(samples_dict,'samples_dict')

        samples_dict_prop = self.make_sample_point_dict(signal_dict, 1)
        #showSignals(samples_dict_prop,"samples_dict_prop")
        ## array with all coordiates

        samp_points = []
        for k in samples_dict.keys():
            samp_points.append(k)
        #showSignals(samp_points,'samp_points')
        ## Compute mean

        mu = []
        for loc in samp_points:
            mu.append(self.o_mu.f_mu(loc))
        #showSignals(mu,'mu')
        #logger.info('mu = \n%s\n' % mu)
        Sigma = self.o_cov.make_cov(samp_points, samp_points)  #works

        logger.debug('Sigma = \n%s\n' % Sigma)
        numpy.set_printoptions(suppress=True)
        #cov1=numpy.sqrt(cov)
        #cov_prior = numpy.diag(cov1)#works
        #logger.info('cov_prior = \n%s\n' % cov_prior)
        ## compute inverse covariance matrix
        #Sigma = numpy.linalg.cholesky(numpy.transpose(cov))
        #logger.info('cholesky = \n%s\n' % Sigma)

        sd_prior = (numpy.diag(numpy.sqrt(Sigma)))  #works
        #logger.info('sd_prior = \n%s\n' % sd_prior)
        #sys.exit()

        sd_prop = {}
        for i, val in enumerate(samples_dict):
            sd_prop[val] = sd_prior[i]

        for i in range(0, self.n_samples - 1, 1):
            for location in samples_dict.keys():
                #location=kk
                # find all signals that condition on 'location'
                signals = signal_dict[location]
                #show(signals, 'signals123')
                #logger.debug('*' * 30)
                #logger.debug('signals = \n%s\n' % str(signals[0]))
                # proposal
                #random.random() *
                samples_dict_prop[location][0] = samples_dict[location][i] +  sd_prop[location] * numpy.random.standard_normal(1)[0]
                #showSignals(samples_dict_prop,'samples_dict_prop[location][0]')
                #logger.debug('samples_dict_prop[location][0] = %s' % samples_dict_prop[location][0])
                log_p_prop = Covarience.log_p_prior(samp_points, samples_dict_prop, 0, mu, Sigma)  #solved
                #logger.debug('log_p_prop1 = %s' % log_p_prop)

                for S in signals:
                    a = self.log_p_of_signal(S, samples_dict_prop, 0)
                    #logger.debug("log_p_of_signal1 = %s" % a)
                    log_p_prop += a
                #logger.debug('log_p_prop = %s' % log_p_prop)

                log_p_old = Covarience.log_p_prior(samp_points, samples_dict, i, mu, Sigma)
                #logger.debug('log_p_old1 = %s' % log_p_old)
                for S in signals:
                    a = self.log_p_of_signal(S, samples_dict, i)
                    #logger.debug("log_p_of_signal2 = %s" % a)
                    log_p_old += a
                #logger.debug('log_p_old = %s' % log_p_old)

                #prvar(log_p_prop)
                aa = log_p_prop - log_p_old
                #logger.debug('log_p_prop-log_p_old = %s' % aa)
                #prvar(aa)
                try:
                    ex = math.exp(log_p_prop - log_p_old)  # if overflowed >1
                    p_acc = min(1, ex)
                except:
                    p_acc = 1

                #logger.debug('p_acc = %s' % p_acc)
                #prvar(p_acc)
                #random.random()
                if  random.random()< p_acc:
                    samples_dict[location][i + 1] = samples_dict_prop[location][0]
                else:
                    samples_dict[location][i + 1] = samples_dict[location][i]
                    samples_dict_prop[location][0] = samples_dict[location][i]

                ## adapt jump distribution, see Roberts G, Rosenthal J (2009), Sec. 3
                ## (adaption is only during burn-in active)
                if self.adaption and i < self.burn_in and i % 50 == 0:

                    tmp = set()
                    for ii in range(0, i+1, 1):
                        tmp.add(samples_dict[location][i])
                    if i == 0:
                        i = 1
                    acc_rate = len(tmp) / (i + 1)
                    #acc_rate = len(set(samples_dict[location][1:i+1]))/(i+1) # current acceptance rate

                    if acc_rate > 0.44:
                        sd_prop[location] *= exp(min(0.1, 1 / sqrt(i)))
                    else:
                        sd_prop[location] /= exp(min(0.1, 1 / sqrt(i)))
                        #
                        #sys.exit()
        ## 3) remove burn-in samples
        if self.burn_in > 0:
            for location in samples_dict.keys():
                for i in range(0, self.burn_in, 1):
                    numpy.delete(samples_dict[location],i)
        return samples_dict

    def log_p_of_signal(self, S, sample_dic, i_sample):
        # find all coordinates on that S is conditioned
        log_p = 0
        R = []
        coor = []

        if len(S.sensor.delta_coor) > 0:

            for delta_coor in S.sensor.delta_coor:
                coor.append(S.position + delta_coor)

            for c in coor:
                R.append(sample_dic[self.getKeyVal(sample_dic, c)][i_sample])

        # # get value of integrated Domain
        if not S.sensor.domain_extent.IsZero():
            domain = Domain(S.position, S.sensor.domain_extent, S.angle)
            I = sample_dic[self.getKeyVal(sample_dic, domain)][i_sample]

        # compute log_p
        if S.sensor.domain_extent.IsZero():  # no domain
            log_p = S.sensor.log_p(S.signalData, R)
        if len(S.sensor.delta_coor) == 0:  # no coordiantes
            log_p = S.sensor.log_p(S.signalData, I)

        if not S.sensor.domain_extent.IsZero:
            if len(S.sensor.delta_coor) > 0:
                log_p = S.sensor.log_p(S.signalData, R, I)

        return log_p

    def getKeyVal(self, dic, key):
        # TODO optimize
        logger.debug(shaw(dic,"dic",0))
        for d in dic.keys():
            if d == key:
                return d

        print "No key found < %s >"%key
        logger.error("No key found < %s >"%key)

    def make_sample_point_dict(self, signal_dict, n_samples):
        '''
        :param signal_dict:
        :type signal_dict: Dict{Location, Vector{Signal}
        :param n_samples:
        :type n_samples:
        :return:
        :rtype:
        '''
        Dic1 = {}
        for c in signal_dict:
            n = numpy.zeros(shape=n_samples)
            Dic1[c] = n
        return Dic1

    def make_signal_dict(self, signals):
        Dic = {}
        for S in signals:
            # construct integration domain and add to Dict
            if not S.sensor.domain_extent.IsZero():
                domain = Domain(S.position, S.sensor.domain_extent, S.angle)

                if domain in Dic:
                    Dic[domain] = Set([S, Dic[domain]])
                else:
                    Dic[domain] = [S]

            if len(S.sensor.delta_coor) > 0:
                coors = []
                for delta_coor in S.sensor.delta_coor:
                    coors.append(S.position + delta_coor)

                for i in range(0, len(coors), 1):
                    if coors[i] in Dic:
                        Dic[coors[i]] = Set([S, Dic[coors[i]]])
                    else:
                        Dic[coors[i]] = [S]
        return Dic







