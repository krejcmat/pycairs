import numpy as np
import logging

from cairs_signal import Signal
from cairs_calib import Gibbs
from cairs_predict import SamplePrediction


logger = logging.getLogger('cairs')


class Prediction():
    def __init__(self, loc_pred, signals,
                 prior_mean, prior_cov,
                 n_sample_calib=20000, burn_in=-1,
                 n_sample_pred=5000, delta=5 * 60):
        """
        :param loc_pred:Vector or array of 'Location's for predictions
        :type loc_pred:
        :param signals: Vector or array of 'Signal's
        :type signals:
        :param prior_mean: Mean instance of Prior
        :type prior_mean:
        :param prior_cov: Covariance instance of Prior
        :type prior_cov:
        ## - optional: -
        :param n_sample_calib:number of MCMC samples for calibration
        :type n_sample_calib:
        :param burn_in:number of calibration samples that are discarded
        :type burn_in:
        :param n_sample_pred:number of samples for prediction
        :type n_sample_pred:
        :param delta:maximal time distance of signals to prediction location in Milli seconds
        :type delta:
        :return:
        :rtype:
        """

        self.loc_pred = loc_pred
        self.signals = signals
        self.prior_mean = prior_mean
        self.prior_mean = prior_mean
        self.prior_cov = prior_cov
        self.n_sample_calib = n_sample_calib
        self.burn_in = burn_in
        self.n_sample_pred = n_sample_pred
        self.delta = delta

    def predict(self):
        logger.info('Module - %s, Class - %s ' % (__name__, self.__class__.__name__))

        loc = np.reshape(self.loc_pred, len(self.loc_pred))
        # print loc
        sig = Signal.find_near_signals(loc, self.signals, self.delta)

        if len(sig) == 0:
            print ("No signals are close enough. Increase 'delta'.")
            return

        print("\n--- calibration ---")
        print(" Consider ", len(sig), " signals")
        gibbs = Gibbs(sig, self.prior_mean, self.prior_cov, self.n_sample_calib, self.burn_in)
        Samp_dict_cal = gibbs.runGibbs()
        print "\n ----prdiction-----"
        print "  %s locations"%len(loc)
        R_dict_pred = SamplePrediction(loc, Samp_dict_cal, self.n_sample_pred, self.prior_mean, self.prior_cov, block_size=500)
        res_predict=R_dict_pred.sample_prediction()
        return  res_predict