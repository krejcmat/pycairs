__author__ = 'matt'
from math import ceil
import numpy as np
import random
from cairs_gauss import Covarience
import logging
logger = logging.getLogger('cairs')
from cairs_helpers import shaw as sh


class SamplePrediction(object):
    # ---------------------------------
    # # sample from density p(R1*, R2*, ..., Rn* | R1, R2, ..., Rn, I1, ..., Ik)
    # #
    # # loc_pred:    vector with locations (Coor, Domain) for predictions
    # # R_dict_cal:  dictionary containg the samples of calibration
    ## n_samples:   number of samples
    ## prior_mean:     mean function of Prior, f(c::Coor)
    ## prior_cov:      covariance function of Prior, f(c1::Coor, c2::Coor)
    ## optional arguments:
    ## block_size:  number of sample points per block.
    ##
    ## Note:
    ## If 'block_size' is small as the number of prediction
    ## locations, only marginals are sampled, i.e. not all dependecies are
    ## visible in the sample!  Set to 'Inf' to sample a realization of the
    ## rain field (might be much slower).
    ##
    ## Scaling: O(n) with number of prediction points
    ##          O(n^3) with n_calib and blocksize
    def __init__(self, loc_pred, R_dict_cal, n_samples, prior_mean, prior_cov, block_size=200):
        """
        :param loc_pred:
        :type loc_pred:
        :param R_dict_cal:
        :type R_dict_cal:
        :param n_samples:
        :type n_samples:
        :param prior_mean:
        :type prior_mean:
        :param prior_cov:
        :type prior_cov:
        :param block_size:
        :type block_size:
        :return:
        :rtype:
        """
        self.loc_pred=loc_pred
        self.R_dict_cal=R_dict_cal
        self.n_sample=n_samples
        self.prior_mean = prior_mean
        self.prior_cov = prior_cov
        self.block_size = block_size
        self.R_dict_pred = {}

    def sample_prediction(self):
        loc_c = self.R_dict_cal.keys()
        loc_pred_cal = []

        # separate locations that have already been used used for calibration
        loc_pred=[]

        for x in self.loc_pred:
            #logger.debug("x = %s"%x)
            if x in loc_c:
                loc_pred_cal.append(x)
            else:
                loc_pred.append(x)
        #if
        self.loc_pred=loc_pred
        #logger.debug(sh(loc_pred_cal,"loc_pred_cal",0))
        #logger.debug(sh(self.loc_pred,'loc_predXX',prt=0))

        # block_size == Inf ? block_size = size(loc_pred,1) : nothing # only one block #TODO
        #logger.debug("self.block_size = %s"%self.block_size)

        #self.block_size = len(self.loc_pred)  #todo check
        n_bloc = int(ceil(float(len(self.loc_pred)) / float(self.block_size)))
        #logger.debug("n_bloc = %s"%n_bloc)
        mu_c = []
        for loc in loc_c:
            mu_c.append(self.prior_mean.f_mu(loc))
        #logger.debug(sh(loc_c,"loc_c",0))
        #logger.debug(sh(mu_c,"mu_c",0))

        Sigma_cc = self.prior_cov.make_cov(loc_c, loc_c)
        #logger.debug(sh(Sigma_cc,"Sigma_cc",0))
        #logger.debug("Sigma_cc typ = %s"%type(Sigma_cc))
        for i in range(0,n_bloc,1):
            tmp=list()
            #aa=(i * self.block_size +1)
            #bb=min((i+1)*self.block_size,len(self.loc_pred))
            #logger.debug("aa = %s"%aa)
            #logger.debug("bb = %s"%bb)
            for idx in range((i * self.block_size),min((i+1)*self.block_size,len(self.loc_pred)),1):
                tmp.append(self.loc_pred[idx])

            #logger.debug("loc_pred[index] %s ="%tmp )

            sampDict=self.sample_prediction_block(tmp,self.R_dict_cal,self.n_sample,mu_c,Sigma_cc)
            #logger.debug("block %s ="%len(sampDict) )
            self.R_dict_pred.update(sampDict.copy())#merge two dict

        n_calib=self.R_dict_cal.values()
        n_calib= len(n_calib[0])

        if len(loc_pred_cal)>0:
            for loc in loc_pred_cal:
                #rand_index=[i for i in range(1,n_calib,1)]
                rand_index = random.randrange(0,n_calib,self.n_sample)
                self.R_dict_pred[loc]=self.R_dict_cal[loc][rand_index]

        return self.R_dict_pred


    def sample_prediction_block(self,loc_pred,R_dict_cal,n_samples,
                                mu_c,Sigma_cc):
        #logger.debug("sample_prediction_block"+"***"*10)
        loc_c=R_dict_cal.keys()
        logger.debug(sh(loc_c,"loc_c",0))
        #logger.debug(sh(set(loc_c),"loc_c set",0))
        mu_p=[]
        for loc in loc_pred:
            mu_p.append(self.prior_mean.f_mu(loc))
        mu_p = np.array(mu_p)

        #logger.debug(sh(loc_c,"loc_c1",0))
        #logger.debug("mu_p1 = %s"%mu_p)

        #logger.debug(sh(loc_pred,"loc_pred1",0))
        #logger.debug("loc_c = %s"%loc_c)

        Sigma_pp = self.prior_cov.make_cov(loc_pred,loc_pred)

        Sigma_pc = self.prior_cov.make_cov(loc_pred,loc_c)
        logger.debug(sh(Sigma_cc,"Sigma_cc1",0))
        logger.debug(sh(Sigma_pc,"Sigma_pc1",0))

        logger.debug(sh(Sigma_pp,'Sigma_pp1',0))
        #Sigma_cond = Sigma_pp - Sigma_pc * inv(Sigma) * Sigma_pc'
        #sig_inv_cc=np.linalg.inv(Sigma_cc)
        #logger.debug(sh(sig_inv_cc,"sig_inv_cc",0))
        tmpsq= np.dot(Sigma_pc,np.linalg.inv(Sigma_cc))
        tmpsq= np.dot(tmpsq,np.transpose(Sigma_pc))
        #logger.debug(sh(tmpsq,"tmpsq1",0))

        Sigma_cond = Sigma_pp - tmpsq
        #logger.debug(sh(Sigma_cond,"Sigma_cond1",0))
        Sigma_cond_chol = np.linalg.cholesky(Sigma_cond)
        logger.debug(sh(Sigma_cond_chol,"Sigma_cond_chol",0))
        tmp=list()
        for i in R_dict_cal.values():
            tmp.append(i[0])
        n_claib= len(tmp)
        n_pred = len(loc_pred)
        #R_array_pred = n_samples+n_pred
        #handler.setLevel(logging.INFO)
        R_array_pred = np.zeros(shape=(n_samples, n_pred), dtype=np.float64)
        #logger.debug("R_array_pred %s"%R_array_pred)
        for i in range(0,n_samples,1):
            rand_index = random.uniform(1,n_claib)
            R=[]
            for c in R_dict_cal.keys():
                R.append(R_dict_cal[c][rand_index])

            R_mu_c = []
            for a, b in zip(R,mu_c):
                R_mu_c.append(long(a)-long(b))
            R_mu_c=np.array(R_mu_c)

            mtmp=np.dot(Sigma_pc,np.linalg.inv(Sigma_cc))
            mean = mu_p + np.dot(mtmp,R_mu_c)
            #logger.debug(sh(mean,"mean",0))
            R_array_pred[i,:] = mean + np.dot(Sigma_cond_chol, np.random.standard_normal(n_pred))

        #logger.debug("R_array_pred %s"%R_array_pred)
        R_dict_pred = {}
        for i in range(0,len(loc_pred)):
            R_dict_pred[loc_pred[i]] = R_array_pred[:,i]
        #logger.debug("R_dict_pred %s"%R_dict_pred)
        return R_dict_pred


    def is_pos_def(self, x):
        return np.all(np.linalg.eigvals(x) > 0)