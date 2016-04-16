#!/usr/bin/env python

#%Module
#% description: Cairs
#% keyword: interpolation
#% keyword: raster
#%End

import sys
from grass.script import core as grass
import os
import numpy as np
import datetime
import logging
from subprocess import call

from scipy.stats import norm

from cairs_coor import Coor, StaticContext
from cairs_sensor import Sensor
from cairs_gauss import Mean, Covarience
from cairs_signal_handling import MySignals
from cairs_prediction import Prediction
from cairs_helpers import shaw,DataPostProcessor
from cairs_grass_interface import GrassInterface

REF_TIME = None
log_p_MWLLen=6
nn=20
ixa=10

coorRG=Coor(5, 6)
coorMW=Coor(4.2, 2)
angle=0.9

def log_p_gauge(S, R):
    mu = 0.1 + R[0] ** 2
    sigma = 0.005

    return norm.logpdf(np.float(S), loc=np.float(mu), scale=np.float(sigma))


def log_p_MWL(S, I):
    mu = I / log_p_MWLLen
    sigma = 0.005

    # # log of normal density, p(S|R)
    # The probability density function for norm is:
    # norm.pdf(x) = exp(-x**2/2)/sqrt(2*pi)
    return norm.logpdf(np.float(S), loc=np.float(mu), scale=np.float(sigma))
    # return scipy.stats.norm(mu, sigma).logpdf(S)


def timeFromToList(from_date, to_date, date_format='%d.%m.%Y %H:%M:%S', **args):
    start = datetime.datetime.strptime(from_date, date_format)
    end = datetime.datetime.strptime(to_date, date_format)

    delta = datetime.timedelta(**args)
    ret = []
    while start <= end:
        # time = start - REF_TIME
        # tm = timedelta(microseconds=time.microseconds)
        ret.append(start)
        start.strftime(date_format)
        start += delta
    return ret


def initLoger():
    logger = logging.getLogger('cairs')
    # logger.setLevel(logging.INFO)
    logging.basicConfig(level=logging.DEBUG)
    handler = logging.FileHandler('cairs_LOGGER.log', mode='w')
    handler.setLevel(logging.DEBUG)
    #formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info("cairs init logger")
    return logger


def main():
    logger = initLoger()
    logger.propagate = False
    st = StaticContext()
    global REF_TIME
    REF_TIME = st.REF_TIME
    path=""
    #dateNow=datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
    dateNow="1"
    name_rain="%srain_field_%s.csv"%(path,dateNow)

    sensor_name="%ssensor_coor_%s.csv"%(path,dateNow)
    pdf_name="%sut_%s.pdf"%(path,datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))



    sensorGauge = Sensor(log_p_gauge)
    sensorMWL = Sensor(log_p_MWL, domain_extent=Coor(log_p_MWLLen, 0))

    mean_GP = Mean(2)
    cov_GP = Covarience(sigma=10,
                        l_spatial=1.5,
                        l_temporal=60,
                        gamma=1)
    sig = MySignals()
    path1 = os.path.join('data', 'Sensor1.csv')
    path2 = os.path.join('data', 'Sensor2.csv')

    sig.add_signal(filecsv=path1,
                   sensor=sensorGauge,
                   position=coorRG,
                   date_format='%d.%m.%Y %H:%M:%S',
                   delim=',')

    sig.add_signal(filecsv=path2,
                   sensor=sensorMWL,
                   position=coorMW,
                   angle=angle,
                   date_format='%d.%m.%Y %H:%M:%S',
                   delim=',')

    #print str(sig.showSig())
   # logger.info(shaw(sig.signals,"sig",0))

    DataPostProcessor.sensor2csv(sig.signals,sensor_name)

    loc_pred = []
    ii = np.linspace(0, ixa, nn)
    jj = np.linspace(0, ixa, nn)

    timetime = timeFromToList('22.11.2013 13:15:00', '22.11.2013 13:24:00', date_format='%d.%m.%Y %H:%M:%S', seconds=60)
    logger.debug(timetime)
    for time in timetime:
        for x in ii:
            for y in jj:
                loc_pred.append(Coor(x, y, time))

    # show(loc_pred,'locpred')
    logger.info("Prediction - start")
    R_pred = Prediction(loc_pred=loc_pred,  # vector or array with locations for predictions
                        signals=sig,  # vector of signals
                        prior_mean=mean_GP,  # mean function of prior
                        prior_cov=cov_GP,  # covariance function of prior
                        n_sample_calib=2000,  # number of iterations of the Gibbs sampler
                        burn_in=500,  # number of removed samples (and length of adaptation)
                        n_sample_pred=600,  # number of samples for predictions
                        delta=90)  # consider all signals within time 'delta'
    prd=R_pred.predict()
    logger.info("Prediction - end")
    logger.debug("bravo")
    #print prd

    DataPostProcessor.summary2csv(prd,name_rain)


    logger.debug(shaw(prd,"result",0,1))

    run="Rscript %scompute_rain_map.r %s %s %s"%(path,name_rain,sensor_name,pdf_name)
    print run
    grass=GrassInterface()
    grass.summary2grass(prd,sig.signals)

if __name__ == "__main__":
    options, flags = grass.parser()
    main()
