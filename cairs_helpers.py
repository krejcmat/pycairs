__author__ = 'matt'
import numpy
from cairs_coor import Coor
from statistics import mean,stdev
import csv
import logging
logger = logging.getLogger('cairs')

def prvar(__x):
    import traceback

    print traceback.extract_stack(limit=2)[0][3][6:][:-1], "=", __x



class DataPostProcessor():

    @staticmethod
    def summary2csv(pred_dict,filename = "rain_field.csv", rel = True):
        logger.info("summary2csv - start")
        items=pred_dict.keys()
        n_coor=[]
        for i in items:
            if type(i) is Coor:
                n_coor.append(i)
        predictions=numpy.zeros(shape=(len(n_coor), 7), dtype=np.float64)

        for i in range(0,len(n_coor),1):
            coor = items[i]
            if type(coor)==Coor:
                Rcoor = pred_dict[coor]

            predictions[i,:] = [coor.x,coor.y,coor.time*1000,mean(Rcoor),
                                stdev(Rcoor),np.percentile(Rcoor, 10),
                                np.percentile(Rcoor, 90)]

        #print predictions

        with open(filename, 'wb') as csvfile:
            csw = csv.writer(csvfile, delimiter=',')
            for row in predictions:
                csw.writerow(row)
        logger.info("summary2csv - end")





    @staticmethod
    def sensor2csv(signal,filename = "sensor_coor.csv", rel = True):
        #positions = numpy.zeros(shape=(0, 4))
        logger.info("sensor2csv - start")
        x=[]
        y=[]
        time=[]
        name=[]
        d=1
        for sig in signal:
            if len(sig.sensor.delta_coor)>0:
                for c in sig. sensor.delta_coor:
                    coor = sig.position + c
                    x.append(coor.x)
                    y.append(coor.y)
                    time.append(coor.time)
                    name.append('coor')
            if not sig.sensor.domain_extent.IsZero():
                for b1 in range(0,2):
                    for b2 in range(0,2):
                        for b3 in range(0,2):
                            coor = Coor.rotate(Coor(sig.position.x + sig.sensor.domain_extent.x*b1,
                                                sig.position.y + sig.sensor.domain_extent.y*b2,
                                                sig.position.time + sig.sensor.domain_extent.time*b3),
                                                sig.position, sig.angle)
                            x.append(coor.x)
                            y.append(coor.y)
                            time.append(coor.time)
                            name.append('domain_%s'%d)
                d+=1

        with open(filename, 'wb') as csvfile:
            csw = csv.writer(csvfile, delimiter=',')
            for xx,yy,timet,namen in zip(x,y,time,name):
                csw.writerow([xx, yy, timet*1000, namen])
        logger.info("sensor2csv - end")

def shaw(sig, name='show',prt=True,res=False):
    stra='********* %s ************\n' % name
    if type(sig) is dict:
        for d, a in sig.iteritems():
            if res:
                stra+= '%s => [%s] \n' % (str(d), str(a))
            else:
                stra+= '%s => [%s] \n' % (str(d), str(a[0]))
    else:
        for d in sig:
            stra+= '%s \n' % (str(d))
        stra+= 'len = %s \n'%len(sig)
        if prt:
            print stra
        else:
            return stra

from pprint import pprint

from pylab import *


def arrayToList(arr):
    if type(arr) == type(array([])):
        return arrayToList(arr.tolist())
    elif type(arr) == type([]):
        return [arrayToList(a) for a in arr]
    else:
        return arr


def prettyArray(arr):
    try:
        pprint(arrayToList(arr))
    except:
        return arr