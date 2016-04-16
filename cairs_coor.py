from abc import ABCMeta
from math import sin, cos
import datetime
# from datetime import timedelta
import numpy
import logging
logger = logging.getLogger('cairs')

REF_TIME = datetime.datetime(1984, 10, 20)
# REF_TIME = datetime.datetime.utcfromtimestamp(0)


class StaticContext():
    def __init__(self):
        self.REF_TIME = REF_TIME


class Location:
    __metaclass__ = ABCMeta


class Coor(Location):
    def __init__(self, x, y, time=REF_TIME):
        self.x = x
        self.y = y
        if time == 0:
            time = REF_TIME
        if type(time) is int or type(time) is float or type(time) is numpy.float64:
            self.time = time
        else:
            self.time = ((time - REF_TIME).total_seconds())
            #self.time*=1000

    def __str__(self):
        return "Coor({0},{1},{2})".format(self.x, self.y, self.time)

    def __eq__(self, other):
        if type(other) is not Coor:
            #logger.debug( "warning < %s >, return False (not equal)"%e)
            return False
        if abs(self.x - other.x)> 1e-6:
            return False
        if abs(self.y - other.y)> 1e-6:
            return False
        if abs(self.time - other.time)> 1e-6:
            return False


        return True

    def __add__(self, other):
        # addition of coordinates
        x = self.x + other.x
        y = self.y + other.y

        time = self.time + other.time
        return Coor(x, y, time)

    def IsZero(self):
        if self.x == 0 and self.y == 0:
            if self.time == 0:
                return True
        return False

    def __sub__(self, other):
        # subtraction of coordinates
        x = self.x - other.x
        y = self.y - other.y
        time = self.time - other.time
        return Coor(x, y, time)

    @staticmethod
    def rotate(coor, origin, angle):
        """ rotate coordinates around 'origin' (in spatial dimensions only)
            angle: rotation angle in rad"""
        if angle == 0:
            return coor
        x_trans = coor.x - origin.x
        y_trans = coor.y - origin.y
        return Coor(cos(angle) * x_trans - sin(angle) * y_trans + origin.x,
                    sin(angle) * x_trans + cos(angle) * y_trans + origin.y,
                    coor.time)


class Domain(Location):
    """An imutable class definning type Domain"""

    def __init__(self, position, extend, angle):
        """Constructor"""
        self.position = position
        self.extend = extend  # extend relative to position, before rotation
        self.angle = angle  # angle, rotated around 'position'

    def __eq__(self, other):
        if type(other) is not Domain:
            return False
        if self.position != other.position:
            return False
        if self.extend != other.extend:
            return False
        if abs(self.angle - other.angle)> 1e-6:
            return False
        return True

    def printAddres(self):
        print self

    def __str__(self):
        return 'Domain( %s, %s, %s )' % (str(self.position), str(self.extend), self.angle)

