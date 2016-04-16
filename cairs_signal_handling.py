__author__ = 'matt'
import csv
from datetime import datetime

from cairs_signal import Signal
from cairs_coor import Coor


class MySignals(object):
    def __init__(self):
        self.signals = []

    def add_signal(self, filecsv, sensor, position, date_format, angle=0, delim=','):
        with open(filecsv, 'rb') as csvfile:
            data = csv.reader(csvfile, delimiter=delim)
            for date, val in data:
                time = datetime.strptime(date, date_format)
                self.signals.append(Signal(val, sensor, Coor(position.x, position.y, time), angle))

    def remove_signal(self, sensor=None, coor=None, angle=None):
        if sensor:
            for i, signal in self.signals:
                if signal.sensor is sensor:
                    del self.signals[i]
                    return
        else:
            for i, signal in self.signals:
                if signal.coor == coor and signal.angle == angle:
                    del self.signals[i]
                    return

    def showSig(self):
        for sig in self.signals:
            print str(sig)